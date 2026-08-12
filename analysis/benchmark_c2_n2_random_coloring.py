#!/usr/bin/env python3
"""Run the complete C2/N2 random-coloring case study for STO-6G and HGBS-5.

For every selected molecular Hamiltonian this script compares:

1. the deterministic ascending signed-coefficient JW baseline
   (``signed_coefficient_lexicographic``);
2. 20 random permutations of the direct JW color groups; and
3. 20 random permutations of the fermionic color groups, followed by the
   induced final JW Pauli ordering.

The default Trotter experiment is first order, total time t=1, 100 steps.
The reported error is

    1 - |<psi_exact | psi_trotter>|.

The random experiment changes color-group order only.  Within-group order is
kept fixed, matching the corrected B2 robustness study.

This script intentionally reuses the validated B2 matrix-free benchmark code
rather than reimplementing the Trotter evolution, graph construction, or BCH
evaluator.  Place it in ``analysis/`` on the QHAT L-sweep branch.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from qhat.analysis import benchmark_b2_coloring_robustness as robustness
    from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
    from qhat.analysis.benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        discover_tensor_paths,
        load_interaction_operator,
        parse_case_metadata,
    )
    from qhat.analysis.benchmark_b2_active_spaces_matrix_free import warm_up_numba
except ImportError:
    import benchmark_b2_coloring_robustness as robustness
    import benchmark_b2_signed_coefficient_baseline as baseline
    from benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        discover_tensor_paths,
        load_interaction_operator,
        parse_case_metadata,
    )
    from benchmark_b2_active_spaces_matrix_free import warm_up_numba


DEFAULT_MOLECULES = ("C-C", "N-N")
DEFAULT_BASES = ("sto-6g", "hgbs-5")
DEFAULT_QUBITS = tuple(range(4, 21, 2))
SIGNED_BASELINE = "signed_coefficient_lexicographic"
RANDOM_SCHEDULE = "random_groups"

SUMMARY_FIELDS = [
    "case_id",
    "molecule",
    "bond_length",
    "basis",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "number_of_fermionic_terms",
    "number_of_pauli_terms",
    "jw_number_of_colors",
    "fermionic_number_of_colors",
    "signed_baseline_one_minus_overlap",
    "jw_random_samples",
    "jw_mean_one_minus_overlap",
    "jw_median_one_minus_overlap",
    "jw_std_one_minus_overlap",
    "jw_minimum_one_minus_overlap",
    "jw_maximum_one_minus_overlap",
    "jw_mean_ratio_to_baseline",
    "jw_schedules_beating_baseline",
    "jw_fraction_beating_baseline",
    "fermionic_random_samples",
    "fermionic_mean_one_minus_overlap",
    "fermionic_median_one_minus_overlap",
    "fermionic_std_one_minus_overlap",
    "fermionic_minimum_one_minus_overlap",
    "fermionic_maximum_one_minus_overlap",
    "fermionic_mean_ratio_to_baseline",
    "fermionic_schedules_beating_baseline",
    "fermionic_fraction_beating_baseline",
    "jw_mean_over_fermionic_mean",
    "jw_pearson_log_bch_squared_vs_log_overlap_error",
    "fermionic_pearson_log_bch_squared_vs_log_overlap_error",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--library",
        type=Path,
        default=Path("hamiltonian_generator/comparable_diatomic_library"),
        help="Root containing the C-C/N-N *.tensors.npz files.",
    )
    parser.add_argument(
        "--molecules",
        nargs="+",
        default=list(DEFAULT_MOLECULES),
        help="Molecules to benchmark. Default: C-C N-N.",
    )
    parser.add_argument(
        "--bases",
        nargs="+",
        default=list(DEFAULT_BASES),
        help="Basis sets to benchmark. Default: sto-6g hgbs-5.",
    )
    parser.add_argument(
        "--qubits",
        type=int,
        nargs="+",
        default=list(DEFAULT_QUBITS),
        help="Active-space qubit counts. Default: 4 6 ... 20.",
    )
    parser.add_argument(
        "--bond-length",
        action="append",
        default=[],
        metavar="MOLECULE=VALUE",
        help=(
            "Optional exact bond-length filter, e.g. "
            "--bond-length C-C=1.50 --bond-length N-N=1.42. "
            "If omitted, each molecule/qubit/basis must have exactly one tensor."
        ),
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=20,
        help="Random color-group permutations per graph level. Default: 20.",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=100,
        help="First-order Trotter steps. Default: 100.",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=1.0,
        dest="evolution_time",
        help="Total evolution time. Default: 1.0.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=125,
        help="Base random seed. Default: 125.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis/c2_n2_case_study"),
        help="Directory for detailed CSVs, summaries, and figures.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
    )
    parser.add_argument(
        "--parallel-threshold",
        type=int,
        default=2**16,
    )
    parser.add_argument("--no-spin-sector", action="store_true")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Append coloring rows and skip schedules already completed.",
    )
    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    if not args.library.exists():
        raise FileNotFoundError(args.library)
    if not args.molecules:
        raise ValueError("--molecules cannot be empty.")
    if not args.bases:
        raise ValueError("--bases cannot be empty.")
    if not args.qubits or any(q <= 0 or q % 2 for q in args.qubits):
        raise ValueError("--qubits must contain positive even integers.")
    if args.samples <= 0:
        raise ValueError("--samples must be positive.")
    if args.steps <= 0:
        raise ValueError("--steps must be positive.")
    if args.evolution_time <= 0.0:
        raise ValueError("--time must be positive.")
    if args.tolerance <= 0.0:
        raise ValueError("--tolerance must be positive.")
    if args.parallel_threshold <= 0:
        raise ValueError("--parallel-threshold must be positive.")


def parse_bond_lengths(values: Sequence[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(
                "--bond-length must use MOLECULE=VALUE, such as C-C=1.50."
            )
        molecule, bond_length = value.split("=", 1)
        molecule = molecule.strip()
        number = float(bond_length)
        if number <= 0.0:
            raise ValueError("Bond lengths must be positive.")
        # Tensor directories/config stubs use two decimal places in this study.
        result[molecule] = f"{number:.2f}"
    return result


def select_cases(
    *,
    library: Path,
    molecules: Sequence[str],
    basis_name: str,
    qubits: Sequence[int],
    bond_lengths: dict[str, str],
) -> list[Path]:
    requested_molecules = set(molecules)
    requested_qubits = set(qubits)
    selected: dict[tuple[str, int], list[Path]] = defaultdict(list)

    tensor_paths = discover_tensor_paths(library, case_pattern=None, limit=None)
    if not tensor_paths:
        raise FileNotFoundError(
            f"No *.tensors.npz files found under {library}. "
            "Run build_c2_n2_active_spaces.py first."
        )

    for tensor_path in tensor_paths:
        try:
            interaction, n_qubits = load_interaction_operator(tensor_path)
            del interaction
            metadata = parse_case_metadata(tensor_path, n_qubits)
        except Exception:
            continue

        if metadata.molecule not in requested_molecules:
            continue
        if metadata.basis.lower() != basis_name.lower():
            continue
        if n_qubits not in requested_qubits:
            continue
        requested_length = bond_lengths.get(metadata.molecule)
        if requested_length is not None and str(metadata.bond_length) != requested_length:
            continue
        selected[(metadata.molecule, n_qubits)].append(tensor_path)

    missing = [
        (molecule, n_qubits)
        for molecule in molecules
        for n_qubits in sorted(requested_qubits)
        if not selected[(molecule, n_qubits)]
    ]
    if missing:
        text = ", ".join(f"{molecule}/{n_qubits}q" for molecule, n_qubits in missing)
        raise FileNotFoundError(
            f"Missing tensor cases for basis {basis_name}: {text}"
        )

    ambiguous = {
        key: paths for key, paths in selected.items() if len(paths) != 1
    }
    if ambiguous:
        lines: list[str] = []
        for (molecule, n_qubits), paths in sorted(ambiguous.items()):
            lines.append(f"{molecule}/{basis_name}/{n_qubits}q:")
            lines.extend(f"  {path}" for path in paths)
        raise ValueError(
            "More than one tensor matches a case. Supply --bond-length filters:\n"
            + "\n".join(lines)
        )

    return [
        selected[(molecule, n_qubits)][0]
        for molecule in molecules
        for n_qubits in sorted(requested_qubits)
    ]


def read_success_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as stream:
        return [
            row
            for row in csv.DictReader(stream)
            if row.get("status") == "success"
        ]


def write_baselines(
    tensor_paths: Sequence[Path],
    args: argparse.Namespace,
    baseline_output: Path,
) -> None:
    """Run the validated deterministic benchmark; summary uses signed baseline."""
    baseline_output.parent.mkdir(parents=True, exist_ok=True)
    completed_cases: set[str] = set()

    append = args.resume and baseline_output.exists() and baseline_output.stat().st_size > 0
    if append:
        for row in read_success_rows(baseline_output):
            if row.get("ordering") == SIGNED_BASELINE:
                completed_cases.add(row["case_id"])

    mode = "a" if append else "w"
    with baseline_output.open(mode, newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=baseline.FIELDNAMES)
        if not append:
            writer.writeheader()
            stream.flush()

        for index, tensor_path in enumerate(tensor_paths, start=1):
            interaction, n_qubits = load_interaction_operator(tensor_path)
            del interaction
            metadata = parse_case_metadata(tensor_path, n_qubits)
            if metadata.case_id in completed_cases:
                print(
                    f"Baseline [{index}/{len(tensor_paths)}] SKIP {metadata.case_id}"
                )
                continue

            print(
                f"Baseline [{index}/{len(tensor_paths)}] RUN  {metadata.case_id}"
            )
            # benchmark_case currently evaluates the full deterministic ordering
            # set.  Keeping that implementation unchanged preserves the validated
            # state evolution.  The corrected summary below uses only the
            # signed_coefficient_lexicographic row as the baseline.
            for row in baseline.benchmark_case(tensor_path, args):
                writer.writerow(row)
            stream.flush()


def run_coloring_benchmarks(
    tensor_paths: Sequence[Path],
    args: argparse.Namespace,
    coloring_output: Path,
) -> None:
    # Exactly the requested experiment: shuffle color-group order only.
    robustness.RANDOM_SCHEDULES = (RANDOM_SCHEDULE,)

    if args.resume:
        (
            completed,
            raw_infidelities,
            current_infidelities,
            raw_bch_norms,
        ) = robustness.load_resume_data(coloring_output)
    else:
        completed = set()
        raw_infidelities = {}
        current_infidelities = {}
        raw_bch_norms = {}

    coloring_output.parent.mkdir(parents=True, exist_ok=True)
    append = (
        args.resume
        and coloring_output.exists()
        and coloring_output.stat().st_size > 0
    )
    mode = "a" if append else "w"

    with coloring_output.open(mode, newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=robustness.FIELDNAMES)
        if not append:
            writer.writeheader()
            stream.flush()

        for index, tensor_path in enumerate(tensor_paths, start=1):
            print()
            print(
                f"Coloring [{index}/{len(tensor_paths)}]: {tensor_path}"
            )
            robustness.benchmark_case(
                tensor_path=tensor_path,
                args=args,
                completed=completed,
                raw_infidelities=raw_infidelities,
                current_infidelities=current_infidelities,
                raw_bch_norms=raw_bch_norms,
                writer=writer,
                output_stream=stream,
            )


def overlap_error(row: dict[str, str]) -> float:
    return 1.0 - float(row["state_overlap_abs"])


def method_statistics(
    rows: Sequence[dict[str, str]],
    signed_baseline: float,
) -> dict[str, float | int]:
    values = np.asarray([overlap_error(row) for row in rows], dtype=float)
    if values.size == 0:
        raise ValueError("No random-group rows were available for a method.")
    count_beating = int(np.sum(values < signed_baseline))
    mean_value = float(np.mean(values))
    return {
        "count": int(values.size),
        "mean": mean_value,
        "median": float(np.median(values)),
        "std": float(np.std(values, ddof=0)),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "mean_ratio_to_baseline": (
            mean_value / signed_baseline if signed_baseline > 0.0 else math.nan
        ),
        "count_beating": count_beating,
        "fraction_beating": count_beating / int(values.size),
    }


def bch_error_correlation(rows: Sequence[dict[str, str]]) -> float:
    x_values: list[float] = []
    y_values: list[float] = []
    for row in rows:
        bch = float(row["bch2_hf_state_norm"])
        error = overlap_error(row)
        if bch > 0.0 and error > 0.0:
            x_values.append(math.log10(bch**2))
            y_values.append(math.log10(error))
    if len(x_values) < 2:
        return math.nan
    return float(np.corrcoef(x_values, y_values)[0, 1])


def write_summary(
    *,
    args: argparse.Namespace,
    basis_name: str,
    baseline_output: Path,
    coloring_output: Path,
    summary_output: Path,
) -> None:
    baseline_rows = read_success_rows(baseline_output)
    coloring_rows = read_success_rows(coloring_output)

    baseline_by_case: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in baseline_rows:
        baseline_by_case[row["case_id"]][row["ordering"]] = row

    coloring_by_case: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in coloring_rows:
        if row.get("schedule") != RANDOM_SCHEDULE:
            continue
        if row.get("graph_level") not in {"jw", "fermionic"}:
            continue
        coloring_by_case[row["case_id"]].append(row)

    summary_output.parent.mkdir(parents=True, exist_ok=True)
    with summary_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()

        for case_id, rows in sorted(
            coloring_by_case.items(),
            key=lambda item: (
                item[1][0]["molecule"],
                int(item[1][0]["n_qubits"]),
            ),
        ):
            baselines = baseline_by_case.get(case_id, {})
            if SIGNED_BASELINE not in baselines:
                raise ValueError(
                    f"Signed ascending baseline missing for {case_id}."
                )

            signed_error = overlap_error(baselines[SIGNED_BASELINE])
            jw_rows = [row for row in rows if row["graph_level"] == "jw"]
            fermionic_rows = [
                row for row in rows if row["graph_level"] == "fermionic"
            ]

            jw_stats = method_statistics(jw_rows, signed_error)
            fermionic_stats = method_statistics(fermionic_rows, signed_error)

            if int(jw_stats["count"]) != args.samples:
                raise ValueError(
                    f"{case_id}: expected {args.samples} JW random-group rows, "
                    f"found {jw_stats['count']}."
                )
            if int(fermionic_stats["count"]) != args.samples:
                raise ValueError(
                    f"{case_id}: expected {args.samples} fermionic random-group rows, "
                    f"found {fermionic_stats['count']}."
                )

            first = rows[0]
            jw_mean = float(jw_stats["mean"])
            fermionic_mean = float(fermionic_stats["mean"])

            writer.writerow(
                {
                    "case_id": case_id,
                    "molecule": first["molecule"],
                    "bond_length": first["bond_length"],
                    "basis": basis_name,
                    "active_occupied": first["active_occupied"],
                    "active_vacant": first["active_vacant"],
                    "n_qubits": first["n_qubits"],
                    "number_of_fermionic_terms": first["number_of_fermionic_terms"],
                    "number_of_pauli_terms": first["number_of_pauli_terms"],
                    "jw_number_of_colors": jw_rows[0]["number_of_colors"],
                    "fermionic_number_of_colors": fermionic_rows[0]["number_of_colors"],
                    "signed_baseline_one_minus_overlap": signed_error,
                    "jw_random_samples": jw_stats["count"],
                    "jw_mean_one_minus_overlap": jw_mean,
                    "jw_median_one_minus_overlap": jw_stats["median"],
                    "jw_std_one_minus_overlap": jw_stats["std"],
                    "jw_minimum_one_minus_overlap": jw_stats["minimum"],
                    "jw_maximum_one_minus_overlap": jw_stats["maximum"],
                    "jw_mean_ratio_to_baseline": jw_stats["mean_ratio_to_baseline"],
                    "jw_schedules_beating_baseline": jw_stats["count_beating"],
                    "jw_fraction_beating_baseline": jw_stats["fraction_beating"],
                    "fermionic_random_samples": fermionic_stats["count"],
                    "fermionic_mean_one_minus_overlap": fermionic_mean,
                    "fermionic_median_one_minus_overlap": fermionic_stats["median"],
                    "fermionic_std_one_minus_overlap": fermionic_stats["std"],
                    "fermionic_minimum_one_minus_overlap": fermionic_stats["minimum"],
                    "fermionic_maximum_one_minus_overlap": fermionic_stats["maximum"],
                    "fermionic_mean_ratio_to_baseline": fermionic_stats[
                        "mean_ratio_to_baseline"
                    ],
                    "fermionic_schedules_beating_baseline": fermionic_stats[
                        "count_beating"
                    ],
                    "fermionic_fraction_beating_baseline": fermionic_stats[
                        "fraction_beating"
                    ],
                    "jw_mean_over_fermionic_mean": (
                        jw_mean / fermionic_mean if fermionic_mean > 0.0 else math.nan
                    ),
                    "jw_pearson_log_bch_squared_vs_log_overlap_error": (
                        bch_error_correlation(jw_rows)
                    ),
                    "fermionic_pearson_log_bch_squared_vs_log_overlap_error": (
                        bch_error_correlation(fermionic_rows)
                    ),
                }
            )


def combine_summaries(summary_paths: Sequence[Path], output: Path) -> pd.DataFrame:
    frames = [pd.read_csv(path) for path in summary_paths]
    combined = pd.concat(frames, ignore_index=True)
    combined = combined.sort_values(
        ["basis", "molecule", "n_qubits"],
        kind="stable",
    ).reset_index(drop=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output, index=False)
    return combined


def _safe_positive(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return np.where(arr > 0.0, arr, np.nan)


def plot_errors(combined: pd.DataFrame, output: Path) -> None:
    molecules = list(dict.fromkeys(combined["molecule"].astype(str)))
    bases = list(dict.fromkeys(combined["basis"].astype(str)))

    fig, axes = plt.subplots(
        len(molecules),
        len(bases),
        figsize=(6.2 * len(bases), 4.6 * len(molecules)),
        squeeze=False,
        sharex=False,
        sharey=False,
    )

    for row_index, molecule in enumerate(molecules):
        for col_index, basis_name in enumerate(bases):
            ax = axes[row_index][col_index]
            group = combined.loc[
                combined["molecule"].eq(molecule)
                & combined["basis"].eq(basis_name)
            ].sort_values("n_qubits")

            if group.empty:
                ax.set_visible(False)
                continue

            x = group["n_qubits"].to_numpy(dtype=int)
            baseline_values = _safe_positive(
                group["signed_baseline_one_minus_overlap"]
            )
            jw_mean = _safe_positive(group["jw_mean_one_minus_overlap"])
            jw_std = group["jw_std_one_minus_overlap"].to_numpy(dtype=float)
            ferm_mean = _safe_positive(
                group["fermionic_mean_one_minus_overlap"]
            )
            ferm_std = group["fermionic_std_one_minus_overlap"].to_numpy(dtype=float)

            ax.plot(x, baseline_values, marker="D", label="JW ascending baseline")
            jw_yerr = np.vstack([np.minimum(jw_std, 0.95 * jw_mean), jw_std])
            ferm_yerr = np.vstack([
                np.minimum(ferm_std, 0.95 * ferm_mean),
                ferm_std,
            ])
            ax.errorbar(
                x,
                jw_mean,
                yerr=jw_yerr,
                marker="o",
                capsize=3,
                label="JW coloring mean +/- std",
            )
            ax.errorbar(
                x,
                ferm_mean,
                yerr=ferm_yerr,
                marker="s",
                capsize=3,
                label="Fermionic coloring mean +/- std",
            )
            ax.set_yscale("log")
            ax.set_xlabel("Active-space qubits")
            ax.set_ylabel(r"$1-|\langle\psi_{exact}|\psi_{Trotter}\rangle|$")
            ax.set_title(f"{molecule.replace('-', '')} / {basis_name.upper()}")
            ax.grid(True, which="both", alpha=0.25)
            ax.legend(fontsize=8)

    fig.suptitle(
        "C2/N2 random coloring compared with the ascending signed-coefficient baseline",
        y=1.01,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_ratios_to_baseline(combined: pd.DataFrame, output: Path) -> None:
    molecules = list(dict.fromkeys(combined["molecule"].astype(str)))
    bases = list(dict.fromkeys(combined["basis"].astype(str)))

    fig, axes = plt.subplots(
        len(molecules),
        len(bases),
        figsize=(6.2 * len(bases), 4.4 * len(molecules)),
        squeeze=False,
    )

    for row_index, molecule in enumerate(molecules):
        for col_index, basis_name in enumerate(bases):
            ax = axes[row_index][col_index]
            group = combined.loc[
                combined["molecule"].eq(molecule)
                & combined["basis"].eq(basis_name)
            ].sort_values("n_qubits")
            if group.empty:
                ax.set_visible(False)
                continue

            x = group["n_qubits"].to_numpy(dtype=int)
            jw_ratio = _safe_positive(group["jw_mean_ratio_to_baseline"])
            ferm_ratio = _safe_positive(group["fermionic_mean_ratio_to_baseline"])

            ax.plot(x, jw_ratio, marker="o", label="JW coloring mean / baseline")
            ax.plot(
                x,
                ferm_ratio,
                marker="s",
                label="Fermionic coloring mean / baseline",
            )
            ax.axhline(1.0, linestyle="--", linewidth=1.2)
            ax.set_yscale("log")
            ax.set_xlabel("Active-space qubits")
            ax.set_ylabel("Mean random-coloring error / ascending baseline error")
            ax.set_title(f"{molecule.replace('-', '')} / {basis_name.upper()}")
            ax.grid(True, which="both", alpha=0.25)
            ax.legend(fontsize=8)

    fig.suptitle(
        "Normalized comparison to the ascending signed-coefficient baseline\n"
        "Below 1 means the random-coloring mean is better than the baseline",
        y=1.03,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def print_compact_summary(combined: pd.DataFrame) -> None:
    columns = [
        "molecule",
        "basis",
        "n_qubits",
        "signed_baseline_one_minus_overlap",
        "jw_mean_one_minus_overlap",
        "fermionic_mean_one_minus_overlap",
        "jw_schedules_beating_baseline",
        "fermionic_schedules_beating_baseline",
        "jw_mean_over_fermionic_mean",
    ]
    print()
    print("Corrected C2/N2 summary")
    print(combined[columns].to_string(index=False))


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.qubits = sorted(set(args.qubits))
    args.bases = [basis.lower() for basis in args.bases]
    bond_lengths = parse_bond_lengths(args.bond_length)

    unknown_bond_molecules = set(bond_lengths).difference(args.molecules)
    if unknown_bond_molecules:
        raise ValueError(
            "Bond-length filter supplied for an unselected molecule: "
            + ", ".join(sorted(unknown_bond_molecules))
        )

    print("C2/N2 coloring case study")
    print(f"Library: {args.library}")
    print(f"Molecules: {args.molecules}")
    print(f"Bases: {args.bases}")
    print(f"Qubits: {args.qubits}")
    print(f"Random group permutations per method/case: {args.samples}")
    print(f"First-order Trotter: t={args.evolution_time}, steps={args.steps}")
    print("Primary metric: 1 - |<psi_exact|psi_trotter>|")
    print(f"Baseline: {SIGNED_BASELINE}")

    warm_up_numba()
    summary_paths: list[Path] = []

    for basis_name in args.bases:
        print()
        print("#" * 96)
        print(f"BASIS: {basis_name}")
        print("#" * 96)

        tensor_paths = select_cases(
            library=args.library,
            molecules=args.molecules,
            basis_name=basis_name,
            qubits=args.qubits,
            bond_lengths=bond_lengths,
        )
        print(f"Selected {len(tensor_paths)} tensor cases for {basis_name}.")

        safe_basis = basis_name.replace("/", "_")
        baseline_output = args.output_dir / f"c2_n2_{safe_basis}_baseline_results.csv"
        coloring_output = args.output_dir / f"c2_n2_{safe_basis}_coloring_results.csv"
        summary_output = args.output_dir / f"c2_n2_{safe_basis}_summary.csv"

        # Some imported helper functions inspect args.output.  Point it at the
        # active detailed file before each phase.
        args.output = baseline_output
        write_baselines(tensor_paths, args, baseline_output)

        args.output = coloring_output
        run_coloring_benchmarks(tensor_paths, args, coloring_output)

        write_summary(
            args=args,
            basis_name=basis_name,
            baseline_output=baseline_output,
            coloring_output=coloring_output,
            summary_output=summary_output,
        )
        summary_paths.append(summary_output)

    combined_path = args.output_dir / "c2_n2_all_summary.csv"
    combined = combine_summaries(summary_paths, combined_path)

    error_figure = args.output_dir / "c2_n2_random_coloring_vs_baseline.png"
    ratio_figure = args.output_dir / "c2_n2_random_coloring_ratio_to_baseline.png"
    plot_errors(combined, error_figure)
    plot_ratios_to_baseline(combined, ratio_figure)
    print_compact_summary(combined)

    print()
    print("Complete")
    print(f"Combined summary: {combined_path}")
    print(f"Error figure:     {error_figure}")
    print(f"Ratio figure:     {ratio_figure}")


if __name__ == "__main__":
    main()
