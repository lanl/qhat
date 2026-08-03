#!/usr/bin/env python3
"""Run the corrected B2-style benchmark for C-C and N-N.

Place this file in ``analysis/`` on the QHAT ``L-sweep`` branch and run it
from the repository root.

The script reuses the validated B2 benchmark implementation. It evaluates:

* raw JW insertion order;
* signed-coefficient / lexicographic baseline;
* current and reversed JW/fermionic color-group schedules;
* random permutations of JW color groups;
* random permutations of fermionic color groups.

Only ``random_groups`` is sampled; within-group order stays fixed. The summary
uses the corrected state-overlap error

    1 - |<psi_exact | psi_trotter>|.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

import numpy as np

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
    "raw_jw_one_minus_overlap",
    "signed_baseline_one_minus_overlap",
    "jw_mean_one_minus_overlap",
    "jw_median_one_minus_overlap",
    "jw_std_one_minus_overlap",
    "jw_minimum_one_minus_overlap",
    "jw_maximum_one_minus_overlap",
    "jw_schedules_beating_baseline",
    "jw_fraction_beating_baseline",
    "fermionic_mean_one_minus_overlap",
    "fermionic_median_one_minus_overlap",
    "fermionic_std_one_minus_overlap",
    "fermionic_minimum_one_minus_overlap",
    "fermionic_maximum_one_minus_overlap",
    "fermionic_schedules_beating_baseline",
    "fermionic_fraction_beating_baseline",
    "jw_pearson_log_bch_squared_vs_log_overlap_error",
    "fermionic_pearson_log_bch_squared_vs_log_overlap_error",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run B2-style random-group coloring benchmarks for comparable "
            "C-C and N-N active spaces."
        )
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=Path("hamiltonian_generator/comparable_diatomic_library"),
        help="Root containing *.tensors.npz files.",
    )
    parser.add_argument(
        "--molecules",
        nargs="+",
        default=["C-C", "N-N"],
        help="Molecules to select. Default: C-C N-N.",
    )
    parser.add_argument(
        "--basis",
        default="sto-6g",
        help="Basis to select. Default: sto-6g.",
    )
    parser.add_argument(
        "--qubits",
        type=int,
        nargs="+",
        default=[12, 16, 18],
        help="Active-space qubit counts. Default: 12 16 18.",
    )
    parser.add_argument(
        "--bond-length",
        action="append",
        default=[],
        metavar="MOLECULE=VALUE",
        help=(
            "Optional filter, for example --bond-length C-C=1.50 "
            "--bond-length N-N=1.42."
        ),
    )
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument(
        "--time",
        type=float,
        default=1.0,
        dest="evolution_time",
    )
    parser.add_argument("--seed", type=int, default=125)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/comparable_diatomic_coloring_results.csv"),
        help="Detailed coloring result CSV.",
    )
    parser.add_argument(
        "--baseline-output",
        type=Path,
        default=Path("analysis/comparable_diatomic_baseline_results.csv"),
        help="Raw JW and signed-baseline result CSV.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("analysis/comparable_diatomic_corrected_summary.csv"),
        help="Corrected 1-|overlap| summary CSV.",
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
        help="Append coloring results and skip completed schedules.",
    )
    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    if not args.library.exists():
        raise FileNotFoundError(args.library)
    if not args.molecules:
        raise ValueError("--molecules cannot be empty.")
    if not args.qubits or any(value <= 0 for value in args.qubits):
        raise ValueError("--qubits must contain positive integers.")
    if args.samples < 0:
        raise ValueError("--samples cannot be negative.")
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
        result[molecule.strip()] = bond_length.strip()
    return result


def select_cases(
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
            "Generate tensors before running this benchmark."
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
        if requested_length is not None:
            if str(metadata.bond_length) != requested_length:
                continue
        selected[(metadata.molecule, n_qubits)].append(tensor_path)

    missing = [
        (molecule, n_qubits)
        for molecule in molecules
        for n_qubits in qubits
        if not selected[(molecule, n_qubits)]
    ]
    if missing:
        text = ", ".join(
            f"{molecule}/{n_qubits}q" for molecule, n_qubits in missing
        )
        raise FileNotFoundError(f"Missing tensor cases: {text}")

    ambiguous = {
        key: paths for key, paths in selected.items() if len(paths) != 1
    }
    if ambiguous:
        lines = []
        for (molecule, n_qubits), paths in sorted(ambiguous.items()):
            lines.append(f"{molecule}/{n_qubits}q:")
            lines.extend(f"  {path}" for path in paths)
        raise ValueError(
            "More than one tensor matches a case. Add --bond-length filters:\n"
            + "\n".join(lines)
        )

    return [
        selected[(molecule, n_qubits)][0]
        for molecule in molecules
        for n_qubits in sorted(qubits)
    ]


def write_baselines(
    tensor_paths: Sequence[Path],
    args: argparse.Namespace,
) -> None:
    args.baseline_output.parent.mkdir(parents=True, exist_ok=True)
    with args.baseline_output.open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=baseline.FIELDNAMES)
        writer.writeheader()
        for index, tensor_path in enumerate(tensor_paths, start=1):
            print()
            print(
                f"Baseline case [{index}/{len(tensor_paths)}]: {tensor_path}"
            )
            for row in baseline.benchmark_case(tensor_path, args):
                writer.writerow(row)
            stream.flush()


def run_coloring_benchmarks(
    tensor_paths: Sequence[Path],
    args: argparse.Namespace,
) -> None:
    # Restrict the imported B2 robustness script to the experiment used in the
    # corrected comparison: random color-group permutations only.
    robustness.RANDOM_SCHEDULES = ("random_groups",)

    if args.resume:
        (
            completed,
            raw_infidelities,
            current_infidelities,
            raw_bch_norms,
        ) = robustness.load_resume_data(args.output)
    else:
        completed = set()
        raw_infidelities = {}
        current_infidelities = {}
        raw_bch_norms = {}

    args.output.parent.mkdir(parents=True, exist_ok=True)
    append = (
        args.resume
        and args.output.exists()
        and args.output.stat().st_size > 0
    )
    mode = "a" if append else "w"

    with args.output.open(mode, newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=robustness.FIELDNAMES)
        if not append:
            writer.writeheader()
            stream.flush()

        for index, tensor_path in enumerate(tensor_paths, start=1):
            print()
            print(
                f"Coloring case [{index}/{len(tensor_paths)}]: {tensor_path}"
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


def read_success_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as stream:
        return [
            row
            for row in csv.DictReader(stream)
            if row.get("status") == "success"
        ]


def overlap_error(row: dict[str, str]) -> float:
    return 1.0 - float(row["state_overlap_abs"])


def method_statistics(
    rows: Sequence[dict[str, str]],
    signed_baseline: float,
) -> dict[str, float | int | str]:
    values = np.asarray([overlap_error(row) for row in rows], dtype=float)
    if values.size == 0:
        return {
            "mean": "",
            "median": "",
            "std": "",
            "minimum": "",
            "maximum": "",
            "count_beating": "",
            "fraction_beating": "",
        }
    count_beating = int(np.sum(values < signed_baseline))
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "std": float(np.std(values, ddof=0)),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "count_beating": count_beating,
        "fraction_beating": count_beating / int(values.size),
    }


def bch_error_correlation(rows: Sequence[dict[str, str]]) -> float | str:
    x_values: list[float] = []
    y_values: list[float] = []
    for row in rows:
        bch = float(row["bch2_hf_state_norm"])
        error = overlap_error(row)
        if bch > 0.0 and error > 0.0:
            x_values.append(math.log10(bch**2))
            y_values.append(math.log10(error))
    if len(x_values) < 2:
        return ""
    return float(np.corrcoef(x_values, y_values)[0, 1])


def write_corrected_summary(args: argparse.Namespace) -> None:
    baseline_rows = read_success_rows(args.baseline_output)
    coloring_rows = read_success_rows(args.output)

    baseline_by_case: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in baseline_rows:
        baseline_by_case[row["case_id"]][row["ordering"]] = row

    coloring_by_case: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in coloring_rows:
        if row["schedule"] != "random_groups":
            continue
        coloring_by_case[row["case_id"]].append(row)

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", newline="", encoding="utf-8") as stream:
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
            if "signed_coefficient_lexicographic" not in baselines:
                print(f"Warning: signed baseline missing for {case_id}")
                continue
            if "jw_raw" not in baselines:
                print(f"Warning: raw JW baseline missing for {case_id}")
                continue

            signed_error = overlap_error(
                baselines["signed_coefficient_lexicographic"]
            )
            raw_error = overlap_error(baselines["jw_raw"])
            jw_rows = [row for row in rows if row["graph_level"] == "jw"]
            fermionic_rows = [
                row for row in rows if row["graph_level"] == "fermionic"
            ]
            jw_stats = method_statistics(jw_rows, signed_error)
            fermionic_stats = method_statistics(
                fermionic_rows, signed_error
            )
            first = rows[0]

            writer.writerow(
                {
                    "case_id": case_id,
                    "molecule": first["molecule"],
                    "bond_length": first["bond_length"],
                    "basis": first["basis"],
                    "active_occupied": first["active_occupied"],
                    "active_vacant": first["active_vacant"],
                    "n_qubits": first["n_qubits"],
                    "number_of_fermionic_terms": first[
                        "number_of_fermionic_terms"
                    ],
                    "number_of_pauli_terms": first["number_of_pauli_terms"],
                    "jw_number_of_colors": jw_rows[0]["number_of_colors"],
                    "fermionic_number_of_colors": fermionic_rows[0][
                        "number_of_colors"
                    ],
                    "raw_jw_one_minus_overlap": raw_error,
                    "signed_baseline_one_minus_overlap": signed_error,
                    "jw_mean_one_minus_overlap": jw_stats["mean"],
                    "jw_median_one_minus_overlap": jw_stats["median"],
                    "jw_std_one_minus_overlap": jw_stats["std"],
                    "jw_minimum_one_minus_overlap": jw_stats["minimum"],
                    "jw_maximum_one_minus_overlap": jw_stats["maximum"],
                    "jw_schedules_beating_baseline": jw_stats[
                        "count_beating"
                    ],
                    "jw_fraction_beating_baseline": jw_stats[
                        "fraction_beating"
                    ],
                    "fermionic_mean_one_minus_overlap": fermionic_stats[
                        "mean"
                    ],
                    "fermionic_median_one_minus_overlap": fermionic_stats[
                        "median"
                    ],
                    "fermionic_std_one_minus_overlap": fermionic_stats["std"],
                    "fermionic_minimum_one_minus_overlap": fermionic_stats[
                        "minimum"
                    ],
                    "fermionic_maximum_one_minus_overlap": fermionic_stats[
                        "maximum"
                    ],
                    "fermionic_schedules_beating_baseline": fermionic_stats[
                        "count_beating"
                    ],
                    "fermionic_fraction_beating_baseline": fermionic_stats[
                        "fraction_beating"
                    ],
                    "jw_pearson_log_bch_squared_vs_log_overlap_error": (
                        bch_error_correlation(jw_rows)
                    ),
                    "fermionic_pearson_log_bch_squared_vs_log_overlap_error": (
                        bch_error_correlation(fermionic_rows)
                    ),
                }
            )


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    bond_lengths = parse_bond_lengths(args.bond_length)

    unknown = set(bond_lengths).difference(args.molecules)
    if unknown:
        raise ValueError(
            "Bond-length filter supplied for an unselected molecule: "
            + ", ".join(sorted(unknown))
        )

    warm_up_numba()
    tensor_paths = select_cases(
        library=args.library,
        molecules=args.molecules,
        basis_name=args.basis,
        qubits=args.qubits,
        bond_lengths=bond_lengths,
    )

    print(f"Selected {len(tensor_paths)} tensor cases")
    print(f"Molecules: {args.molecules}")
    print(f"Basis: {args.basis}")
    print(f"Qubits: {sorted(args.qubits)}")
    print(f"Random-group samples per graph: {args.samples}")
    print("Primary summary metric: 1 - |<psi_exact|psi_trotter>|")

    write_baselines(tensor_paths, args)
    run_coloring_benchmarks(tensor_paths, args)
    write_corrected_summary(args)

    print()
    print("Complete")
    print(f"Baseline CSV: {args.baseline_output}")
    print(f"Detailed CSV: {args.output}")
    print(f"Summary CSV:  {args.summary}")


if __name__ == "__main__":
    main()
