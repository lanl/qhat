#!/usr/bin/env python3
"""Build and optionally execute a BCH-cancellation validation panel.

The panel is held out from the eight cases used to formulate the current BCH
cancellation hypothesis.  Cases are paired at the same active-space size, with
one case where signed fermionic ordering beats the strongest deterministic JW
ordering and one negative control where it does not.

By default this script only writes an auditable manifest and prints the exact
benchmark commands.  Pass ``--execute`` under the repository's Python 3.11
environment to run the existing structure-ablation benchmark.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PERFORMANCE = Path(
    "analysis/fermionic_aware_performance/fermionic_aware_case_performance.csv"
)
DEFAULT_OUTDIR = Path("analysis/cancellation_hypothesis_validation")

DISCOVERY_CASE_IDS = frozenset(
    {
        "F-F_1.28_hgbs-5_as-016-002",
        "Li-Li_2.66_hgbs-5_as-006-012",
        "Be-Be_2.04_hgbs-5_as-008-012",
        "B-B_1.70_hgbs-5_as-008-010",
        "H2O_s-1.00_hgbs-5_as-008-010",
        "NH3_s-1.00_hgbs-5_as-008-010",
        "BeH2_s-1.00_hgbs-5_as-006-012",
        "O-O_1.26_hgbs-5_as-008-010",
    }
)


@dataclass(frozen=True)
class PanelCase:
    case_id: str
    matched_pair: str
    expected_outcome: str
    panel_tier: str


# The four-case pilot is intentionally affordable enough for an interactive
# validation.  The full panel adds larger, more discriminating held-out cases.
PANEL_CASES = (
    PanelCase(
        "Li-Li_2.66_hgbs-5_as-006-006",
        "12q_pair_a",
        "favorable",
        "pilot",
    ),
    PanelCase(
        "H2O_s-1.00_hgbs-5_as-004-008",
        "12q_pair_a",
        "negative_control",
        "pilot",
    ),
    PanelCase(
        "Be-H_1.47_hgbs-5_as-005-007",
        "12q_heteronuclear",
        "favorable",
        "pilot",
    ),
    PanelCase(
        "Li-H_1.81_hgbs-5_as-004-008",
        "12q_heteronuclear",
        "negative_control",
        "pilot",
    ),
    PanelCase(
        "Be-Be_2.04_hgbs-5_as-008-006",
        "14q_diatomic_a",
        "favorable",
        "full",
    ),
    PanelCase(
        "B-B_1.70_hgbs-5_as-006-008",
        "14q_diatomic_a",
        "negative_control",
        "full",
    ),
    PanelCase(
        "H2O_s-1.00_hgbs-5_as-006-008",
        "14q_mixed",
        "favorable",
        "full",
    ),
    PanelCase(
        "F-F_1.28_hgbs-5_as-012-002",
        "14q_mixed",
        "negative_control",
        "full",
    ),
    PanelCase(
        "F-F_1.28_hgbs-5_as-014-002",
        "16q_diatomic",
        "favorable",
        "full",
    ),
    PanelCase(
        "B-B_1.70_hgbs-5_as-006-010",
        "16q_diatomic",
        "negative_control",
        "full",
    ),
    PanelCase(
        "NH3_s-1.00_hgbs-5_as-006-010",
        "16q_polyatomic",
        "favorable",
        "full",
    ),
    PanelCase(
        "Li-H_1.81_hgbs-5_as-004-012",
        "16q_polyatomic",
        "negative_control",
        "full",
    ),
    PanelCase(
        "CH4_s-1.00_hgbs-5_as-008-010",
        "18q_polyatomic",
        "favorable",
        "full",
    ),
    PanelCase(
        "Be-H_1.47_hgbs-5_as-005-013",
        "18q_polyatomic",
        "negative_control",
        "full",
    ),
    PanelCase(
        "B-H_1.29_hgbs-5_as-006-012",
        "18q_heteronuclear",
        "favorable",
        "full",
    ),
    PanelCase(
        "Li-H_1.81_hgbs-5_as-004-014",
        "18q_heteronuclear",
        "negative_control",
        "full",
    ),
    PanelCase(
        "H2O_s-1.00_hgbs-5_as-008-012",
        "20q_polyatomic",
        "favorable",
        "full",
    ),
    PanelCase(
        "B-H_1.29_hgbs-5_as-006-014",
        "20q_polyatomic",
        "negative_control",
        "full",
    ),
    PanelCase(
        "B-B_1.70_hgbs-5_as-010-010",
        "20q_diatomic",
        "favorable",
        "full",
    ),
    PanelCase(
        "Li-H_1.81_hgbs-5_as-004-016",
        "20q_diatomic",
        "negative_control",
        "full",
    ),
)

DEFAULT_SCHEDULES = (
    "fermionic_signed_reference",
    "fermionic_magnitude_reference",
    "jw_signed_baseline",
    "jw_magnitude_baseline",
    "signed_parent_descendants_round_robin",
    "signed_parent_blocks_randomized",
    "signed_parent_within_randomized",
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--performance", type=Path, default=DEFAULT_PERFORMANCE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--panel", choices=("pilot", "full"), default="pilot")
    parser.add_argument("--steps", nargs="+", type=int, default=[100])
    parser.add_argument("--times", nargs="+", type=float, default=[1.0])
    parser.add_argument("--samples", type=int)
    parser.add_argument("--seed", type=int, default=20260822)
    parser.add_argument("--schedules", nargs="+", default=list(DEFAULT_SCHEDULES))
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Run the benchmark. Without this flag only the manifest is written.",
    )
    return parser.parse_args()


def selected_panel(panel: str) -> tuple[PanelCase, ...]:
    if panel == "pilot":
        return tuple(case for case in PANEL_CASES if case.panel_tier == "pilot")
    if panel == "full":
        return PANEL_CASES
    raise ValueError(f"Unknown panel: {panel}")


def build_manifest(
    performance: pd.DataFrame,
    panel: str,
    repo_root: Path = REPO_ROOT,
) -> pd.DataFrame:
    specs = selected_panel(panel)
    case_ids = [case.case_id for case in specs]
    if len(case_ids) != len(set(case_ids)):
        raise ValueError("Validation panel contains duplicate case IDs.")
    overlap = DISCOVERY_CASE_IDS.intersection(case_ids)
    if overlap:
        raise ValueError(f"Validation panel overlaps discovery cases: {sorted(overlap)}")

    indexed = performance.set_index("case_id", drop=False)
    missing = sorted(set(case_ids).difference(indexed.index))
    if missing:
        raise ValueError(f"Performance table is missing panel cases: {missing}")

    records: list[dict[str, object]] = []
    for spec in specs:
        source = indexed.loc[spec.case_id]
        if isinstance(source, pd.DataFrame):
            raise ValueError(f"Duplicate performance rows for {spec.case_id}.")
        tensor = Path(str(source["tensor_path"]))
        absolute_tensor = tensor if tensor.is_absolute() else repo_root / tensor
        advantage = float(source["fermionic_advantage_factor"])
        observed_outcome = "favorable" if advantage > 1.0 else "negative_control"
        if observed_outcome != spec.expected_outcome:
            raise ValueError(
                f"Panel label mismatch for {spec.case_id}: expected "
                f"{spec.expected_outcome}, observed {observed_outcome}."
            )
        if not bool(source["valid_comparison"]):
            raise ValueError(f"Panel case is not a valid comparison: {spec.case_id}")
        records.append(
            {
                "case_id": spec.case_id,
                "matched_pair": spec.matched_pair,
                "expected_outcome": spec.expected_outcome,
                "panel_tier": spec.panel_tier,
                "molecule": source["molecule"],
                "basis": source["basis"],
                "active_occupied": int(source["active_occupied"]),
                "active_vacant": int(source["active_vacant"]),
                "n_qubits": int(source["n_qubits"]),
                "fermionic_advantage_factor": advantage,
                "historical_jw_signed_one_minus_overlap": float(
                    source["jw_signed_one_minus_overlap"]
                ),
                "historical_jw_magnitude_one_minus_overlap": float(
                    source["jw_magnitude_one_minus_overlap"]
                ),
                "historical_fermionic_one_minus_overlap": float(
                    source["fermionic_aware_one_minus_overlap"]
                ),
                "historical_fermionic_source_csv": source[
                    "fermionic_aware_source_csv"
                ],
                "tensor_path": str(tensor),
                "tensor_exists": absolute_tensor.exists(),
                "bch_discovery_case": False,
            }
        )
    manifest = pd.DataFrame(records)
    if not manifest["tensor_exists"].all():
        unavailable = manifest.loc[~manifest["tensor_exists"], "case_id"].tolist()
        raise FileNotFoundError(f"Panel tensors are unavailable: {unavailable}")
    pair_sizes = manifest.groupby("matched_pair")["n_qubits"].nunique()
    if not (pair_sizes == 1).all():
        raise ValueError("Every matched pair must use one active-space size.")
    pair_outcomes = manifest.groupby("matched_pair")["expected_outcome"].nunique()
    if not (pair_outcomes == 2).all():
        raise ValueError("Every matched pair must contain both outcome classes.")
    return manifest


def benchmark_command(
    python: Path,
    manifest: pd.DataFrame,
    output: Path,
    steps: Sequence[int],
    evolution_time: float,
    samples: int,
    seed: int,
    schedules: Sequence[str],
) -> list[str]:
    command = [
        str(python),
        "analysis/benchmark_fermionic_structure_ablation.py",
    ]
    for tensor in manifest["tensor_path"]:
        command.extend(("--tensor", str(tensor)))
    command.extend(("--steps", *(str(value) for value in steps)))
    command.extend(("--time", str(evolution_time)))
    command.extend(("--samples", str(samples)))
    command.extend(("--seed", str(seed)))
    command.extend(("--schedules", *schedules))
    command.extend(("--output", str(output)))
    return command


def write_commands(path: Path, commands: Sequence[Sequence[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter=" ", quoting=csv.QUOTE_MINIMAL)
        for command in commands:
            writer.writerow(command)


def validate_runtime(python: Path) -> None:
    if not python.exists():
        raise FileNotFoundError(python)
    check = subprocess.run(
        [str(python), "-c", "import numba, openfermion, scipy"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if check.returncode:
        raise RuntimeError(
            f"{python} cannot import the benchmark dependencies:\n"
            f"{check.stderr.strip()}"
        )


def main() -> None:
    args = parse_arguments()
    samples = args.samples if args.samples is not None else (5 if args.panel == "pilot" else 20)
    if samples < 0:
        raise ValueError("--samples cannot be negative.")
    if any(value <= 0 for value in args.steps):
        raise ValueError("--steps must be positive.")
    if any(value <= 0.0 for value in args.times):
        raise ValueError("--times must be positive.")

    performance_path = (
        args.performance
        if args.performance.is_absolute()
        else REPO_ROOT / args.performance
    )
    performance = pd.read_csv(performance_path)
    manifest = build_manifest(performance, args.panel)
    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.outdir / f"{args.panel}_panel_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    output = args.outdir / f"{args.panel}_ablation_results.csv"
    commands = [
        benchmark_command(
            args.python,
            manifest,
            output,
            args.steps,
            evolution_time,
            samples,
            args.seed,
            args.schedules,
        )
        for evolution_time in args.times
    ]
    command_path = args.outdir / f"{args.panel}_benchmark_commands.txt"
    write_commands(command_path, commands)

    print(f"Panel: {args.panel} ({len(manifest)} cases)")
    print(manifest[["case_id", "matched_pair", "expected_outcome", "fermionic_advantage_factor"]].to_string(index=False))
    print(f"Manifest: {manifest_path}")
    print(f"Commands: {command_path}")
    print(f"Benchmark output: {output}")

    if not args.execute:
        print("Dry run only. Add --execute to launch the benchmark.")
        return

    validate_runtime(args.python)
    environment = os.environ.copy()
    environment.setdefault(
        "MPLCONFIGDIR",
        str(Path(tempfile.gettempdir()) / "qhat-matplotlib-cache"),
    )
    environment.setdefault(
        "NUMBA_CACHE_DIR",
        str(Path(tempfile.gettempdir()) / "qhat-numba-cache"),
    )
    for command in commands:
        print("Executing:", " ".join(command), flush=True)
        subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=environment,
            check=True,
        )


if __name__ == "__main__":
    main()
