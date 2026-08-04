#!/usr/bin/env python3
"""Generate one Trotter-error figure for every L-sweep Hamiltonian case.

This is a batch driver for ``analysis/visualize_L_sweep_trotter.py``.
It loads the benchmark CSV files once, discovers all unique ``case_id``
values, and calls ``plot_case_grid`` for every case.

Example
-------
python analysis/visualize_all_L_sweep_cases.py \
    --csv \
        analysis/l_sweep_trotter_state_t1.csv \
        analysis/l_sweep_trotter_state_t2.csv \
        analysis/l_sweep_trotter_state_t5.csv \
        analysis/l_sweep_trotter_state_dense_t1.csv \
        analysis/l_sweep_trotter_state_dense_t2.csv \
        analysis/l_sweep_trotter_state_dense_t5.csv \
    --metric state_infidelity \
    --expected-cases 153 \
    --output-dir analysis/l_sweep_all_case_figures
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from visualize_L_sweep_trotter import (
    default_floors,
    load_results,
    plot_case_grid,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create one visualization for every L-sweep benchmark case."
    )
    parser.add_argument(
        "--csv",
        type=Path,
        nargs="+",
        required=True,
        help="Benchmark CSV files to combine before plotting.",
    )
    parser.add_argument(
        "--metric",
        choices=("state_infidelity", "state_vector_2norm_error"),
        default="state_infidelity",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis/l_sweep_all_case_figures"),
    )
    parser.add_argument(
        "--expected-cases",
        type=int,
        default=153,
        help="Expected number of valid L-sweep cases.",
    )
    parser.add_argument(
        "--plot-floor",
        type=float,
        default=None,
        help=(
            "Positive plotting floor for logarithmic y-axes. "
            "Defaults to the value used by visualize_L_sweep_trotter.py."
        ),
    )
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Do not regenerate a figure when its output file already exists.",
    )
    return parser.parse_args()


def expected_output_path(
    output_root: Path,
    molecule: str,
    case_id: str,
    metric: str,
) -> Path:
    safe_case = "".join(
        character if character.isalnum() or character in "_.-" else "_"
        for character in case_id
    )
    return output_root / "cases" / molecule / f"case_{safe_case}_{metric}.png"


def main() -> None:
    args = parse_args()

    _, default_plot_floor = default_floors(args.metric)
    plot_floor = (
        args.plot_floor if args.plot_floor is not None else default_plot_floor
    )

    data = load_results(args.csv)
    case_metadata = (
        data[
            [
                "case_id",
                "molecule",
                "basis",
                "bond_length",
                "active_occupied",
                "active_vacant",
                "n_qubits",
            ]
        ]
        .drop_duplicates("case_id")
        .sort_values(
            [
                "molecule",
                "n_qubits",
                "basis",
                "bond_length",
                "active_occupied",
                "active_vacant",
                "case_id",
            ]
        )
        .reset_index(drop=True)
    )

    number_of_cases = len(case_metadata)
    print(f"Loaded {len(data):,} successful benchmark rows.")
    print(f"Discovered {number_of_cases} unique Hamiltonian cases.")

    if number_of_cases != args.expected_cases:
        print(
            "WARNING: expected "
            f"{args.expected_cases} cases but discovered {number_of_cases}."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    global_times = sorted(
        float(value) for value in data["evolution_time"].dropna().unique()
    )
    global_formulas = sorted(
        int(value) for value in data["formula_order"].dropna().unique()
    )
    global_steps = sorted(
        int(value) for value in data["trotter_steps"].dropna().unique()
    )
    global_orderings = sorted(data["ordering"].dropna().unique())
    global_factorization_levels = sorted(
        data["factorization_level"].dropna().unique()
    )

    expected_rows_per_case = (
        len(global_times)
        * len(global_formulas)
        * len(global_steps)
        * len(global_orderings)
    )

    manifest_rows: list[dict[str, object]] = []
    generated = 0
    skipped = 0
    failed = 0

    for index, metadata in case_metadata.iterrows():
        case_id = str(metadata["case_id"])
        molecule = str(metadata["molecule"])
        case_output_dir = args.output_dir / "cases" / molecule
        case_output_dir.mkdir(parents=True, exist_ok=True)

        output_path = expected_output_path(
            args.output_dir,
            molecule,
            case_id,
            args.metric,
        )

        subset = data[data["case_id"].eq(case_id)].copy()
        observed_rows = len(subset)
        duplicate_key_columns = [
            "case_id",
            "evolution_time",
            "formula_order",
            "trotter_steps",
            "ordering",
        ]
        unique_rows = len(subset.drop_duplicates(duplicate_key_columns))
        missing_rows = max(expected_rows_per_case - unique_rows, 0)

        row: dict[str, object] = {
            "case_id": case_id,
            "molecule": molecule,
            "basis": metadata["basis"],
            "bond_length": metadata["bond_length"],
            "active_occupied": int(metadata["active_occupied"]),
            "active_vacant": int(metadata["active_vacant"]),
            "n_qubits": int(metadata["n_qubits"]),
            "evolution_times": json.dumps(
                sorted(
                    float(value)
                    for value in subset["evolution_time"].dropna().unique()
                )
            ),
            "formula_orders": json.dumps(
                sorted(
                    int(value)
                    for value in subset["formula_order"].dropna().unique()
                )
            ),
            "trotter_steps": json.dumps(
                sorted(
                    int(value)
                    for value in subset["trotter_steps"].dropna().unique()
                )
            ),
            "orderings": json.dumps(
                sorted(str(value) for value in subset["ordering"].dropna().unique())
            ),
            "factorization_levels": json.dumps(
                sorted(
                    str(value)
                    for value in subset["factorization_level"].dropna().unique()
                )
            ),
            "observed_rows": observed_rows,
            "unique_grid_rows": unique_rows,
            "expected_rows_from_global_grid": expected_rows_per_case,
            "missing_rows_from_global_grid": missing_rows,
            "grid_complete": missing_rows == 0,
            "status": "",
            "error_message": "",
            "output_path": str(output_path),
        }

        print(f"[{index + 1}/{number_of_cases}] {case_id}")

        if args.skip_existing and output_path.exists():
            row["status"] = "skipped_existing"
            skipped += 1
            manifest_rows.append(row)
            continue

        try:
            created_path = plot_case_grid(
                data=data,
                case_id=case_id,
                metric=args.metric,
                plot_floor=plot_floor,
                output_dir=case_output_dir,
                dpi=args.dpi,
            )
            row["status"] = "generated"
            row["output_path"] = str(created_path)
            generated += 1
        except Exception as error:
            row["status"] = "failed"
            row["error_message"] = f"{type(error).__name__}: {error}"
            failed += 1
            print(f"  FAILED: {row['error_message']}")

        manifest_rows.append(row)

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = (
        args.output_dir / f"all_case_plot_manifest_{args.metric}.csv"
    )
    manifest.to_csv(manifest_path, index=False)

    print()
    print("=" * 72)
    print(f"Cases discovered:       {number_of_cases}")
    print(
        "Factorization levels:   "
        + ", ".join(str(value) for value in global_factorization_levels)
    )
    print(f"Figures generated:      {generated}")
    print(f"Existing figures kept:  {skipped}")
    print(f"Failures:               {failed}")
    print(f"Manifest:               {manifest_path}")
    print(f"Output root:            {args.output_dir}")
    print("=" * 72)


if __name__ == "__main__":
    main()
