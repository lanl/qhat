#!/usr/bin/env python3
"""Compare signed-ascending and magnitude-descending JW baselines.

The script scans deterministic result CSVs below ``analysis/`` and joins rows
using the exact tensor path and numerical settings.  It writes one auditable
case-level table and one molecule-level summary.  Very small errors can be
classified below a configurable numerical floor instead of being interpreted
as meaningful wins or losses.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


JW_SIGNED = "signed_coefficient_lexicographic"
JW_MAGNITUDE = "jw_magnitude_descending_lexicographic"
FERMIONIC_SIGNED = "fermionic_signed_coefficient_lexicographic"
ORDERINGS = (JW_SIGNED, JW_MAGNITUDE, FERMIONIC_SIGNED)

CASE_FIELDS = [
    "case_id",
    "tensor_path",
    "molecule",
    "bond_length",
    "basis",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "trotter_steps",
    "evolution_time",
    "coefficient_tolerance",
    "jw_signed_one_minus_overlap",
    "jw_magnitude_one_minus_overlap",
    "jw_magnitude_ratio_to_signed",
    "jw_baseline_winner",
    "jw_improvement_factor",
    "below_numerical_floor",
    "fermionic_signed_one_minus_overlap",
    "strongest_jw_one_minus_overlap",
    "strongest_jw_ratio_to_fermionic",
    "fermionic_vs_strongest_jw_winner",
    "jw_signed_source_csv",
    "jw_magnitude_source_csv",
    "fermionic_signed_source_csv",
]

SUMMARY_FIELDS = [
    "molecule",
    "case_pairs_total",
    "case_pairs_above_floor",
    "case_pairs_below_floor",
    "jw_magnitude_wins",
    "jw_signed_wins",
    "jw_ties",
    "jw_magnitude_win_fraction",
    "median_jw_magnitude_ratio_to_signed",
    "geometric_mean_jw_magnitude_ratio_to_signed",
    "fermionic_comparisons_above_floor",
    "fermionic_wins_vs_strongest_jw",
    "strongest_jw_wins_vs_fermionic",
    "fermionic_ties_vs_strongest_jw",
    "fermionic_win_fraction_vs_strongest_jw",
    "median_strongest_jw_ratio_to_fermionic",
    "numerical_error_floor",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("analysis"))
    parser.add_argument(
        "--case-output",
        type=Path,
        default=Path("analysis/jw_magnitude_descending_case_comparison.csv"),
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=Path("analysis/jw_magnitude_descending_molecule_summary.csv"),
    )
    parser.add_argument("--error-floor", type=float, default=1.0e-12)
    return parser.parse_args()


def result_key(row: dict[str, str]) -> tuple[str, int, float, float]:
    return (
        row["tensor_path"],
        int(row["trotter_steps"]),
        float(row["evolution_time"]),
        float(row.get("coefficient_tolerance") or 1.0e-12),
    )


def close_enough(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=1.0e-6, abs_tol=1.0e-14)


def winner(left: float, right: float, left_name: str, right_name: str) -> str:
    if close_enough(left, right):
        return "tie"
    return left_name if left < right else right_name


def read_deterministic_rows(
    input_dir: Path,
    excluded_paths: set[Path],
) -> dict[tuple[str, int, float, float], dict[str, tuple[dict[str, str], Path]]]:
    grouped: dict[
        tuple[str, int, float, float],
        dict[str, tuple[dict[str, str], Path]],
    ] = defaultdict(dict)

    for path in sorted(input_dir.rglob("*.csv")):
        if path.resolve() in excluded_paths:
            continue
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            fields = set(reader.fieldnames or [])
            required = {
                "status",
                "tensor_path",
                "ordering",
                "trotter_steps",
                "evolution_time",
                "one_minus_overlap",
            }
            if not required.issubset(fields):
                continue

            for row in reader:
                ordering = row.get("ordering", "")
                if row.get("status") != "success" or ordering not in ORDERINGS:
                    continue
                try:
                    key = result_key(row)
                    error = float(row["one_minus_overlap"])
                except (KeyError, TypeError, ValueError):
                    continue
                if not math.isfinite(error) or error < 0.0:
                    continue

                existing = grouped[key].get(ordering)
                if existing is not None:
                    existing_error = float(existing[0]["one_minus_overlap"])
                    if not close_enough(error, existing_error):
                        raise ValueError(
                            "Conflicting duplicate deterministic results for "
                            f"{key}, {ordering}: {existing[1]} vs {path}"
                        )
                    continue
                grouped[key][ordering] = (row, path)

    return grouped


def optional_error(
    schedules: dict[str, tuple[dict[str, str], Path]],
    ordering: str,
) -> float | None:
    entry = schedules.get(ordering)
    return float(entry[0]["one_minus_overlap"]) if entry else None


def optional_source(
    schedules: dict[str, tuple[dict[str, str], Path]],
    ordering: str,
) -> str:
    entry = schedules.get(ordering)
    return str(entry[1]) if entry else ""


def build_case_rows(
    grouped: dict[
        tuple[str, int, float, float],
        dict[str, tuple[dict[str, str], Path]],
    ],
    error_floor: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for key, schedules in grouped.items():
        if JW_SIGNED not in schedules or JW_MAGNITUDE not in schedules:
            continue

        signed_row = schedules[JW_SIGNED][0]
        signed = float(signed_row["one_minus_overlap"])
        magnitude = float(schedules[JW_MAGNITUDE][0]["one_minus_overlap"])
        below_floor = min(signed, magnitude) <= error_floor
        ratio = magnitude / signed if signed > 0.0 else math.nan
        jw_winner = (
            "below_numerical_floor"
            if below_floor
            else winner(magnitude, signed, "jw_magnitude", "jw_signed")
        )

        fermionic = optional_error(schedules, FERMIONIC_SIGNED)
        strongest_jw = min(signed, magnitude)
        fermionic_ratio = math.nan
        fermionic_winner = "not_available"
        if fermionic is not None:
            if min(fermionic, strongest_jw) <= error_floor:
                fermionic_winner = "below_numerical_floor"
            elif fermionic > 0.0:
                fermionic_ratio = strongest_jw / fermionic
                fermionic_winner = winner(
                    fermionic,
                    strongest_jw,
                    "fermionic_signed",
                    "strongest_jw",
                )

        rows.append(
            {
                "case_id": signed_row.get("case_id", ""),
                "tensor_path": key[0],
                "molecule": signed_row.get("molecule", ""),
                "bond_length": signed_row.get("bond_length", ""),
                "basis": signed_row.get("basis", ""),
                "active_occupied": signed_row.get("active_occupied", ""),
                "active_vacant": signed_row.get("active_vacant", ""),
                "n_qubits": signed_row.get("n_qubits", ""),
                "trotter_steps": key[1],
                "evolution_time": key[2],
                "coefficient_tolerance": key[3],
                "jw_signed_one_minus_overlap": signed,
                "jw_magnitude_one_minus_overlap": magnitude,
                "jw_magnitude_ratio_to_signed": ratio,
                "jw_baseline_winner": jw_winner,
                "jw_improvement_factor": (
                    max(signed, magnitude) / min(signed, magnitude)
                    if min(signed, magnitude) > 0.0
                    else math.nan
                ),
                "below_numerical_floor": below_floor,
                "fermionic_signed_one_minus_overlap": (
                    fermionic if fermionic is not None else ""
                ),
                "strongest_jw_one_minus_overlap": strongest_jw,
                "strongest_jw_ratio_to_fermionic": fermionic_ratio,
                "fermionic_vs_strongest_jw_winner": fermionic_winner,
                "jw_signed_source_csv": optional_source(schedules, JW_SIGNED),
                "jw_magnitude_source_csv": optional_source(
                    schedules,
                    JW_MAGNITUDE,
                ),
                "fermionic_signed_source_csv": optional_source(
                    schedules,
                    FERMIONIC_SIGNED,
                ),
            }
        )

    rows.sort(
        key=lambda row: (
            str(row["molecule"]),
            str(row["basis"]),
            str(row["bond_length"]),
            int(row["active_occupied"]),
            int(row["active_vacant"]),
            int(row["trotter_steps"]),
            float(row["evolution_time"]),
            str(row["tensor_path"]),
        )
    )
    return rows


def geometric_mean(values: Iterable[float]) -> float:
    values = list(values)
    return math.exp(sum(math.log(value) for value in values) / len(values))


def summarize_group(
    molecule: str,
    rows: list[dict[str, Any]],
    error_floor: float,
) -> dict[str, Any]:
    valid = [row for row in rows if not row["below_numerical_floor"]]
    ratios = [float(row["jw_magnitude_ratio_to_signed"]) for row in valid]
    magnitude_wins = sum(row["jw_baseline_winner"] == "jw_magnitude" for row in valid)
    signed_wins = sum(row["jw_baseline_winner"] == "jw_signed" for row in valid)
    jw_ties = sum(row["jw_baseline_winner"] == "tie" for row in valid)

    fermionic_valid = [
        row
        for row in rows
        if row["fermionic_vs_strongest_jw_winner"]
        in {"fermionic_signed", "strongest_jw", "tie"}
    ]
    fermionic_ratios = [
        float(row["strongest_jw_ratio_to_fermionic"])
        for row in fermionic_valid
    ]
    fermionic_wins = sum(
        row["fermionic_vs_strongest_jw_winner"] == "fermionic_signed"
        for row in fermionic_valid
    )
    strongest_jw_wins = sum(
        row["fermionic_vs_strongest_jw_winner"] == "strongest_jw"
        for row in fermionic_valid
    )
    fermionic_ties = sum(
        row["fermionic_vs_strongest_jw_winner"] == "tie"
        for row in fermionic_valid
    )

    return {
        "molecule": molecule,
        "case_pairs_total": len(rows),
        "case_pairs_above_floor": len(valid),
        "case_pairs_below_floor": len(rows) - len(valid),
        "jw_magnitude_wins": magnitude_wins,
        "jw_signed_wins": signed_wins,
        "jw_ties": jw_ties,
        "jw_magnitude_win_fraction": magnitude_wins / len(valid) if valid else "",
        "median_jw_magnitude_ratio_to_signed": (
            statistics.median(ratios) if ratios else ""
        ),
        "geometric_mean_jw_magnitude_ratio_to_signed": (
            geometric_mean(ratios) if ratios else ""
        ),
        "fermionic_comparisons_above_floor": len(fermionic_valid),
        "fermionic_wins_vs_strongest_jw": fermionic_wins,
        "strongest_jw_wins_vs_fermionic": strongest_jw_wins,
        "fermionic_ties_vs_strongest_jw": fermionic_ties,
        "fermionic_win_fraction_vs_strongest_jw": (
            fermionic_wins / len(fermionic_valid) if fermionic_valid else ""
        ),
        "median_strongest_jw_ratio_to_fermionic": (
            statistics.median(fermionic_ratios) if fermionic_ratios else ""
        ),
        "numerical_error_floor": error_floor,
    }


def build_summary_rows(
    case_rows: list[dict[str, Any]],
    error_floor: float,
) -> list[dict[str, Any]]:
    by_molecule: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in case_rows:
        by_molecule[str(row["molecule"])].append(row)
    summaries = [summarize_group("ALL", case_rows, error_floor)]
    summaries.extend(
        summarize_group(molecule, rows, error_floor)
        for molecule, rows in sorted(by_molecule.items())
    )
    return summaries


def write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_arguments()
    if args.error_floor < 0.0:
        raise ValueError("--error-floor must be non-negative")
    excluded = {args.case_output.resolve(), args.summary_output.resolve()}
    grouped = read_deterministic_rows(args.input_dir, excluded)
    case_rows = build_case_rows(grouped, args.error_floor)
    summary_rows = build_summary_rows(case_rows, args.error_floor)
    write_csv(args.case_output, CASE_FIELDS, case_rows)
    write_csv(args.summary_output, SUMMARY_FIELDS, summary_rows)

    overall = summary_rows[0]
    print(f"Matched signed/JW-magnitude cases: {len(case_rows)}")
    print(
        "Above numerical floor: "
        f"{overall['case_pairs_above_floor']} "
        f"(floor={args.error_floor:.1e})"
    )
    print(
        "JW magnitude wins / signed wins / ties: "
        f"{overall['jw_magnitude_wins']} / "
        f"{overall['jw_signed_wins']} / {overall['jw_ties']}"
    )
    print(
        "Fermionic signed wins vs strongest JW: "
        f"{overall['fermionic_wins_vs_strongest_jw']} / "
        f"{overall['fermionic_comparisons_above_floor']}"
    )
    print(f"Case output:    {args.case_output}")
    print(f"Summary output: {args.summary_output}")


if __name__ == "__main__":
    main()
