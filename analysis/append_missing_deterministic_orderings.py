#!/usr/bin/env python3
"""Run only missing deterministic ordering variants and checkpoint them.

This driver never runs randomized color schedules and never rewrites an
existing plot.  It evaluates only the requested deterministic orderings and
stores them in an add-on CSV.  Existing raw-JW rows may be supplied through
reference CSVs so ratio columns can be populated without evolving raw JW again.

Run from the QHAT repository root on the L-sweep branch.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
from typing import Any, Iterable, Sequence

try:
    from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
    from qhat.analysis.benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        load_interaction_operator,
        parse_case_metadata,
    )
    from qhat.analysis.benchmark_b2_active_spaces_matrix_free import warm_up_numba
    from qhat.analysis.benchmark_comparable_diatomics import (
        parse_bond_lengths,
        select_cases,
    )
except ImportError:
    import benchmark_b2_signed_coefficient_baseline as baseline
    from benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        load_interaction_operator,
        parse_case_metadata,
    )
    from benchmark_b2_active_spaces_matrix_free import warm_up_numba
    from benchmark_comparable_diatomics import parse_bond_lengths, select_cases


DEFAULT_MISSING_ORDERINGS = (
    "jw_magnitude_descending_lexicographic",
    "fermionic_signed_coefficient_lexicographic",
    "fermionic_magnitude_descending_lexicographic",
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate only deterministic ordering rows that are not already "
            "present in an add-on CSV."
        )
    )
    parser.add_argument(
        "--library",
        type=Path,
        required=True,
        help="Root containing the selected *.tensors.npz files.",
    )
    parser.add_argument(
        "--molecules",
        nargs="+",
        required=True,
        help="Exact molecule names, for example B-B Li-Li Be-Be F-F.",
    )
    parser.add_argument(
        "--basis",
        required=True,
        help="Exact basis name, for example sto-6g or hgbs-5.",
    )
    parser.add_argument(
        "--qubits",
        nargs="+",
        type=int,
        required=True,
    )
    parser.add_argument(
        "--bond-length",
        action="append",
        default=[],
        metavar="MOLECULE=VALUE",
    )
    parser.add_argument(
        "--orderings",
        nargs="+",
        default=list(DEFAULT_MISSING_ORDERINGS),
        help="Ordering names to evaluate. Defaults to the three missing rows.",
    )
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument(
        "--time",
        type=float,
        default=1.0,
        dest="evolution_time",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help=(
            "Add-on CSV. Existing rows are preserved; completed ordering rows "
            "are skipped."
        ),
    )
    parser.add_argument(
        "--reference-csv",
        type=Path,
        action="append",
        default=[],
        help=(
            "Existing baseline CSV used only to retrieve raw-JW metrics. "
            "Repeat as needed. These files are never modified."
        ),
    )
    parser.add_argument(
        "--reference-glob",
        default="analysis/*baseline_results.csv",
        help=(
            "Glob for additional read-only baseline CSVs. Set to an empty "
            "string to disable."
        ),
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
        "--dry-run",
        action="store_true",
        help="Print missing rows without running Trotter evolution.",
    )
    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    if not args.library.exists():
        raise FileNotFoundError(args.library)
    if not args.molecules:
        raise ValueError("--molecules cannot be empty.")
    if not args.qubits or any(value <= 0 for value in args.qubits):
        raise ValueError("--qubits must contain positive integers.")
    if args.steps <= 0:
        raise ValueError("--steps must be positive.")
    if args.evolution_time <= 0.0:
        raise ValueError("--time must be positive.")
    if args.tolerance <= 0.0:
        raise ValueError("--tolerance must be positive.")
    if args.parallel_threshold <= 0:
        raise ValueError("--parallel-threshold must be positive.")

    unknown = set(args.orderings).difference(baseline.ORDERING_NAMES)
    if unknown:
        raise ValueError(
            "Unsupported ordering names: " + ", ".join(sorted(unknown))
        )
    if "jw_raw" in args.orderings or "signed_coefficient_lexicographic" in args.orderings:
        print(
            "Warning: a previously completed ordering was explicitly selected; "
            "it will run only when absent from the add-on CSV."
        )


def read_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            return []
        missing_fields = set(baseline.FIELDNAMES).difference(reader.fieldnames)
        if missing_fields:
            raise ValueError(
                f"{path} is missing expected columns: "
                + ", ".join(sorted(missing_fields))
            )
        return list(reader)


def write_rows_atomic(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=baseline.FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def float_matches(value: Any, expected: float, atol: float = 1.0e-14) -> bool:
    try:
        return math.isclose(
            float(value),
            float(expected),
            rel_tol=0.0,
            abs_tol=atol,
        )
    except (TypeError, ValueError):
        return False


def row_matches_settings(row: dict[str, Any], args: argparse.Namespace) -> bool:
    try:
        steps_match = int(row.get("trotter_steps", -1)) == args.steps
    except (TypeError, ValueError):
        return False
    return (
        steps_match
        and float_matches(row.get("evolution_time"), args.evolution_time)
        and float_matches(row.get("coefficient_tolerance"), args.tolerance)
    )


def row_matches_case(
    row: dict[str, Any],
    metadata: Any,
    args: argparse.Namespace,
) -> bool:
    try:
        n_qubits_match = int(row.get("n_qubits", -1)) == metadata.n_qubits
        occupied_match = (
            int(row.get("active_occupied", -1)) == metadata.active_occupied
        )
        vacant_match = (
            int(row.get("active_vacant", -1)) == metadata.active_vacant
        )
    except (TypeError, ValueError):
        return False

    return (
        row.get("status") == "success"
        and str(row.get("molecule", "")) == str(metadata.molecule)
        and str(row.get("basis", "")).lower() == str(metadata.basis).lower()
        and n_qubits_match
        and occupied_match
        and vacant_match
        and float_matches(row.get("bond_length"), metadata.bond_length)
        and row_matches_settings(row, args)
    )


def existing_orderings_for_case(
    rows: Sequence[dict[str, Any]],
    metadata: Any,
    args: argparse.Namespace,
) -> set[str]:
    return {
        str(row.get("ordering", ""))
        for row in rows
        if row_matches_case(row, metadata, args)
    }


def remove_failed_or_duplicate_rows(
    rows: Iterable[dict[str, Any]],
    metadata: Any,
    args: argparse.Namespace,
    orderings: set[str],
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for row in rows:
        same_metadata = (
            str(row.get("molecule", "")) == str(metadata.molecule)
            and str(row.get("basis", "")).lower() == str(metadata.basis).lower()
            and str(row.get("n_qubits", "")) == str(metadata.n_qubits)
            and str(row.get("active_occupied", ""))
            == str(metadata.active_occupied)
            and str(row.get("active_vacant", "")) == str(metadata.active_vacant)
            and float_matches(row.get("bond_length"), metadata.bond_length)
            and row_matches_settings(row, args)
        )
        if same_metadata and str(row.get("ordering", "")) in orderings:
            continue
        result.append(row)
    return result


def optional_float(row: dict[str, Any], name: str) -> float | None:
    try:
        value = float(row.get(name, ""))
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def raw_reference_for_case(
    rows: Sequence[dict[str, Any]],
    metadata: Any,
    args: argparse.Namespace,
) -> dict[str, float]:
    for row in rows:
        if not row_matches_case(row, metadata, args):
            continue
        if row.get("ordering") != "jw_raw":
            continue

        result: dict[str, float] = {}
        for name in (
            "one_minus_overlap",
            "state_infidelity",
            "bch2_hf_state_norm",
        ):
            value = optional_float(row, name)
            if value is not None:
                result[name] = value

        if "one_minus_overlap" not in result:
            overlap = optional_float(row, "state_overlap_abs")
            if overlap is not None:
                result["one_minus_overlap"] = 1.0 - overlap
        return result
    return {}


def discover_reference_paths(args: argparse.Namespace) -> list[Path]:
    paths: list[Path] = []
    seen: set[Path] = set()

    candidates = list(args.reference_csv)
    if args.reference_glob:
        candidates.extend(sorted(Path(".").glob(args.reference_glob)))

    for path in candidates:
        resolved = path.resolve()
        if resolved == args.output.resolve():
            continue
        if path.exists() and path.stat().st_size > 0 and resolved not in seen:
            seen.add(resolved)
            paths.append(path)
    return paths


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    bond_lengths = parse_bond_lengths(args.bond_length)

    tensor_paths = select_cases(
        library=args.library,
        molecules=args.molecules,
        basis_name=args.basis,
        qubits=args.qubits,
        bond_lengths=bond_lengths,
    )

    output_rows = read_rows(args.output)
    reference_paths = discover_reference_paths(args)
    reference_rows: list[dict[str, Any]] = []
    for path in reference_paths:
        reference_rows.extend(read_rows(path))

    print("=" * 100)
    print("Missing-only deterministic ordering run")
    print("=" * 100)
    print(f"Library:              {args.library}")
    print(f"Molecules:            {args.molecules}")
    print(f"Basis:                {args.basis}")
    print(f"Selected cases:       {len(tensor_paths)}")
    print(f"Add-on output:        {args.output}")
    print(f"Read-only references: {len(reference_paths)} CSV files")
    print("Requested orderings:")
    for ordering in args.orderings:
        print(f"  - {ordering}")

    plans: list[tuple[Path, Any, list[str], dict[str, float]]] = []
    all_reference_rows = [*output_rows, *reference_rows]

    for tensor_path in tensor_paths:
        interaction, n_qubits = load_interaction_operator(tensor_path)
        del interaction
        metadata = parse_case_metadata(tensor_path, n_qubits)
        completed = existing_orderings_for_case(output_rows, metadata, args)
        missing = [
            ordering
            for ordering in args.orderings
            if ordering not in completed
        ]
        raw_reference = raw_reference_for_case(
            all_reference_rows,
            metadata,
            args,
        )
        plans.append((tensor_path, metadata, missing, raw_reference))

    for tensor_path, metadata, missing, raw_reference in plans:
        print()
        print(
            f"{metadata.molecule} {metadata.basis} {metadata.n_qubits}q "
            f"({metadata.active_occupied}+{metadata.active_vacant})"
        )
        if missing:
            print("  Missing: " + ", ".join(missing))
        else:
            print("  Complete in add-on CSV; nothing to run.")
        if raw_reference:
            print("  Found existing raw-JW reference metrics.")
        else:
            print("  Raw-JW reference not found; ratio columns will be blank.")

    if args.dry_run:
        return

    warm_up_numba()

    for case_index, (tensor_path, metadata, missing, raw_reference) in enumerate(
        plans,
        start=1,
    ):
        if not missing:
            continue

        print()
        print("-" * 100)
        print(
            f"[{case_index}/{len(plans)}] "
            f"{metadata.molecule} {metadata.basis} {metadata.n_qubits}q"
        )

        missing_set = set(missing)
        output_rows = remove_failed_or_duplicate_rows(
            output_rows,
            metadata,
            args,
            missing_set,
        )

        try:
            new_rows = baseline.benchmark_case(
                tensor_path=tensor_path,
                args=args,
                ordering_names=missing,
                raw_reference=raw_reference,
            )
        except Exception as error:
            print(f"  FAILED: {type(error).__name__}: {error}")
            raise

        output_rows.extend(new_rows)
        write_rows_atomic(args.output, output_rows)
        print(f"  Checkpointed {len(new_rows)} new rows to {args.output}")

    write_rows_atomic(args.output, output_rows)
    print()
    print(f"Done. Existing plots and reference CSVs were not modified.")
    print(f"Add-on CSV: {args.output}")


if __name__ == "__main__":
    main()
