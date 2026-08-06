#!/usr/bin/env python3
"""Append-only benchmark for magnitude ordering inside fixed color groups.

This script is intentionally separate from the existing QHAT L-sweep benchmark
files. It does not modify or overwrite the current coloring, random-schedule,
baseline, or plotting outputs.

For each selected tensor case, it constructs the same deterministic greedy-color
partitions used by ``benchmark_b2_coloring_robustness.py`` and evaluates only two
new schedules:

1. ``jw:magnitude_within_groups``
   Keep the JW color-group order fixed and sort Pauli terms inside each group by
   decreasing absolute final JW coefficient. Exact magnitude ties are broken by
   ascending dense Pauli-string lexicographic order.

2. ``fermionic:magnitude_within_groups``
   Keep the fermionic color-group order fixed and sort complete Hermitian
   fermionic terms inside each group by decreasing coefficient magnitude. Exact
   magnitude ties are broken by the term's component-key lexicographic order.
   The sorted fermionic order then induces a permutation of the final combined
   JW Pauli strings using the repository's existing first-occurrence rule.

The output CSV is append-only. Existing successful rows with the same case,
graph level, schedule, Trotter steps, and evolution time are skipped.

Run this file from the QHAT repository root after placing it in ``analysis/``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from openfermion import get_fermion_operator, jordan_wigner

try:
    from qhat.analysis.benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        build_fermionic_noncommutation_graph,
        build_hermitian_fermion_terms,
        build_pauli_noncommutation_graph,
        clean_fermion_operator,
        discover_tensor_paths,
        load_interaction_operator,
        parse_case_metadata,
        real_coefficient,
        validate_pauli_order,
    )
except ImportError:
    from benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        build_fermionic_noncommutation_graph,
        build_hermitian_fermion_terms,
        build_pauli_noncommutation_graph,
        clean_fermion_operator,
        discover_tensor_paths,
        load_interaction_operator,
        parse_case_metadata,
        real_coefficient,
        validate_pauli_order,
    )

try:
    from qhat.analysis.benchmark_b2_active_spaces_matrix_free import (
        build_hartree_fock_state,
        compare_states,
        compile_ordered_terms,
        exact_reference_state,
        evolve_trotter_state,
        number_sector_basis_indices,
        warm_up_numba,
    )
except ImportError:
    from benchmark_b2_active_spaces_matrix_free import (
        build_hartree_fock_state,
        compare_states,
        compile_ordered_terms,
        exact_reference_state,
        evolve_trotter_state,
        number_sector_basis_indices,
        warm_up_numba,
    )

try:
    from qhat.analysis.benchmark_b2_coloring_robustness import (
        HFCommutatorEvaluator,
        deterministic_color_groups,
        graph_density,
        induced_pauli_order_indices,
        precompute_fermion_to_pauli_indices,
    )
except ImportError:
    from benchmark_b2_coloring_robustness import (
        HFCommutatorEvaluator,
        deterministic_color_groups,
        graph_density,
        induced_pauli_order_indices,
        precompute_fermion_to_pauli_indices,
    )


PauliKey = tuple[tuple[int, str], ...]
SCHEDULE_NAME = "magnitude_within_groups"
GRAPH_LEVELS = ("jw", "fermionic")

FIELDNAMES = [
    "status",
    "error_message",
    "case_id",
    "tensor_path",
    "molecule",
    "bond_length",
    "basis",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "number_of_fermionic_terms",
    "number_of_pauli_terms",
    "graph_level",
    "schedule",
    "ordering_definition",
    "graph_vertices",
    "graph_edges",
    "graph_density",
    "number_of_colors",
    "group_sizes",
    "group_order",
    "original_within_group_order_hash",
    "sorted_within_group_order_hash",
    "pauli_order_hash",
    "trotter_steps",
    "trotter_dt",
    "evolution_time",
    "nominal_exponential_count",
    "bch2_hf_state_norm",
    "leading_state_error_norm_estimate",
    "exact_sector_dimension",
    "exact_build_time_seconds",
    "exact_evolution_time_seconds",
    "schedule_build_time_seconds",
    "bch_evaluation_time_seconds",
    "trotter_runtime_seconds",
    "state_overlap_abs",
    "one_minus_overlap",
    "state_infidelity",
    "state_vector_2norm_error",
    "phase_aligned_state_2norm_error",
    "particle_number_leakage",
    "spin_sector_leakage",
    "coefficient_tolerance",
]

ResumeKey = tuple[str, str, str, int, float]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate only magnitude-within-fixed-color-group schedules for "
            "JW and fermionic coloring, using an append-only output CSV."
        )
    )
    parser.add_argument(
        "--library",
        type=Path,
        required=True,
        help="Root containing the selected *.tensors.npz files.",
    )
    parser.add_argument(
        "--molecule",
        required=True,
        help="Exact molecule metadata label, for example Li-Li or F-F.",
    )
    parser.add_argument(
        "--basis",
        required=True,
        help="Basis metadata label, for example sto-6g or hgbs-5.",
    )
    parser.add_argument(
        "--bond-length",
        required=True,
        help="Bond length to select, for example 2.66.",
    )
    parser.add_argument(
        "--qubits",
        type=int,
        nargs="+",
        required=True,
        help="Requested active-space qubit counts.",
    )
    parser.add_argument(
        "--graph-level",
        choices=("both", "jw", "fermionic"),
        default="both",
        help="Evaluate both graph levels or only one. Default: both.",
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
        "--output",
        type=Path,
        default=Path("analysis/magnitude_within_color_groups_results.csv"),
        help=(
            "Separate append-only result CSV. Existing successful rows are "
            "skipped."
        ),
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help="Coefficient and commutator tolerance.",
    )
    parser.add_argument(
        "--parallel-threshold",
        type=int,
        default=2**16,
        help="Use the parallel Pauli kernel at or above this vector size.",
    )
    parser.add_argument(
        "--no-spin-sector",
        action="store_true",
        help="Use only particle-number restriction for the exact reference.",
    )
    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    if not args.library.exists():
        raise FileNotFoundError(args.library)
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


def bond_lengths_match(left: Any, right: str) -> bool:
    try:
        return math.isclose(
            float(left),
            float(right),
            rel_tol=0.0,
            abs_tol=1.0e-10,
        )
    except (TypeError, ValueError):
        return str(left) == str(right)


def select_cases(
    library: Path,
    molecule: str,
    basis: str,
    bond_length: str,
    requested_qubits: Sequence[int],
) -> list[Path]:
    requested = set(requested_qubits)
    selected: dict[int, list[Path]] = defaultdict(list)

    tensor_paths = discover_tensor_paths(library, case_pattern=None, limit=None)
    if not tensor_paths:
        raise FileNotFoundError(
            f"No *.tensors.npz files were found under {library}."
        )

    for tensor_path in tensor_paths:
        try:
            interaction, n_qubits = load_interaction_operator(tensor_path)
            del interaction
            metadata = parse_case_metadata(tensor_path, n_qubits)
        except Exception:
            continue

        if metadata.molecule != molecule:
            continue
        if str(metadata.basis).lower() != basis.lower():
            continue
        if not bond_lengths_match(metadata.bond_length, bond_length):
            continue
        if n_qubits in requested:
            selected[n_qubits].append(tensor_path)

    missing = sorted(requested.difference(selected))
    if missing:
        raise FileNotFoundError(
            "No matching tensor file for "
            f"{molecule}/{basis} at bond length {bond_length} and qubits "
            f"{missing}."
        )

    ambiguous = {
        n_qubits: paths
        for n_qubits, paths in selected.items()
        if len(paths) != 1
    }
    if ambiguous:
        details = "\n".join(
            f"  {n_qubits}: {[str(path) for path in paths]}"
            for n_qubits, paths in sorted(ambiguous.items())
        )
        raise ValueError(
            "Expected exactly one tensor file for each requested qubit count. "
            "Separate duplicate active-space families into different runs:\n"
            f"{details}"
        )

    return [selected[n_qubits][0] for n_qubits in sorted(requested)]


def dense_pauli_string(pauli_key: PauliKey, n_qubits: int) -> str:
    operators = dict(pauli_key)
    return "".join(
        operators.get(qubit, "I")
        for qubit in range(n_qubits)
    )


def flatten(groups: Sequence[Sequence[int]]) -> list[int]:
    return [int(node) for group in groups for node in group]


def hash_nested_integers(values: Sequence[Sequence[int]]) -> str:
    payload = json.dumps(
        [list(map(int, group)) for group in values],
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def hash_integer_order(values: Sequence[int]) -> str:
    payload = ",".join(str(int(value)) for value in values)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def fermionic_lexicographic_key(term: Any) -> tuple[Any, ...]:
    return tuple(sorted(term.component_keys))


def fermionic_coefficient_magnitude(term: Any) -> float:
    coefficients = [
        complex(term.operator.terms[key])
        for key in term.component_keys
        if key in term.operator.terms
    ]
    if not coefficients:
        coefficients = [complex(value) for value in term.operator.terms.values()]
    if not coefficients:
        raise ValueError("A Hermitian fermionic term has no coefficients.")
    return max(abs(value) for value in coefficients)


def sorted_jw_groups_by_magnitude(
    groups: Sequence[Sequence[int]],
    raw_pauli_keys: Sequence[PauliKey],
    final_coefficients: dict[PauliKey, complex],
    n_qubits: int,
    tolerance: float,
) -> list[list[int]]:
    def sort_key(node: int) -> tuple[float, str]:
        pauli_key = raw_pauli_keys[node]
        coefficient = real_coefficient(
            final_coefficients[pauli_key],
            tolerance,
        )
        return (
            -abs(coefficient),
            dense_pauli_string(pauli_key, n_qubits),
        )

    return [sorted(map(int, group), key=sort_key) for group in groups]


def sorted_fermionic_groups_by_magnitude(
    groups: Sequence[Sequence[int]],
    hermitian_terms: Sequence[Any],
) -> list[list[int]]:
    def sort_key(node: int) -> tuple[float, tuple[Any, ...]]:
        term = hermitian_terms[node]
        return (
            -fermionic_coefficient_magnitude(term),
            fermionic_lexicographic_key(term),
        )

    return [sorted(map(int, group), key=sort_key) for group in groups]


def resume_key(
    case_id: str,
    graph_level: str,
    schedule: str,
    steps: int,
    evolution_time: float,
) -> ResumeKey:
    return (
        case_id,
        graph_level,
        schedule,
        steps,
        evolution_time,
    )


def load_completed(output_path: Path) -> set[ResumeKey]:
    completed: set[ResumeKey] = set()
    if not output_path.exists() or output_path.stat().st_size == 0:
        return completed

    with output_path.open("r", newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != FIELDNAMES:
            raise ValueError(
                f"Existing output schema does not match this script: {output_path}"
            )
        for row in reader:
            if row.get("status") != "success":
                continue
            try:
                completed.add(
                    resume_key(
                        case_id=row["case_id"],
                        graph_level=row["graph_level"],
                        schedule=row["schedule"],
                        steps=int(row["trotter_steps"]),
                        evolution_time=float(row["evolution_time"]),
                    )
                )
            except (KeyError, TypeError, ValueError):
                continue
    return completed


def blank_row() -> dict[str, Any]:
    return {field: "" for field in FIELDNAMES}


def add_metadata(row: dict[str, Any], metadata: Any, tensor_path: Path) -> None:
    row.update(
        {
            "case_id": metadata.case_id,
            "tensor_path": str(tensor_path),
            "molecule": metadata.molecule,
            "bond_length": metadata.bond_length,
            "basis": metadata.basis,
            "active_occupied": metadata.active_occupied,
            "active_vacant": metadata.active_vacant,
            "n_qubits": metadata.n_qubits,
        }
    )


def ordering_definition(graph_level: str) -> str:
    if graph_level == "jw":
        return (
            "fixed deterministic JW color groups and group order; final JW "
            "Pauli terms sorted inside each group by decreasing absolute "
            "coefficient magnitude, with dense Pauli lexicographic tie-breaking"
        )
    if graph_level == "fermionic":
        return (
            "fixed deterministic fermionic color groups and group order; "
            "complete Hermitian fermionic terms sorted inside each group by "
            "decreasing coefficient magnitude, with component-key "
            "lexicographic tie-breaking; final JW Pauli order induced by first "
            "occurrence"
        )
    raise ValueError(f"Unsupported graph level: {graph_level!r}")


def selected_graph_levels(value: str) -> tuple[str, ...]:
    if value == "both":
        return GRAPH_LEVELS
    return (value,)


def benchmark_case(
    tensor_path: Path,
    args: argparse.Namespace,
    completed: set[ResumeKey],
    writer: csv.DictWriter,
    output_stream: Any,
) -> None:
    interaction, n_qubits = load_interaction_operator(tensor_path)
    metadata = parse_case_metadata(tensor_path, n_qubits)

    requested_levels = selected_graph_levels(args.graph_level)
    pending_levels = [
        level
        for level in requested_levels
        if resume_key(
            metadata.case_id,
            level,
            SCHEDULE_NAME,
            args.steps,
            args.evolution_time,
        )
        not in completed
    ]
    if not pending_levels:
        print("    All requested magnitude-within-group rows already exist; skipping.")
        return

    fermion_hamiltonian = clean_fermion_operator(
        get_fermion_operator(interaction),
        args.tolerance,
    )
    full_jw_hamiltonian = jordan_wigner(fermion_hamiltonian)
    full_jw_hamiltonian.compress(abs_tol=args.tolerance)
    final_coefficients = {
        key: coefficient
        for key, coefficient in full_jw_hamiltonian.terms.items()
        if key != () and abs(coefficient) > args.tolerance
    }
    raw_pauli_keys = list(final_coefficients)
    raw_index_by_key = {
        key: index for index, key in enumerate(raw_pauli_keys)
    }

    print("    Building complete Hermitian fermionic terms...")
    hermitian_terms = build_hermitian_fermion_terms(
        fermion_hamiltonian,
        args.tolerance,
    )

    print("    Building fixed noncommutation graphs and color groups...")
    jw_graph, _ = build_pauli_noncommutation_graph(raw_pauli_keys)
    fermionic_graph, _ = build_fermionic_noncommutation_graph(
        hermitian_terms,
        args.tolerance,
    )
    original_groups = {
        "jw": deterministic_color_groups(jw_graph),
        "fermionic": deterministic_color_groups(fermionic_graph),
    }
    sorted_groups = {
        "jw": sorted_jw_groups_by_magnitude(
            original_groups["jw"],
            raw_pauli_keys,
            final_coefficients,
            n_qubits,
            args.tolerance,
        ),
        "fermionic": sorted_fermionic_groups_by_magnitude(
            original_groups["fermionic"],
            hermitian_terms,
        ),
    }

    fermion_to_pauli_indices = precompute_fermion_to_pauli_indices(
        hermitian_terms=hermitian_terms,
        final_coefficients=final_coefficients,
        raw_index_by_key=raw_index_by_key,
        tolerance=args.tolerance,
    )

    pauli_orders: dict[str, list[int]] = {
        "jw": flatten(sorted_groups["jw"]),
        "fermionic": induced_pauli_order_indices(
            fermionic_node_order=flatten(sorted_groups["fermionic"]),
            fermion_to_pauli_indices=fermion_to_pauli_indices,
            number_of_pauli_terms=len(raw_pauli_keys),
        ),
    }

    for level, indices in pauli_orders.items():
        validate_pauli_order(
            f"{level}_{SCHEDULE_NAME}",
            [raw_pauli_keys[index] for index in indices],
            raw_pauli_keys,
        )

    print("    Precomputing the state-dependent BCH diagnostic...")
    bch_evaluator = HFCommutatorEvaluator.build(
        pauli_keys=raw_pauli_keys,
        coefficients=final_coefficients,
        pauli_graph=jw_graph,
        n_qubits=n_qubits,
        n_electrons=metadata.active_occupied,
        tolerance=args.tolerance,
    )

    print("    Building the exact symmetry-sector reference state...")
    (
        exact_sector_state,
        exact_basis_indices,
        exact_build_time,
        exact_evolution_time,
    ) = exact_reference_state(
        fermion_hamiltonian=fermion_hamiltonian,
        n_qubits=n_qubits,
        n_electrons=metadata.active_occupied,
        evolution_time=args.evolution_time,
        tolerance=args.tolerance,
        spin_preserving=not args.no_spin_sector,
    )
    number_basis_indices = number_sector_basis_indices(
        n_qubits,
        metadata.active_occupied,
    )
    initial_state = build_hartree_fock_state(
        n_qubits,
        metadata.active_occupied,
    )
    compiled_raw_terms = compile_ordered_terms(
        raw_pauli_keys,
        final_coefficients,
        n_qubits,
        args.tolerance,
    )

    graph_by_level = {
        "jw": jw_graph,
        "fermionic": fermionic_graph,
    }

    for level in requested_levels:
        key = resume_key(
            metadata.case_id,
            level,
            SCHEDULE_NAME,
            args.steps,
            args.evolution_time,
        )
        if key in completed:
            print(f"    SKIP {level}:{SCHEDULE_NAME}")
            continue

        print(f"    RUN  {level}:{SCHEDULE_NAME}")
        row = blank_row()
        add_metadata(row, metadata, tensor_path)
        graph = graph_by_level[level]
        groups = original_groups[level]
        magnitude_groups = sorted_groups[level]
        pauli_order_indices = pauli_orders[level]

        row.update(
            {
                "status": "success",
                "number_of_fermionic_terms": len(hermitian_terms),
                "number_of_pauli_terms": len(raw_pauli_keys),
                "graph_level": level,
                "schedule": SCHEDULE_NAME,
                "ordering_definition": ordering_definition(level),
                "graph_vertices": graph.number_of_nodes(),
                "graph_edges": graph.number_of_edges(),
                "graph_density": graph_density(graph),
                "number_of_colors": len(groups),
                "group_sizes": json.dumps(
                    [len(group) for group in groups],
                    separators=(",", ":"),
                ),
                "group_order": json.dumps(
                    list(range(len(groups))),
                    separators=(",", ":"),
                ),
                "original_within_group_order_hash": hash_nested_integers(groups),
                "sorted_within_group_order_hash": hash_nested_integers(
                    magnitude_groups
                ),
                "pauli_order_hash": hash_integer_order(pauli_order_indices),
                "trotter_steps": args.steps,
                "trotter_dt": args.evolution_time / args.steps,
                "evolution_time": args.evolution_time,
                "exact_sector_dimension": exact_sector_state.size,
                "exact_build_time_seconds": exact_build_time,
                "exact_evolution_time_seconds": exact_evolution_time,
                "coefficient_tolerance": args.tolerance,
            }
        )

        try:
            schedule_start = time.perf_counter()
            ordered_terms = [
                compiled_raw_terms[index]
                for index in pauli_order_indices
            ]
            row["schedule_build_time_seconds"] = (
                time.perf_counter() - schedule_start
            )

            bch_start = time.perf_counter()
            bch_norm = bch_evaluator.evaluate(pauli_order_indices)
            row["bch_evaluation_time_seconds"] = (
                time.perf_counter() - bch_start
            )
            row["bch2_hf_state_norm"] = bch_norm
            row["leading_state_error_norm_estimate"] = (
                0.5
                * args.evolution_time**2
                * bch_norm
                / args.steps
            )

            trotter_start = time.perf_counter()
            approximate_state, nominal_exponential_count = evolve_trotter_state(
                initial_state=initial_state,
                terms=ordered_terms,
                formula_order=1,
                trotter_steps=args.steps,
                evolution_time=args.evolution_time,
                parallel_threshold=args.parallel_threshold,
            )
            row["trotter_runtime_seconds"] = (
                time.perf_counter() - trotter_start
            )
            row["nominal_exponential_count"] = nominal_exponential_count

            metrics = compare_states(
                exact_sector_state=exact_sector_state,
                exact_basis_indices=exact_basis_indices,
                approximate_full_state=approximate_state,
                number_basis_indices=number_basis_indices,
            )
            metric_fields = (
                "state_overlap_abs",
                "state_infidelity",
                "state_vector_2norm_error",
                "phase_aligned_state_2norm_error",
                "particle_number_leakage",
                "spin_sector_leakage",
            )
            for field in metric_fields:
                row[field] = metrics.get(field, "")

            overlap = float(metrics["state_overlap_abs"])
            row["one_minus_overlap"] = max(0.0, 1.0 - overlap)

            writer.writerow(row)
            output_stream.flush()
            completed.add(key)

            print(
                "      "
                f"1-overlap={row['one_minus_overlap']:.6e}  "
                f"infidelity={float(row['state_infidelity']):.6e}  "
                f"BCH_HF={bch_norm:.6e}"
            )
        except Exception as error:
            row["status"] = "failed"
            row["error_message"] = f"{type(error).__name__}: {error}"
            writer.writerow(row)
            output_stream.flush()
            print(f"      FAILED: {type(error).__name__}: {error}")


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    warm_up_numba()

    tensor_paths = select_cases(
        library=args.library,
        molecule=args.molecule,
        basis=args.basis,
        bond_length=args.bond_length,
        requested_qubits=args.qubits,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    completed = load_completed(args.output)
    output_exists = args.output.exists() and args.output.stat().st_size > 0

    print("=" * 96)
    print("Magnitude ordering inside fixed color groups")
    print("=" * 96)
    print(f"Molecule:          {args.molecule}")
    print(f"Basis:             {args.basis}")
    print(f"Bond length:       {args.bond_length}")
    print(f"Qubits:            {sorted(args.qubits)}")
    print(f"Graph level:       {args.graph_level}")
    print(f"Trotter steps:     {args.steps}")
    print(f"Evolution time:    {args.evolution_time}")
    print(f"Append-only CSV:   {args.output}")
    print(f"Completed rows:    {len(completed)}")

    with args.output.open("a", newline="", encoding="utf-8") as output_stream:
        writer = csv.DictWriter(output_stream, fieldnames=FIELDNAMES)
        if not output_exists:
            writer.writeheader()
            output_stream.flush()

        for index, tensor_path in enumerate(tensor_paths, start=1):
            print()
            print("-" * 96)
            print(f"[{index}/{len(tensor_paths)}] {tensor_path}")
            try:
                benchmark_case(
                    tensor_path=tensor_path,
                    args=args,
                    completed=completed,
                    writer=writer,
                    output_stream=output_stream,
                )
            except Exception as error:
                print(f"    FAILED CASE: {type(error).__name__}: {error}")
                failure = blank_row()
                failure.update(
                    {
                        "status": "failed",
                        "error_message": f"{type(error).__name__}: {error}",
                        "case_id": tensor_path.name.removesuffix(".tensors.npz"),
                        "tensor_path": str(tensor_path),
                        "molecule": args.molecule,
                        "bond_length": args.bond_length,
                        "basis": args.basis,
                        "trotter_steps": args.steps,
                        "trotter_dt": args.evolution_time / args.steps,
                        "evolution_time": args.evolution_time,
                        "coefficient_tolerance": args.tolerance,
                    }
                )
                writer.writerow(failure)
                output_stream.flush()

    print()
    print(f"Finished. Results were appended to: {args.output}")


if __name__ == "__main__":
    main()
