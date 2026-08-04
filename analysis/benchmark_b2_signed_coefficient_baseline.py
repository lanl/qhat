#!/usr/bin/env python3
"""Evaluate Reuben's clarified fixed baseline for the B2 robustness study.

Baseline definition
-------------------
Order the final combined Jordan-Wigner Pauli terms by:

1. ascending signed coefficient;
2. ascending dense Pauli-string lexicographic order for exact coefficient ties.

The raw Jordan-Wigner insertion order is also evaluated as a consistency check.
This script does not rebuild or rerun the randomized coloring schedules. It writes
one small CSV containing the two deterministic orderings for each selected B2
active space.

Place this file in ``analysis/`` on the QHAT ``L-sweep`` branch and run it from
the repository root.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
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
        HermitianFermionTerm,
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
        HermitianFermionTerm,
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
    )
except ImportError:
    from benchmark_b2_coloring_robustness import HFCommutatorEvaluator


PauliKey = tuple[tuple[int, str], ...]

ORDERING_NAMES = (
    "jw_raw",
    "signed_coefficient_lexicographic",
)

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
    "number_of_pauli_terms",
    "ordering",
    "ordering_definition",
    "pauli_order_hash",
    "first_pauli_string",
    "first_coefficient",
    "last_pauli_string",
    "last_coefficient",
    "trotter_steps",
    "trotter_dt",
    "evolution_time",
    "nominal_exponential_count",
    "bch2_hf_state_norm",
    "bch_squared_ratio_to_raw_jw",
    "exact_sector_dimension",
    "exact_build_time_seconds",
    "exact_evolution_time_seconds",
    "trotter_runtime_seconds",
    "state_overlap_abs",
    "one_minus_overlap",
    "state_infidelity",
    "state_vector_2norm_error",
    "phase_aligned_state_2norm_error",
    "particle_number_leakage",
    "spin_sector_leakage",
    "one_minus_overlap_ratio_to_raw_jw",
    "state_infidelity_ratio_to_raw_jw",
    "coefficient_tolerance",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate raw JW and ascending signed-coefficient/lexicographic "
            "baseline orderings for selected B2/STO-6G active spaces."
        )
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=Path("hamiltonian_generator/b2_active_space_library"),
        help="Root containing B2 *.tensors.npz files.",
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
        default="1.70",
        help="Select one B2 bond-length directory. Default: 1.70.",
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
        default=Path("analysis/b2_signed_coefficient_baseline_results.csv"),
        help="Output CSV for the deterministic baseline calculations.",
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


def select_cases(
    library: Path,
    requested_qubits: Sequence[int],
    bond_length: str,
) -> list[Path]:
    """Select exactly one B2/STO-6G tensor file for each requested size."""
    requested = set(requested_qubits)
    selected: dict[int, list[Path]] = defaultdict(list)

    for path in discover_tensor_paths(library, case_pattern=None, limit=None):
        try:
            interaction, n_qubits = load_interaction_operator(path)
            del interaction
            metadata = parse_case_metadata(path, n_qubits)
        except Exception:
            continue

        if metadata.molecule != "B-B":
            continue
        if metadata.basis.lower() != "sto-6g":
            continue
        if str(metadata.bond_length) != str(bond_length):
            continue
        if n_qubits in requested:
            selected[n_qubits].append(path)

    missing = sorted(requested.difference(selected))
    if missing:
        raise FileNotFoundError(
            "No matching B2/STO-6G tensor files for active-space qubits "
            f"{missing} at bond length {bond_length}."
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
            "Expected exactly one tensor file per active-space size:\n"
            f"{details}"
        )

    return [selected[n_qubits][0] for n_qubits in sorted(requested)]


def dense_pauli_string(pauli_key: PauliKey, n_qubits: int) -> str:
    """Return an n-character Pauli string in qubit order 0, 1, ..., n-1."""
    operators = dict(pauli_key)
    return "".join(operators.get(qubit, "I") for qubit in range(n_qubits))


def signed_coefficient_lexicographic_order(
    raw_pauli_keys: Sequence[PauliKey],
    final_coefficients: dict[PauliKey, complex],
    n_qubits: int,
    tolerance: float,
) -> list[PauliKey]:
    """Sort by ascending signed coefficient, then dense Pauli string.

    The coefficient comparison uses the validated real coefficient directly.
    The lexicographic key is used only when coefficients are exactly equal as
    floating-point values, matching the clarified baseline wording.
    """
    return sorted(
        raw_pauli_keys,
        key=lambda key: (
            real_coefficient(final_coefficients[key], tolerance),
            dense_pauli_string(key, n_qubits),
        ),
    )


def hash_pauli_order(pauli_keys: Sequence[PauliKey], n_qubits: int) -> str:
    payload = "\n".join(
        dense_pauli_string(pauli_key, n_qubits)
        for pauli_key in pauli_keys
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def safe_ratio(numerator: float, denominator: float) -> str | float:
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        return ""
    if denominator <= 0.0:
        return ""
    return numerator / denominator


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


def ordering_definition(ordering_name: str) -> str:
    if ordering_name == "jw_raw":
        return "OpenFermion/QHAT final JW insertion order"
    if ordering_name == "signed_coefficient_lexicographic":
        return (
            "ascending signed coefficient; exact ties broken by ascending "
            "dense Pauli-string lexicographic order"
        )
    raise ValueError(f"Unsupported ordering {ordering_name!r}.")


def benchmark_case(
    tensor_path: Path,
    args: argparse.Namespace,
) -> list[dict[str, Any]]:
    interaction, n_qubits = load_interaction_operator(tensor_path)
    metadata = parse_case_metadata(tensor_path, n_qubits)

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
    signed_baseline_keys = signed_coefficient_lexicographic_order(
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        n_qubits=n_qubits,
        tolerance=args.tolerance,
    )

    validate_pauli_order(
        "signed_coefficient_lexicographic",
        signed_baseline_keys,
        raw_pauli_keys,
    )

    orderings = {
        "jw_raw": list(raw_pauli_keys),
        "signed_coefficient_lexicographic": signed_baseline_keys,
    }

    raw_index_by_key = {
        key: index for index, key in enumerate(raw_pauli_keys)
    }
    ordering_indices = {
        name: [raw_index_by_key[key] for key in keys]
        for name, keys in orderings.items()
    }

    print("    Building Pauli noncommutation graph and BCH evaluator...")
    jw_graph, _ = build_pauli_noncommutation_graph(raw_pauli_keys)
    bch_evaluator = HFCommutatorEvaluator.build(
        pauli_keys=raw_pauli_keys,
        coefficients=final_coefficients,
        pauli_graph=jw_graph,
        n_qubits=n_qubits,
        n_electrons=metadata.active_occupied,
        tolerance=args.tolerance,
    )

    print("    Building exact symmetry-sector reference...")
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

    rows: list[dict[str, Any]] = []
    raw_one_minus_overlap: float | None = None
    raw_infidelity: float | None = None
    raw_bch_norm: float | None = None

    for ordering_name in ORDERING_NAMES:
        pauli_keys = orderings[ordering_name]
        pauli_order_indices = ordering_indices[ordering_name]

        print(f"    Running {ordering_name}...")
        row = blank_row()
        add_metadata(row, metadata, tensor_path)
        row.update(
            {
                "status": "success",
                "number_of_pauli_terms": len(raw_pauli_keys),
                "ordering": ordering_name,
                "ordering_definition": ordering_definition(ordering_name),
                "pauli_order_hash": hash_pauli_order(pauli_keys, n_qubits),
                "first_pauli_string": dense_pauli_string(
                    pauli_keys[0], n_qubits
                ),
                "first_coefficient": real_coefficient(
                    final_coefficients[pauli_keys[0]], args.tolerance
                ),
                "last_pauli_string": dense_pauli_string(
                    pauli_keys[-1], n_qubits
                ),
                "last_coefficient": real_coefficient(
                    final_coefficients[pauli_keys[-1]], args.tolerance
                ),
                "trotter_steps": args.steps,
                "trotter_dt": args.evolution_time / args.steps,
                "evolution_time": args.evolution_time,
                "exact_sector_dimension": exact_sector_state.size,
                "exact_build_time_seconds": exact_build_time,
                "exact_evolution_time_seconds": exact_evolution_time,
                "coefficient_tolerance": args.tolerance,
            }
        )

        bch_norm = bch_evaluator.evaluate(pauli_order_indices)
        row["bch2_hf_state_norm"] = bch_norm

        ordered_terms = [
            compiled_raw_terms[index]
            for index in pauli_order_indices
        ]
        trotter_start = time.perf_counter()
        approximate_state, nominal_exponential_count = evolve_trotter_state(
            initial_state=initial_state,
            terms=ordered_terms,
            formula_order=1,
            trotter_steps=args.steps,
            evolution_time=args.evolution_time,
            parallel_threshold=args.parallel_threshold,
        )
        row["trotter_runtime_seconds"] = time.perf_counter() - trotter_start
        row["nominal_exponential_count"] = nominal_exponential_count

        metrics = compare_states(
            exact_sector_state=exact_sector_state,
            exact_basis_indices=exact_basis_indices,
            approximate_full_state=approximate_state,
            number_basis_indices=number_basis_indices,
        )
        row.update(metrics)

        overlap = float(metrics["state_overlap_abs"])
        one_minus_overlap = 1.0 - overlap
        infidelity = float(metrics["state_infidelity"])
        row["one_minus_overlap"] = one_minus_overlap

        if ordering_name == "jw_raw":
            raw_one_minus_overlap = one_minus_overlap
            raw_infidelity = infidelity
            raw_bch_norm = bch_norm
            row["one_minus_overlap_ratio_to_raw_jw"] = 1.0
            row["state_infidelity_ratio_to_raw_jw"] = 1.0
            row["bch_squared_ratio_to_raw_jw"] = 1.0
        else:
            if (
                raw_one_minus_overlap is None
                or raw_infidelity is None
                or raw_bch_norm is None
            ):
                raise RuntimeError("Raw JW must be evaluated first.")
            row["one_minus_overlap_ratio_to_raw_jw"] = safe_ratio(
                one_minus_overlap,
                raw_one_minus_overlap,
            )
            row["state_infidelity_ratio_to_raw_jw"] = safe_ratio(
                infidelity,
                raw_infidelity,
            )
            bch_ratio = safe_ratio(bch_norm, raw_bch_norm)
            row["bch_squared_ratio_to_raw_jw"] = (
                bch_ratio**2 if isinstance(bch_ratio, float) else ""
            )

        rows.append(row)
        print(
            "      "
            f"1-overlap={one_minus_overlap:.6e}  "
            f"infidelity={infidelity:.6e}  "
            f"BCH_HF={bch_norm:.6e}"
        )

    return rows


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    warm_up_numba()

    tensor_paths = select_cases(
        library=args.library,
        requested_qubits=args.qubits,
        bond_length=args.bond_length,
    )

    print("=" * 96)
    print("B2 signed-coefficient baseline study")
    print("=" * 96)
    print(f"Cases:             {len(tensor_paths)}")
    print(f"Active qubits:     {sorted(args.qubits)}")
    print(f"Bond length:       {args.bond_length}")
    print(f"Trotter steps:     {args.steps}")
    print(f"Evolution time:    {args.evolution_time}")
    print(f"Output:            {args.output}")
    print(
        "Baseline: ascending signed coefficient; exact coefficient ties "
        "use ascending dense Pauli-string lexicographic order."
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, Any]] = []

    for case_index, tensor_path in enumerate(tensor_paths, start=1):
        print()
        print("-" * 96)
        print(f"[{case_index}/{len(tensor_paths)}] {tensor_path}")
        try:
            all_rows.extend(benchmark_case(tensor_path, args))
        except Exception as error:
            print(f"    FAILED: {type(error).__name__}: {error}")
            failure = blank_row()
            failure.update(
                {
                    "status": "failed",
                    "error_message": f"{type(error).__name__}: {error}",
                    "case_id": tensor_path.name.removesuffix(".tensors.npz"),
                    "tensor_path": str(tensor_path),
                    "trotter_steps": args.steps,
                    "trotter_dt": args.evolution_time / args.steps,
                    "evolution_time": args.evolution_time,
                    "coefficient_tolerance": args.tolerance,
                }
            )
            all_rows.append(failure)

    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(all_rows)

    successful = sum(row["status"] == "success" for row in all_rows)
    failed = sum(row["status"] == "failed" for row in all_rows)

    print()
    print("=" * 96)
    print(f"Successful rows: {successful}")
    print(f"Failed rows:     {failed}")
    print(f"Results:         {args.output}")
    print("=" * 96)

def fermionic_term_order_indices(
    hermitian_terms,
    ordering_method: str,
    tolerance: float,
) -> list[int]:
    """
    Return indices of complete Hermitian fermionic terms.

    Supported methods:
      - signed_ascending:
          increasing signed coefficient, then fermionic lexicographic order
      - magnitude_descending:
          decreasing absolute coefficient, then fermionic lexicographic order
    """

    def canonical_component_key(term):
        # One deterministic representative of the Hermitian pair.
        return min(term.component_keys)

    def signed_weight(term):
        key = canonical_component_key(term)
        return real_coefficient(
            term.operator.terms[key],
            tolerance,
        )

    if ordering_method == "signed_ascending":
        return sorted(
            range(len(hermitian_terms)),
            key=lambda index: (
                signed_weight(hermitian_terms[index]),
                canonical_component_key(hermitian_terms[index]),
            ),
        )

    if ordering_method == "magnitude_descending":
        return sorted(
            range(len(hermitian_terms)),
            key=lambda index: (
                -abs(signed_weight(hermitian_terms[index])),
                canonical_component_key(hermitian_terms[index]),
            ),
        )

    raise ValueError(
        f"Unsupported fermionic ordering method: {ordering_method}"
    )

if __name__ == "__main__":
    main()