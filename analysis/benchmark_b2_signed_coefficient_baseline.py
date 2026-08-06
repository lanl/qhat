#!/usr/bin/env python3
"""Evaluate deterministic JW- and fermionic-informed Trotter orderings.

The benchmark evaluates five orderings of the same identity-free final
Jordan--Wigner Hamiltonian:

    jw_raw
        Raw OpenFermion Jordan--Wigner insertion order.

    signed_coefficient_lexicographic
        Final JW Pauli terms sorted by increasing signed coefficient, with
        ascending dense Pauli-string lexicographic order used to break exact
        coefficient ties.  This is the existing signed-coefficient baseline.

    jw_magnitude_descending_lexicographic
        Final JW Pauli terms sorted by decreasing absolute coefficient
        magnitude, with ascending dense Pauli-string lexicographic tie-breaking.

    fermionic_signed_coefficient_lexicographic
        Complete Hermitian fermionic terms sorted by increasing signed
        canonical coefficient, with fermionic lexicographic tie-breaking.  The
        ordered fermionic terms induce an ordering of the final combined JW
        Pauli factors using the existing first-occurrence rule.

    fermionic_magnitude_descending_lexicographic
        Complete Hermitian fermionic terms sorted by decreasing absolute
        canonical coefficient magnitude, with fermionic lexicographic
        tie-breaking.  The ordered fermionic terms induce an ordering of the
        final combined JW Pauli factors using the same first-occurrence rule.

The standalone command-line interface selects B2/STO-6G cases.  The generic
``benchmark_case`` function is also imported by
``benchmark_comparable_diatomics.py`` and can therefore be used for the other
molecules and bases in the L-sweep study.
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
        induced_pauli_order_indices,
        precompute_fermion_to_pauli_indices,
    )
except ImportError:
    from benchmark_b2_coloring_robustness import (
        HFCommutatorEvaluator,
        induced_pauli_order_indices,
        precompute_fermion_to_pauli_indices,
    )


PauliKey = tuple[tuple[int, str], ...]
FermionKey = tuple[tuple[int, int], ...]

ORDERING_NAMES = (
    "jw_raw",
    "signed_coefficient_lexicographic",
    "jw_magnitude_descending_lexicographic",
    "fermionic_signed_coefficient_lexicographic",
    "fermionic_magnitude_descending_lexicographic",
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
            "Evaluate raw JW, JW coefficient orderings, and fermionic-induced "
            "coefficient orderings for selected B2/STO-6G active spaces."
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
        default=Path("analysis/b2_deterministic_ordering_results.csv"),
        help="Output CSV for the deterministic ordering calculations.",
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
    """Sort final JW Pauli terms by increasing signed coefficient."""

    return sorted(
        raw_pauli_keys,
        key=lambda key: (
            real_coefficient(final_coefficients[key], tolerance),
            dense_pauli_string(key, n_qubits),
        ),
    )


def magnitude_descending_lexicographic_order(
    raw_pauli_keys: Sequence[PauliKey],
    final_coefficients: dict[PauliKey, complex],
    n_qubits: int,
    tolerance: float,
) -> list[PauliKey]:
    """Sort final JW Pauli terms by decreasing absolute coefficient."""

    return sorted(
        raw_pauli_keys,
        key=lambda key: (
            -abs(real_coefficient(final_coefficients[key], tolerance)),
            dense_pauli_string(key, n_qubits),
        ),
    )


def fermionic_term_lexicographic_key(
    term: HermitianFermionTerm,
) -> tuple[FermionKey, ...]:
    """Return a deterministic lexicographic key for one Hermitian term."""

    return tuple(sorted(term.component_keys))


def fermionic_term_signed_weight(
    term: HermitianFermionTerm,
    tolerance: float,
) -> float:
    """Return the signed coefficient of the canonical fermionic component.

    A complete Hermitian fermionic term may contain a monomial and its adjoint.
    The lexicographically smallest component is used as the deterministic
    representative.  The two components have the same absolute magnitude, so
    this convention is unambiguous for magnitude ordering and deterministic for
    signed ordering.
    """

    canonical_key = fermionic_term_lexicographic_key(term)[0]
    return real_coefficient(term.operator.terms[canonical_key], tolerance)


def fermionic_term_order_indices(
    hermitian_terms: Sequence[HermitianFermionTerm],
    ordering_method: str,
    tolerance: float,
) -> list[int]:
    """Return an ordering of complete Hermitian fermionic-term indices."""

    if ordering_method == "signed_ascending":
        return sorted(
            range(len(hermitian_terms)),
            key=lambda index: (
                fermionic_term_signed_weight(
                    hermitian_terms[index],
                    tolerance,
                ),
                fermionic_term_lexicographic_key(hermitian_terms[index]),
            ),
        )

    if ordering_method == "magnitude_descending":
        return sorted(
            range(len(hermitian_terms)),
            key=lambda index: (
                -abs(
                    fermionic_term_signed_weight(
                        hermitian_terms[index],
                        tolerance,
                    )
                ),
                fermionic_term_lexicographic_key(hermitian_terms[index]),
            ),
        )

    raise ValueError(
        f"Unsupported fermionic ordering method: {ordering_method!r}"
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


def add_metadata(
    row: dict[str, Any],
    metadata: Any,
    tensor_path: Path,
) -> None:
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
            "final JW Pauli terms ordered by ascending signed coefficient; "
            "exact ties broken by ascending dense Pauli-string "
            "lexicographic order"
        )

    if ordering_name == "jw_magnitude_descending_lexicographic":
        return (
            "final JW Pauli terms ordered by decreasing absolute coefficient "
            "magnitude; exact ties broken by ascending dense Pauli-string "
            "lexicographic order"
        )

    if ordering_name == "fermionic_signed_coefficient_lexicographic":
        return (
            "complete Hermitian fermionic terms ordered by ascending signed "
            "canonical coefficient, with fermionic lexicographic "
            "tie-breaking; final JW Pauli order induced by first occurrence"
        )

    if ordering_name == "fermionic_magnitude_descending_lexicographic":
        return (
            "complete Hermitian fermionic terms ordered by decreasing "
            "absolute canonical coefficient magnitude, with fermionic "
            "lexicographic tie-breaking; final JW Pauli order induced by "
            "first occurrence"
        )

    raise ValueError(f"Unsupported ordering {ordering_name!r}.")


def build_deterministic_orderings(
    fermion_hamiltonian: Any,
    raw_pauli_keys: Sequence[PauliKey],
    final_coefficients: dict[PauliKey, complex],
    n_qubits: int,
    tolerance: float,
) -> dict[str, list[PauliKey]]:
    """Build and validate all five deterministic Pauli-factor orders."""

    raw_pauli_keys = list(raw_pauli_keys)
    raw_index_by_key = {
        key: index
        for index, key in enumerate(raw_pauli_keys)
    }

    jw_signed_keys = signed_coefficient_lexicographic_order(
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        n_qubits=n_qubits,
        tolerance=tolerance,
    )
    jw_magnitude_keys = magnitude_descending_lexicographic_order(
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        n_qubits=n_qubits,
        tolerance=tolerance,
    )

    hermitian_terms = build_hermitian_fermion_terms(
        fermion_hamiltonian,
        tolerance,
    )
    fermion_to_pauli_indices = precompute_fermion_to_pauli_indices(
        hermitian_terms=hermitian_terms,
        final_coefficients=final_coefficients,
        raw_index_by_key=raw_index_by_key,
        tolerance=tolerance,
    )

    fermionic_signed_nodes = fermionic_term_order_indices(
        hermitian_terms=hermitian_terms,
        ordering_method="signed_ascending",
        tolerance=tolerance,
    )
    fermionic_magnitude_nodes = fermionic_term_order_indices(
        hermitian_terms=hermitian_terms,
        ordering_method="magnitude_descending",
        tolerance=tolerance,
    )

    fermionic_signed_indices = induced_pauli_order_indices(
        fermionic_node_order=fermionic_signed_nodes,
        fermion_to_pauli_indices=fermion_to_pauli_indices,
        number_of_pauli_terms=len(raw_pauli_keys),
    )
    fermionic_magnitude_indices = induced_pauli_order_indices(
        fermionic_node_order=fermionic_magnitude_nodes,
        fermion_to_pauli_indices=fermion_to_pauli_indices,
        number_of_pauli_terms=len(raw_pauli_keys),
    )

    orderings = {
        "jw_raw": list(raw_pauli_keys),
        "signed_coefficient_lexicographic": jw_signed_keys,
        "jw_magnitude_descending_lexicographic": jw_magnitude_keys,
        "fermionic_signed_coefficient_lexicographic": [
            raw_pauli_keys[index]
            for index in fermionic_signed_indices
        ],
        "fermionic_magnitude_descending_lexicographic": [
            raw_pauli_keys[index]
            for index in fermionic_magnitude_indices
        ],
    }

    if set(orderings) != set(ORDERING_NAMES):
        raise RuntimeError(
            "Ordering dictionary does not match ORDERING_NAMES."
        )

    for ordering_name, pauli_keys in orderings.items():
        validate_pauli_order(
            ordering_name,
            pauli_keys,
            raw_pauli_keys,
        )

    return orderings


def benchmark_case(
    tensor_path: Path,
    args: argparse.Namespace,
    ordering_names: Sequence[str] | None = None,
    raw_reference: dict[str, float] | None = None,
) -> list[dict[str, Any]]:
    """Evaluate selected deterministic orderings for one tensor case.

    When ``ordering_names`` is omitted, all orderings in ``ORDERING_NAMES``
    are evaluated, preserving the original behavior.  A missing-only driver
    can pass only the new ordering names so raw JW and the existing signed
    baseline are not evolved again.

    ``raw_reference`` may contain ``one_minus_overlap``, ``state_infidelity``,
    and ``bch2_hf_state_norm`` from an existing raw-JW row.  These values are
    used only to populate ratio columns; the selected ordering results do not
    depend on them.
    """
    selected_orderings = list(
        ORDERING_NAMES if ordering_names is None else ordering_names
    )
    if not selected_orderings:
        raise ValueError("At least one ordering must be selected.")
    if len(selected_orderings) != len(set(selected_orderings)):
        raise ValueError("Selected ordering names must be unique.")
    unknown_orderings = set(selected_orderings).difference(ORDERING_NAMES)
    if unknown_orderings:
        raise ValueError(
            "Unsupported ordering names: "
            + ", ".join(sorted(unknown_orderings))
        )

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
    if not raw_pauli_keys:
        raise ValueError("The identity-free JW Hamiltonian has no Pauli terms.")

    print("    Building deterministic JW and fermionic-induced orders...")
    orderings = build_deterministic_orderings(
        fermion_hamiltonian=fermion_hamiltonian,
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        n_qubits=n_qubits,
        tolerance=args.tolerance,
    )

    raw_index_by_key = {
        key: index
        for index, key in enumerate(raw_pauli_keys)
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
    raw_reference = raw_reference or {}

    def optional_reference_value(name: str) -> float | None:
        value = raw_reference.get(name)
        if value is None:
            return None
        try:
            number = float(value)
        except (TypeError, ValueError):
            return None
        return number if math.isfinite(number) else None

    raw_one_minus_overlap = optional_reference_value(
        "one_minus_overlap"
    )
    raw_infidelity = optional_reference_value(
        "state_infidelity"
    )
    raw_bch_norm = optional_reference_value(
        "bch2_hf_state_norm"
    )

    for ordering_name in selected_orderings:
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
                    pauli_keys[0],
                    n_qubits,
                ),
                "first_coefficient": real_coefficient(
                    final_coefficients[pauli_keys[0]],
                    args.tolerance,
                ),
                "last_pauli_string": dense_pauli_string(
                    pauli_keys[-1],
                    n_qubits,
                ),
                "last_coefficient": real_coefficient(
                    final_coefficients[pauli_keys[-1]],
                    args.tolerance,
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
            if raw_one_minus_overlap is not None:
                row["one_minus_overlap_ratio_to_raw_jw"] = safe_ratio(
                    one_minus_overlap,
                    raw_one_minus_overlap,
                )
            if raw_infidelity is not None:
                row["state_infidelity_ratio_to_raw_jw"] = safe_ratio(
                    infidelity,
                    raw_infidelity,
                )
            if raw_bch_norm is not None:
                bch_ratio = safe_ratio(bch_norm, raw_bch_norm)
                row["bch_squared_ratio_to_raw_jw"] = (
                    bch_ratio**2
                    if isinstance(bch_ratio, float)
                    else ""
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
    print("B2 deterministic JW and fermionic coefficient-ordering study")
    print("=" * 96)
    print(f"Cases:             {len(tensor_paths)}")
    print(f"Active qubits:     {sorted(args.qubits)}")
    print(f"Bond length:       {args.bond_length}")
    print(f"Trotter steps:     {args.steps}")
    print(f"Evolution time:    {args.evolution_time}")
    print(f"Output:            {args.output}")
    print("Orderings:")
    for ordering_name in ORDERING_NAMES:
        print(f"  - {ordering_name}")

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


if __name__ == "__main__":
    main()
