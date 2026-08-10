#!/usr/bin/env python3
"""Exact cross-time Trotter-error coherence study for Be2 across basis sets.

This benchmark is designed to test the mechanism suggested by the Be2
trajectory study:

    A schedule can have a larger *instantaneous* local Trotter defect at every
    time and still have a smaller final error because local defects accumulate
    coherently and can cancel after subsequent Trotter propagation.

For a fixed first-order Trotter slice S and exact one-step propagator U,

    S^N - U^N = sum_{j=0}^{N-1} S^(N-1-j) (S-U) U^j.

For the exact state psi_j = U^j psi_0, define the local defect

    delta_j = (S-U) psi_j.

The final state-vector error is therefore exactly

    Delta_N = sum_j S^(N-1-j) delta_j.

Because S is unitary,

    ||Delta_N||^2
      = sum_j ||delta_j||^2
        + 2 sum_{m<n} Re <v_m|v_n>,

where v_j = S^(N-1-j) delta_j.  The second term is the exact cross-time
interference term.  A strongly negative value means destructive cancellation
between defects generated at different Trotter steps.

The script also computes an exact telescoping decomposition of the final
complex overlap deviation.  If psi_N is the exact final state,

    <psi_N|psi_Trotter> - 1 = sum_j a_j,

with

    a_j = < (S^dagger)^(N-1-j) psi_N | delta_j >.

This gives two complementary cancellation diagnostics:

  1. vector_cancellation_ratio = ||Delta_N||^2 / sum_j ||delta_j||^2
  2. overlap_contribution_cancellation_ratio = |sum_j a_j| / sum_j |a_j|

Both are exact identities for the implemented product formula; no BCH
truncation is used in these two decompositions.

Recommended use: compare Be2/HGBS-5/20q and Be2/STO-6G/20q, especially
fermionic_signed_reference versus signed_parent_descendants_round_robin.

Run from the QHAT repository root after placing this file in analysis/.
The benchmark is append/resume safe at the completed (case, schedule, steps,
time) level.  Trace rows are only appended after a schedule finishes, so an
interrupted schedule is cleanly recomputed on restart.
"""

from __future__ import annotations

import argparse
import csv
import math
import time
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from openfermion import get_fermion_operator, jordan_wigner
try:
    from openfermion import get_number_preserving_sparse_operator
except ImportError:  # Compatibility with older OpenFermion exports.
    from openfermion.linalg import get_number_preserving_sparse_operator
from scipy.sparse.linalg import expm_multiply

try:
    from qhat.analysis import benchmark_fermionic_structure_ablation as ablation
    from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
    from qhat.analysis.benchmark_b2_active_spaces_matrix_free import (
        build_hartree_fock_state,
        compile_ordered_terms,
        evolve_trotter_state,
        remove_fermionic_identity,
        spin_sector_basis_indices,
        warm_up_numba,
    )
except ImportError:
    import benchmark_fermionic_structure_ablation as ablation
    import benchmark_b2_signed_coefficient_baseline as baseline
    from benchmark_b2_active_spaces_matrix_free import (
        build_hartree_fock_state,
        compile_ordered_terms,
        evolve_trotter_state,
        remove_fermionic_identity,
        spin_sector_basis_indices,
        warm_up_numba,
    )


REFERENCE = "fermionic_signed_reference"
FMAG = "fermionic_magnitude_reference"
JW = "jw_signed_baseline"
ROUND = "signed_parent_descendants_round_robin"
SCHEDULES = (REFERENCE, FMAG, JW, ROUND)

SUMMARY_FIELDS = [
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
    "schedule",
    "pauli_order_hash",
    "trotter_steps",
    "dt",
    "evolution_time",
    "nominal_exponential_count",
    "exact_sector_dimension",
    "actual_overlap_real",
    "actual_overlap_imag",
    "actual_overlap_abs",
    "actual_one_minus_overlap",
    "actual_infidelity",
    "actual_raw_state_error_norm2",
    "actual_phase_aligned_error_norm2",
    "sum_local_defect_norm2",
    "sum_local_phase_aligned_error_norm2",
    "vector_interference_term",
    "vector_interference_fraction",
    "vector_cancellation_ratio",
    "overlap_deviation_abs",
    "sum_abs_overlap_contributions",
    "overlap_contribution_cancellation_ratio",
    "net_real_overlap_loss",
    "sum_abs_real_overlap_loss_contributions",
    "real_overlap_loss_cancellation_ratio",
    "reconstructed_overlap_real",
    "reconstructed_overlap_imag",
    "reconstructed_overlap_abs",
    "overlap_reconstruction_abs_error",
    "backward_exact_recovery_2norm_error",
    "final_trotter_norm_error",
    "mean_local_defect_norm",
    "max_local_defect_norm",
    "min_local_defect_norm",
    "runtime_seconds",
    "coefficient_tolerance",
]

TRACE_FIELDS = [
    "case_id",
    "basis",
    "schedule",
    "trotter_steps",
    "dt",
    "evolution_time",
    "step_index",
    "step_start_time",
    "step_end_time",
    "future_trotter_slices",
    "local_defect_norm",
    "local_defect_norm2",
    "local_overlap_abs",
    "local_one_minus_overlap",
    "overlap_contribution_real",
    "overlap_contribution_imag",
    "overlap_contribution_abs",
    "signed_real_overlap_loss_contribution",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Measure exact cross-time cancellation of local first-order "
            "Trotter defects for Be2 ordering schedules."
        )
    )
    parser.add_argument(
        "--tensor",
        type=Path,
        action="append",
        required=True,
        help="Tensor file; repeat for HGBS-5 and STO-6G.",
    )
    parser.add_argument(
        "--steps",
        type=int,
        nargs="+",
        default=[50, 100, 200],
        help="Trotter step counts. Default: 50 100 200.",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=1.0,
        help="Total evolution time. Default: 1.0.",
    )
    parser.add_argument(
        "--schedules",
        nargs="+",
        choices=SCHEDULES,
        default=[REFERENCE, ROUND],
        help=(
            "Schedules to decompose. Default is the focused comparison: "
            "fermionic signed and round robin."
        ),
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=Path("analysis/be2_cross_time_coherence_summary.csv"),
    )
    parser.add_argument(
        "--trace-output",
        type=Path,
        default=Path("analysis/be2_cross_time_coherence_trace.csv"),
    )
    parser.add_argument(
        "--parallel-threshold",
        type=int,
        default=1 << 18,
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0e-12,
    )
    parser.add_argument(
        "--reconstruction-tolerance",
        type=float,
        default=5.0e-8,
        help=(
            "Maximum allowed |overlap_direct-overlap_telescoped|. "
            "Default: 5e-8."
        ),
    )
    return parser.parse_args()


def blank_summary() -> dict[str, Any]:
    return {field: "" for field in SUMMARY_FIELDS}


def completed_keys(path: Path) -> set[tuple[str, str, int, float]]:
    keys: set[tuple[str, str, int, float]] = set()
    if not path.exists() or path.stat().st_size == 0:
        return keys
    with path.open("r", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row.get("status") != "success":
                continue
            try:
                keys.add(
                    (
                        row["case_id"],
                        row["schedule"],
                        int(row["trotter_steps"]),
                        float(row["evolution_time"]),
                    )
                )
            except (KeyError, TypeError, ValueError):
                continue
    return keys


def safe_ratio(numerator: float, denominator: float) -> float | str:
    if not np.isfinite(denominator) or abs(denominator) <= 0.0:
        return ""
    return numerator / denominator


def build_deterministic_orders(
    raw_pauli_keys: Sequence[ablation.PauliKey],
    final_coefficients: dict[ablation.PauliKey, complex],
    n_qubits: int,
    fermion_hamiltonian: Any,
    tolerance: float,
) -> tuple[dict[str, tuple[int, ...]], int]:
    """Build the same four deterministic schedules used by the ablation."""
    raw_index_by_key = {key: i for i, key in enumerate(raw_pauli_keys)}
    hermitian_terms = ablation.build_hermitian_fermion_terms(
        fermion_hamiltonian,
        tolerance,
    )
    fermion_to_pauli_indices = ablation.precompute_fermion_to_pauli_indices(
        hermitian_terms=hermitian_terms,
        final_coefficients=final_coefficients,
        raw_index_by_key=raw_index_by_key,
        tolerance=tolerance,
    )
    signed_parent_order = baseline.fermionic_term_order_indices(
        hermitian_terms=hermitian_terms,
        ordering_method="signed_ascending",
        tolerance=tolerance,
    )
    magnitude_parent_order = baseline.fermionic_term_order_indices(
        hermitian_terms=hermitian_terms,
        ordering_method="magnitude_descending",
        tolerance=tolerance,
    )
    reference_buckets, _, fallback = ablation.build_reference_parent_buckets(
        signed_parent_order,
        fermion_to_pauli_indices,
        len(raw_pauli_keys),
    )
    reference = ablation.induced_pauli_order_indices(
        fermionic_node_order=signed_parent_order,
        fermion_to_pauli_indices=fermion_to_pauli_indices,
        number_of_pauli_terms=len(raw_pauli_keys),
    )
    magnitude = ablation.induced_pauli_order_indices(
        fermionic_node_order=magnitude_parent_order,
        fermion_to_pauli_indices=fermion_to_pauli_indices,
        number_of_pauli_terms=len(raw_pauli_keys),
    )
    jw_keys = baseline.signed_coefficient_lexicographic_order(
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        n_qubits=n_qubits,
        tolerance=tolerance,
    )
    jw = [raw_index_by_key[key] for key in jw_keys]
    round_robin, _ = ablation.round_robin_descendants(reference_buckets, fallback)

    orders = {
        REFERENCE: tuple(map(int, reference)),
        FMAG: tuple(map(int, magnitude)),
        JW: tuple(map(int, jw)),
        ROUND: tuple(map(int, round_robin)),
    }
    for name, indices in orders.items():
        ablation.validate_pauli_order(
            name,
            [raw_pauli_keys[index] for index in indices],
            raw_pauli_keys,
        )
    return orders, len(hermitian_terms)


def build_exact_spin_problem(
    fermion_hamiltonian: Any,
    n_qubits: int,
    n_electrons: int,
    tolerance: float,
) -> tuple[Any, np.ndarray, np.ndarray]:
    reference_determinant = np.array(
        [index < n_electrons for index in range(n_qubits)],
        dtype=bool,
    )
    sparse_hamiltonian = get_number_preserving_sparse_operator(
        remove_fermionic_identity(fermion_hamiltonian, tolerance),
        num_qubits=n_qubits,
        num_electrons=n_electrons,
        spin_preserving=True,
        reference_determinant=reference_determinant,
        excitation_level=None,
    )
    basis_indices = spin_sector_basis_indices(n_qubits, n_electrons)
    if sparse_hamiltonian.shape[0] != basis_indices.size:
        raise ValueError(
            "Exact spin-sector dimension mismatch: "
            f"{sparse_hamiltonian.shape[0]} != {basis_indices.size}."
        )
    initial_sector = np.zeros(sparse_hamiltonian.shape[0], dtype=np.complex128)
    initial_sector[0] = 1.0
    return sparse_hamiltonian, basis_indices, initial_sector


def scatter_sector_state(
    sector_state: np.ndarray,
    basis_indices: np.ndarray,
    full_dimension: int,
) -> np.ndarray:
    full = np.zeros(full_dimension, dtype=np.complex128)
    full[basis_indices] = sector_state
    return full


def apply_slice(
    state: np.ndarray,
    terms: Sequence[tuple[float, np.uint64, np.uint64, int]],
    dt: float,
    parallel_threshold: int,
) -> np.ndarray:
    result, _ = evolve_trotter_state(
        initial_state=state,
        terms=terms,
        formula_order=1,
        trotter_steps=1,
        evolution_time=dt,
        parallel_threshold=parallel_threshold,
    )
    return result


def direct_final_metrics(
    exact_final_sector: np.ndarray,
    spin_basis_indices: np.ndarray,
    trotter_final_full: np.ndarray,
) -> tuple[complex, float, float, float, float, float, float]:
    overlap = np.vdot(exact_final_sector, trotter_final_full[spin_basis_indices])
    overlap_abs = min(1.0, max(0.0, float(abs(overlap))))
    one_minus = max(0.0, 1.0 - overlap_abs)
    infidelity = max(0.0, 1.0 - overlap_abs**2)

    exact_final_full = scatter_sector_state(
        exact_final_sector,
        spin_basis_indices,
        trotter_final_full.size,
    )
    raw_error = trotter_final_full - exact_final_full
    raw_norm2 = float(np.vdot(raw_error, raw_error).real)

    if abs(overlap) > 0.0:
        phase = np.exp(-1j * np.angle(overlap))
    else:
        phase = 1.0 + 0.0j
    aligned_error = phase * trotter_final_full - exact_final_full
    aligned_norm2 = float(np.vdot(aligned_error, aligned_error).real)
    norm_error = abs(float(np.linalg.norm(trotter_final_full)) - 1.0)
    return (
        overlap,
        overlap_abs,
        one_minus,
        infidelity,
        raw_norm2,
        aligned_norm2,
        norm_error,
    )


def run_one_schedule(
    *,
    metadata: Any,
    tensor_path: Path,
    schedule: str,
    terms: Sequence[tuple[float, np.uint64, np.uint64, int]],
    order_hash: str,
    sparse_hamiltonian: Any,
    spin_basis_indices: np.ndarray,
    initial_sector: np.ndarray,
    initial_full: np.ndarray,
    number_of_fermionic_terms: int,
    number_of_pauli_terms: int,
    n_qubits: int,
    steps: int,
    evolution_time: float,
    parallel_threshold: int,
    tolerance: float,
    reconstruction_tolerance: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    run_start = time.perf_counter()
    dt = evolution_time / steps
    inverse_terms = tuple(reversed(terms))

    # Exact final state in the fixed-N, fixed-Sz sector.
    exact_final_sector = np.asarray(
        expm_multiply((-1j * evolution_time) * sparse_hamiltonian, initial_sector),
        dtype=np.complex128,
    )
    exact_final_full = scatter_sector_state(
        exact_final_sector,
        spin_basis_indices,
        initial_full.size,
    )

    # Direct product-formula result.  This is the final quantity being explained.
    trotter_final_full, nominal_count = evolve_trotter_state(
        initial_state=initial_full,
        terms=terms,
        formula_order=1,
        trotter_steps=steps,
        evolution_time=evolution_time,
        parallel_threshold=parallel_threshold,
    )
    (
        direct_overlap,
        direct_overlap_abs,
        direct_one_minus,
        direct_infidelity,
        raw_final_norm2,
        aligned_final_norm2,
        final_norm_error,
    ) = direct_final_metrics(
        exact_final_sector,
        spin_basis_indices,
        trotter_final_full,
    )

    # Backward telescoping pass.
    # bra_full at step j equals (S^dagger)^(N-1-j) |psi_exact(T)>.
    bra_full = exact_final_full.copy()
    exact_next_sector = exact_final_sector.copy()
    contributions: list[complex] = []
    trace_rows: list[dict[str, Any]] = []
    local_norms: list[float] = []
    sum_local_defect_norm2 = 0.0
    sum_local_phase_aligned_norm2 = 0.0

    for j in range(steps - 1, -1, -1):
        exact_current_sector = np.asarray(
            expm_multiply((+1j * dt) * sparse_hamiltonian, exact_next_sector),
            dtype=np.complex128,
        )
        exact_current_full = scatter_sector_state(
            exact_current_sector,
            spin_basis_indices,
            initial_full.size,
        )
        exact_next_full = scatter_sector_state(
            exact_next_sector,
            spin_basis_indices,
            initial_full.size,
        )

        trotter_from_exact = apply_slice(
            exact_current_full,
            terms,
            dt,
            parallel_threshold,
        )
        local_defect = trotter_from_exact - exact_next_full
        local_norm2 = float(np.vdot(local_defect, local_defect).real)
        local_norm = math.sqrt(max(0.0, local_norm2))
        local_norms.append(local_norm)
        sum_local_defect_norm2 += local_norm2

        local_overlap = np.vdot(
            exact_next_sector,
            trotter_from_exact[spin_basis_indices],
        )
        local_overlap_abs = min(1.0, max(0.0, float(abs(local_overlap))))
        local_one_minus = max(0.0, 1.0 - local_overlap_abs)
        # For normalized states, the phase-aligned squared 2-norm error is
        # exactly 2(1-|overlap|), up to floating-point norm drift.
        sum_local_phase_aligned_norm2 += 2.0 * local_one_minus

        contribution = np.vdot(bra_full, local_defect)
        contributions.append(complex(contribution))

        trace_rows.append(
            {
                "case_id": metadata.case_id,
                "basis": metadata.basis,
                "schedule": schedule,
                "trotter_steps": steps,
                "dt": dt,
                "evolution_time": evolution_time,
                "step_index": j,
                "step_start_time": j * dt,
                "step_end_time": (j + 1) * dt,
                "future_trotter_slices": steps - 1 - j,
                "local_defect_norm": local_norm,
                "local_defect_norm2": local_norm2,
                "local_overlap_abs": local_overlap_abs,
                "local_one_minus_overlap": local_one_minus,
                "overlap_contribution_real": float(contribution.real),
                "overlap_contribution_imag": float(contribution.imag),
                "overlap_contribution_abs": float(abs(contribution)),
                "signed_real_overlap_loss_contribution": float(-contribution.real),
            }
        )

        # Move the final exact bra one Trotter slice backward: S^dagger.
        # S^dagger is the reversed product with -dt.
        bra_full = apply_slice(
            bra_full,
            inverse_terms,
            -dt,
            parallel_threshold,
        )
        exact_next_sector = exact_current_sector

    # Restore natural time order in the trace CSV.
    trace_rows.sort(key=lambda row: int(row["step_index"]))

    contribution_sum = sum(contributions, 0.0 + 0.0j)
    reconstructed_overlap = 1.0 + contribution_sum
    reconstruction_error = float(abs(reconstructed_overlap - direct_overlap))
    if reconstruction_error > reconstruction_tolerance:
        raise RuntimeError(
            "Telescoping overlap reconstruction failed: "
            f"|reconstructed-direct|={reconstruction_error:.3e} > "
            f"{reconstruction_tolerance:.3e}.  This usually indicates an "
            "incorrect inverse product ordering or a convention mismatch."
        )

    backward_recovery = float(np.linalg.norm(exact_next_sector - initial_sector))
    vector_interference = raw_final_norm2 - sum_local_defect_norm2
    vector_interference_fraction = safe_ratio(
        vector_interference,
        sum_local_defect_norm2,
    )
    vector_cancellation_ratio = safe_ratio(
        raw_final_norm2,
        sum_local_defect_norm2,
    )

    sum_abs_contrib = float(sum(abs(value) for value in contributions))
    overlap_cancel_ratio = safe_ratio(abs(contribution_sum), sum_abs_contrib)
    net_real_loss = float(-contribution_sum.real)
    sum_abs_real_loss = float(sum(abs(value.real) for value in contributions))
    real_loss_cancel_ratio = safe_ratio(abs(contribution_sum.real), sum_abs_real_loss)

    summary = blank_summary()
    summary.update(
        {
            "status": "success",
            "case_id": metadata.case_id,
            "tensor_path": str(tensor_path),
            "molecule": metadata.molecule,
            "bond_length": metadata.bond_length,
            "basis": metadata.basis,
            "active_occupied": metadata.active_occupied,
            "active_vacant": metadata.active_vacant,
            "n_qubits": n_qubits,
            "number_of_fermionic_terms": number_of_fermionic_terms,
            "number_of_pauli_terms": number_of_pauli_terms,
            "schedule": schedule,
            "pauli_order_hash": order_hash,
            "trotter_steps": steps,
            "dt": dt,
            "evolution_time": evolution_time,
            "nominal_exponential_count": nominal_count,
            "exact_sector_dimension": initial_sector.size,
            "actual_overlap_real": float(direct_overlap.real),
            "actual_overlap_imag": float(direct_overlap.imag),
            "actual_overlap_abs": direct_overlap_abs,
            "actual_one_minus_overlap": direct_one_minus,
            "actual_infidelity": direct_infidelity,
            "actual_raw_state_error_norm2": raw_final_norm2,
            "actual_phase_aligned_error_norm2": aligned_final_norm2,
            "sum_local_defect_norm2": sum_local_defect_norm2,
            "sum_local_phase_aligned_error_norm2": sum_local_phase_aligned_norm2,
            "vector_interference_term": vector_interference,
            "vector_interference_fraction": vector_interference_fraction,
            "vector_cancellation_ratio": vector_cancellation_ratio,
            "overlap_deviation_abs": float(abs(contribution_sum)),
            "sum_abs_overlap_contributions": sum_abs_contrib,
            "overlap_contribution_cancellation_ratio": overlap_cancel_ratio,
            "net_real_overlap_loss": net_real_loss,
            "sum_abs_real_overlap_loss_contributions": sum_abs_real_loss,
            "real_overlap_loss_cancellation_ratio": real_loss_cancel_ratio,
            "reconstructed_overlap_real": float(reconstructed_overlap.real),
            "reconstructed_overlap_imag": float(reconstructed_overlap.imag),
            "reconstructed_overlap_abs": float(abs(reconstructed_overlap)),
            "overlap_reconstruction_abs_error": reconstruction_error,
            "backward_exact_recovery_2norm_error": backward_recovery,
            "final_trotter_norm_error": final_norm_error,
            "mean_local_defect_norm": float(np.mean(local_norms)),
            "max_local_defect_norm": float(np.max(local_norms)),
            "min_local_defect_norm": float(np.min(local_norms)),
            "runtime_seconds": time.perf_counter() - run_start,
            "coefficient_tolerance": tolerance,
        }
    )
    return summary, trace_rows


def append_summary_row(
    writer: csv.DictWriter,
    stream: Any,
    row: dict[str, Any],
) -> None:
    writer.writerow({field: row.get(field, "") for field in SUMMARY_FIELDS})
    stream.flush()


def append_trace_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=TRACE_FIELDS)
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in TRACE_FIELDS})
        stream.flush()


def main() -> int:
    args = parse_args()
    if args.time <= 0.0:
        raise ValueError("--time must be positive.")
    steps_values = sorted(set(int(value) for value in args.steps))
    if any(value <= 0 for value in steps_values):
        raise ValueError("All --steps values must be positive integers.")

    warm_up_numba()

    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    completed = completed_keys(args.summary_output)
    write_header = (
        not args.summary_output.exists() or args.summary_output.stat().st_size == 0
    )

    with args.summary_output.open("a", newline="", encoding="utf-8") as summary_stream:
        summary_writer = csv.DictWriter(summary_stream, fieldnames=SUMMARY_FIELDS)
        if write_header:
            summary_writer.writeheader()
            summary_stream.flush()

        for tensor_path in args.tensor:
            interaction, n_qubits = ablation.load_interaction_operator(tensor_path)
            metadata = ablation.parse_case_metadata(tensor_path, n_qubits)
            fermion_hamiltonian = ablation.clean_fermion_operator(
                get_fermion_operator(interaction),
                args.tolerance,
            )
            jw_hamiltonian = jordan_wigner(fermion_hamiltonian)
            jw_hamiltonian.compress(abs_tol=args.tolerance)
            final_coefficients = {
                key: coefficient
                for key, coefficient in jw_hamiltonian.terms.items()
                if key != () and abs(coefficient) > args.tolerance
            }
            raw_pauli_keys = list(final_coefficients)
            if not raw_pauli_keys:
                raise ValueError(f"No identity-free JW Pauli terms for {tensor_path}.")

            print("=" * 72)
            print(f"Case:   {metadata.case_id}")
            print(f"Basis:  {metadata.basis}")
            print(f"Tensor: {tensor_path}")
            print(f"Qubits: {n_qubits}")
            print("=" * 72)

            orders, number_of_fermionic_terms = build_deterministic_orders(
                raw_pauli_keys,
                final_coefficients,
                n_qubits,
                fermion_hamiltonian,
                args.tolerance,
            )
            compiled_raw_terms = compile_ordered_terms(
                raw_pauli_keys,
                final_coefficients,
                n_qubits,
                args.tolerance,
            )
            ordered_terms = {
                name: tuple(compiled_raw_terms[index] for index in orders[name])
                for name in args.schedules
            }

            print("Building exact fixed-N/fixed-Sz Hamiltonian...")
            sparse_hamiltonian, spin_basis_indices, initial_sector = (
                build_exact_spin_problem(
                    fermion_hamiltonian,
                    n_qubits,
                    metadata.active_occupied,
                    args.tolerance,
                )
            )
            initial_full = build_hartree_fock_state(
                n_qubits,
                metadata.active_occupied,
            )
            print(f"Exact sector dimension: {initial_sector.size}")
            print(f"Pauli terms: {len(raw_pauli_keys)}")

            for steps in steps_values:
                for schedule in args.schedules:
                    key = (metadata.case_id, schedule, steps, float(args.time))
                    if key in completed:
                        print(
                            f"SKIP {metadata.basis} steps={steps} "
                            f"schedule={schedule}"
                        )
                        continue

                    print(
                        f"RUN  {metadata.basis} steps={steps} "
                        f"dt={args.time/steps:.6g} schedule={schedule}"
                    )
                    try:
                        summary, trace = run_one_schedule(
                            metadata=metadata,
                            tensor_path=tensor_path,
                            schedule=schedule,
                            terms=ordered_terms[schedule],
                            order_hash=ablation.hash_integer_order(orders[schedule]),
                            sparse_hamiltonian=sparse_hamiltonian,
                            spin_basis_indices=spin_basis_indices,
                            initial_sector=initial_sector,
                            initial_full=initial_full,
                            number_of_fermionic_terms=number_of_fermionic_terms,
                            number_of_pauli_terms=len(raw_pauli_keys),
                            n_qubits=n_qubits,
                            steps=steps,
                            evolution_time=float(args.time),
                            parallel_threshold=args.parallel_threshold,
                            tolerance=args.tolerance,
                            reconstruction_tolerance=args.reconstruction_tolerance,
                        )
                    except Exception as exc:
                        failed = blank_summary()
                        failed.update(
                            {
                                "status": "failed",
                                "error_message": f"{type(exc).__name__}: {exc}",
                                "case_id": metadata.case_id,
                                "tensor_path": str(tensor_path),
                                "molecule": metadata.molecule,
                                "bond_length": metadata.bond_length,
                                "basis": metadata.basis,
                                "active_occupied": metadata.active_occupied,
                                "active_vacant": metadata.active_vacant,
                                "n_qubits": n_qubits,
                                "schedule": schedule,
                                "trotter_steps": steps,
                                "dt": args.time / steps,
                                "evolution_time": args.time,
                                "coefficient_tolerance": args.tolerance,
                            }
                        )
                        append_summary_row(summary_writer, summary_stream, failed)
                        raise

                    # Commit trace only after the entire schedule passes the exact
                    # telescoping reconstruction check.
                    append_trace_rows(args.trace_output, trace)
                    append_summary_row(summary_writer, summary_stream, summary)
                    completed.add(key)
                    print(
                        "  error={:.3e}  D={:.3e}  Q/D={:.3e}  "
                        "overlap-cancel={:.3e}  recon={:.3e}".format(
                            float(summary["actual_one_minus_overlap"]),
                            float(summary["sum_local_defect_norm2"]),
                            float(summary["vector_cancellation_ratio"]),
                            float(summary["overlap_contribution_cancellation_ratio"]),
                            float(summary["overlap_reconstruction_abs_error"]),
                        )
                    )

    print()
    print(f"Summary: {args.summary_output}")
    print(f"Trace:   {args.trace_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
