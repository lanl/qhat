#!/usr/bin/env python3
"""Trajectory-aware local Trotter-error case study for the Be2 crossover.

Purpose
-------
The existing HGBS-5 Be2/20q fixed-dt study shows that signed fermionic ordering
is best at short total time, while round-robin is better by t=1.  This script
asks *why* by measuring the one-step local Trotter error on the exact evolving
state |psi_exact(t)>.

For every checkpoint t and deterministic schedule pi, it compares

    exp(-i H dt) |psi_exact(t)>

with

    S_pi(dt) |psi_exact(t)>,

where S_pi is one first-order Trotter slice using the same final identity-free
JW Pauli Hamiltonian.  The primary diagnostic is the phase-aligned state-vector
error divided by dt^2.  The quantity

    2 * ||delta psi_local|| / dt^2

is an empirical trajectory-aware proxy for the leading BCH action norm.  At
small dt and t=0 it should be close to the existing HF-state BCH diagnostic,
but here it can be followed along the exact trajectory without explicitly
constructing a giant commutator operator.

The script is append-only: successful checkpoint/schedule rows are skipped.
Run it from the QHAT repository root after placing it in analysis/ next to
benchmark_fermionic_structure_ablation.py.
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
        compile_ordered_terms,
        evolve_trotter_state,
        number_sector_basis_indices,
        remove_fermionic_identity,
        spin_sector_basis_indices,
        warm_up_numba,
    )
except ImportError:
    import benchmark_fermionic_structure_ablation as ablation
    import benchmark_b2_signed_coefficient_baseline as baseline
    from benchmark_b2_active_spaces_matrix_free import (
        compile_ordered_terms,
        evolve_trotter_state,
        number_sector_basis_indices,
        remove_fermionic_identity,
        spin_sector_basis_indices,
        warm_up_numba,
    )


REFERENCE = "fermionic_signed_reference"
FMAG = "fermionic_magnitude_reference"
JW = "jw_signed_baseline"
ROUND = "signed_parent_descendants_round_robin"
SCHEDULES = (REFERENCE, FMAG, JW, ROUND)

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
    "schedule",
    "checkpoint_time",
    "local_dt",
    "pauli_order_hash",
    "exact_sector_dimension",
    "exact_checkpoint_advance_seconds",
    "exact_local_target_seconds",
    "local_trotter_seconds",
    "local_overlap_abs",
    "local_one_minus_overlap",
    "local_infidelity",
    "local_phase_aligned_2norm_error",
    "local_phase_aligned_2norm_error_over_dt2",
    "trajectory_bch_proxy_2err_over_dt2",
    "trajectory_bch_proxy_ratio_to_signed",
    "spin_sector_leakage_after_one_step",
    "particle_number_leakage_after_one_step",
    "coefficient_tolerance",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Measure one-step local first-order Trotter error on the exact "
            "evolving state to explain the Be2 schedule crossover."
        )
    )
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument(
        "--checkpoints",
        type=float,
        nargs="+",
        default=[0.0, 0.1, 0.25, 0.5, 0.6, 0.7, 0.75, 0.8, 1.0],
        help="Exact-trajectory checkpoint times. Default targets the Be2 crossover.",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=0.01,
        help="Single local Trotter-step duration. Default: 0.01.",
    )
    parser.add_argument(
        "--schedules",
        nargs="+",
        choices=SCHEDULES,
        default=list(SCHEDULES),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/be2_local_error_trajectory_hgbs5.csv"),
    )
    parser.add_argument(
        "--parallel-threshold",
        type=int,
        default=1 << 18,
        help="State-vector size threshold for parallel Pauli kernels.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0e-12,
    )
    parser.add_argument(
        "--no-spin-sector",
        action="store_true",
        help="Use only fixed particle number for exact evolution.",
    )
    return parser.parse_args()


def blank_row() -> dict[str, Any]:
    return {field: "" for field in FIELDNAMES}


def completed_keys(path: Path) -> set[tuple[str, str, float, float]]:
    keys: set[tuple[str, str, float, float]] = set()
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
                        float(row["checkpoint_time"]),
                        float(row["local_dt"]),
                    )
                )
            except (KeyError, TypeError, ValueError):
                continue
    return keys


def read_reference_proxies(path: Path) -> dict[tuple[str, float, float], float]:
    refs: dict[tuple[str, float, float], float] = {}
    if not path.exists() or path.stat().st_size == 0:
        return refs
    with path.open("r", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row.get("status") != "success" or row.get("schedule") != REFERENCE:
                continue
            try:
                refs[
                    (
                        row["case_id"],
                        float(row["checkpoint_time"]),
                        float(row["local_dt"]),
                    )
                ] = float(row["trajectory_bch_proxy_2err_over_dt2"])
            except (KeyError, TypeError, ValueError):
                continue
    return refs


def safe_ratio(numerator: float, denominator: float) -> float | str:
    if not np.isfinite(denominator) or denominator <= 0.0:
        return ""
    return numerator / denominator


def build_deterministic_orders(
    raw_pauli_keys: Sequence[ablation.PauliKey],
    final_coefficients: dict[ablation.PauliKey, complex],
    n_qubits: int,
    fermion_hamiltonian: Any,
    tolerance: float,
) -> tuple[dict[str, tuple[int, ...]], int]:
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


def build_exact_sector_problem(
    fermion_hamiltonian: Any,
    n_qubits: int,
    n_electrons: int,
    tolerance: float,
    spin_preserving: bool,
) -> tuple[Any, np.ndarray, np.ndarray]:
    reference_determinant = np.array(
        [index < n_electrons for index in range(n_qubits)],
        dtype=bool,
    )
    fermion_without_identity = remove_fermionic_identity(
        fermion_hamiltonian,
        tolerance,
    )
    sparse_hamiltonian = get_number_preserving_sparse_operator(
        fermion_without_identity,
        num_qubits=n_qubits,
        num_electrons=n_electrons,
        spin_preserving=spin_preserving,
        reference_determinant=reference_determinant,
        excitation_level=None,
    )
    if spin_preserving:
        basis_indices = spin_sector_basis_indices(n_qubits, n_electrons)
    else:
        basis_indices = number_sector_basis_indices(n_qubits, n_electrons)
    if sparse_hamiltonian.shape[0] != basis_indices.size:
        raise ValueError(
            "Exact sector dimension mismatch: "
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


def local_metrics(
    exact_target_sector: np.ndarray,
    exact_basis_indices: np.ndarray,
    approximate_full: np.ndarray,
    number_basis_indices: np.ndarray,
) -> dict[str, float]:
    approximate_sector = approximate_full[exact_basis_indices]
    overlap = np.vdot(exact_target_sector, approximate_sector)
    overlap_abs = min(1.0, max(0.0, float(abs(overlap))))

    # Direct phase alignment is more stable than sqrt(2 - 2|overlap|) when
    # the local error is extremely small.
    if abs(overlap) > 0.0:
        phase = np.exp(-1j * np.angle(overlap))
    else:
        phase = 1.0 + 0.0j
    aligned = phase * approximate_full
    exact_full = np.zeros_like(approximate_full)
    exact_full[exact_basis_indices] = exact_target_sector
    phase_error = float(np.linalg.norm(exact_full - aligned))

    spin_weight = float(np.sum(np.abs(approximate_sector) ** 2).real)
    number_weight = float(
        np.sum(np.abs(approximate_full[number_basis_indices]) ** 2).real
    )
    spin_weight = min(1.0, max(0.0, spin_weight))
    number_weight = min(1.0, max(0.0, number_weight))
    return {
        "local_overlap_abs": overlap_abs,
        "local_one_minus_overlap": max(0.0, 1.0 - overlap_abs),
        "local_infidelity": max(0.0, 1.0 - overlap_abs**2),
        "local_phase_aligned_2norm_error": phase_error,
        "spin_sector_leakage_after_one_step": 1.0 - spin_weight,
        "particle_number_leakage_after_one_step": 1.0 - number_weight,
    }


def append_row(
    writer: csv.DictWriter,
    stream: Any,
    row: dict[str, Any],
) -> None:
    writer.writerow({field: row.get(field, "") for field in FIELDNAMES})
    stream.flush()


def main() -> int:
    args = parse_args()
    if args.dt <= 0.0:
        raise ValueError("--dt must be positive.")
    checkpoints = sorted(set(float(value) for value in args.checkpoints))
    if not checkpoints or checkpoints[0] < 0.0:
        raise ValueError("Checkpoint times must be nonnegative.")

    interaction, n_qubits = ablation.load_interaction_operator(args.tensor)
    metadata = ablation.parse_case_metadata(args.tensor, n_qubits)
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
        raise ValueError("The identity-free JW Hamiltonian has no Pauli terms.")

    print(f"Case: {metadata.case_id}")
    print(f"Tensor: {args.tensor}")
    print(f"Qubits: {n_qubits}; checkpoints={checkpoints}; local dt={args.dt}")

    print("Building deterministic schedules...")
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
        name: [compiled_raw_terms[index] for index in orders[name]]
        for name in args.schedules
    }

    print("Building exact symmetry-sector Hamiltonian once...")
    exact_build_start = time.perf_counter()
    sparse_hamiltonian, exact_basis_indices, exact_sector_state = (
        build_exact_sector_problem(
            fermion_hamiltonian,
            n_qubits,
            metadata.active_occupied,
            args.tolerance,
            spin_preserving=not args.no_spin_sector,
        )
    )
    print(
        "Exact sector dimension: "
        f"{exact_sector_state.size}; build={time.perf_counter()-exact_build_start:.3f}s"
    )
    number_basis_indices = number_sector_basis_indices(
        n_qubits,
        metadata.active_occupied,
    )
    full_dimension = 1 << n_qubits

    warm_up_numba()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    completed = completed_keys(args.output)
    reference_proxies = read_reference_proxies(args.output)
    write_header = not args.output.exists() or args.output.stat().st_size == 0

    with args.output.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        if write_header:
            writer.writeheader()
            stream.flush()

        current_time = 0.0
        current_sector = exact_sector_state.copy()

        for checkpoint in checkpoints:
            if checkpoint < current_time - 1.0e-14:
                raise RuntimeError("Checkpoints must be processed in increasing order.")

            advance_seconds = 0.0
            delta = checkpoint - current_time
            if delta > 1.0e-14:
                start = time.perf_counter()
                current_sector = np.asarray(
                    expm_multiply((-1j * delta) * sparse_hamiltonian, current_sector),
                    dtype=np.complex128,
                )
                advance_seconds = time.perf_counter() - start
                current_time = checkpoint

            start = time.perf_counter()
            exact_target_sector = np.asarray(
                expm_multiply((-1j * args.dt) * sparse_hamiltonian, current_sector),
                dtype=np.complex128,
            )
            target_seconds = time.perf_counter() - start
            local_initial_full = scatter_sector_state(
                current_sector,
                exact_basis_indices,
                full_dimension,
            )

            checkpoint_reference_proxy = reference_proxies.get(
                (metadata.case_id, checkpoint, args.dt)
            )

            # Always process reference first so same-checkpoint ratios are known.
            schedule_order = [REFERENCE] + [
                name for name in args.schedules if name != REFERENCE
            ]
            schedule_order = [name for name in schedule_order if name in args.schedules]

            for schedule in schedule_order:
                key = (metadata.case_id, schedule, checkpoint, args.dt)
                if key in completed:
                    print(f"SKIP t={checkpoint:g} {schedule}")
                    continue

                print(f"RUN  t={checkpoint:g} {schedule}")
                row = blank_row()
                row.update(
                    {
                        "status": "success",
                        "case_id": metadata.case_id,
                        "tensor_path": str(args.tensor),
                        "molecule": metadata.molecule,
                        "bond_length": metadata.bond_length,
                        "basis": metadata.basis,
                        "active_occupied": metadata.active_occupied,
                        "active_vacant": metadata.active_vacant,
                        "n_qubits": n_qubits,
                        "number_of_fermionic_terms": number_of_fermionic_terms,
                        "number_of_pauli_terms": len(raw_pauli_keys),
                        "schedule": schedule,
                        "checkpoint_time": checkpoint,
                        "local_dt": args.dt,
                        "pauli_order_hash": ablation.hash_integer_order(orders[schedule]),
                        "exact_sector_dimension": exact_sector_state.size,
                        "exact_checkpoint_advance_seconds": advance_seconds,
                        "exact_local_target_seconds": target_seconds,
                        "coefficient_tolerance": args.tolerance,
                    }
                )

                try:
                    start = time.perf_counter()
                    approximate_full, _ = evolve_trotter_state(
                        initial_state=local_initial_full,
                        terms=ordered_terms[schedule],
                        formula_order=1,
                        trotter_steps=1,
                        evolution_time=args.dt,
                        parallel_threshold=args.parallel_threshold,
                    )
                    row["local_trotter_seconds"] = time.perf_counter() - start
                    metrics = local_metrics(
                        exact_target_sector,
                        exact_basis_indices,
                        approximate_full,
                        number_basis_indices,
                    )
                    row.update(metrics)
                    phase_error = float(metrics["local_phase_aligned_2norm_error"])
                    scaled = phase_error / (args.dt**2)
                    proxy = 2.0 * scaled
                    row["local_phase_aligned_2norm_error_over_dt2"] = scaled
                    row["trajectory_bch_proxy_2err_over_dt2"] = proxy

                    if schedule == REFERENCE:
                        checkpoint_reference_proxy = proxy
                        reference_proxies[(metadata.case_id, checkpoint, args.dt)] = proxy
                        row["trajectory_bch_proxy_ratio_to_signed"] = 1.0
                    elif checkpoint_reference_proxy is not None:
                        row["trajectory_bch_proxy_ratio_to_signed"] = safe_ratio(
                            proxy,
                            checkpoint_reference_proxy,
                        )
                except Exception as exc:  # Preserve append-only diagnostics.
                    row["status"] = "failed"
                    row["error_message"] = f"{type(exc).__name__}: {exc}"
                    append_row(writer, stream, row)
                    raise

                append_row(writer, stream, row)
                completed.add(key)

    print(f"Wrote: {args.output}")
    print()
    print("Primary plot to make from the CSV:")
    print("  x = checkpoint_time")
    print("  y = trajectory_bch_proxy_2err_over_dt2 (log scale)")
    print("If round-robin crosses below signed fermionic as t increases, the")
    print("Be2 reversal is explained by state-dependent local error changing")
    print("along the trajectory. If it does not cross, the reversal instead")
    print("points to time-integrated vector-direction/cancellation effects.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())