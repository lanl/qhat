#!/usr/bin/env python3
"""Matrix-free B2 active-space Trotter-ordering benchmark.

This script is designed for the QHAT ``L-sweep`` branch. It reuses the branch's
existing ordering construction from ``analysis/benchmark_L_sweep_trotter.py``
so that the three compared methods are unchanged:

    jw_raw
    jw_coloring
    fermionic_coloring

Unlike the existing benchmark, this file never constructs a dense 2**n by 2**n
Hamiltonian or unitary. Trotter evolution is applied directly to a state vector
with Numba-compiled Pauli actions. The exact reference is computed in the fixed
particle-number and fixed-Sz determinant sector with OpenFermion's
``get_number_preserving_sparse_operator`` and SciPy ``expm_multiply``.

Place this file in ``analysis/`` and run it from the QHAT repository root.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from numba import njit, prange
from scipy.sparse.linalg import expm_multiply

from openfermion import FermionOperator, get_fermion_operator, jordan_wigner

try:
    from openfermion import get_number_preserving_sparse_operator
except ImportError:  # Compatibility with older OpenFermion exports.
    from openfermion.linalg import get_number_preserving_sparse_operator

try:
    from qhat.analysis.benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        SUPPORTED_ORDERINGS,
        build_hermitian_fermion_terms,
        build_orderings,
        clean_fermion_operator,
        discover_tensor_paths,
        load_interaction_operator,
        parse_case_metadata,
        real_coefficient,
    )
except ImportError:
    # This path works when the file is placed directly in analysis/ and run as
    # ``python analysis/benchmark_b2_active_spaces_matrix_free.py``.
    from benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        SUPPORTED_ORDERINGS,
        build_hermitian_fermion_terms,
        build_orderings,
        clean_fermion_operator,
        discover_tensor_paths,
        load_interaction_operator,
        parse_case_metadata,
        real_coefficient,
    )


PauliKey = tuple[tuple[int, str], ...]

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
    "ordering",
    "formula_order",
    "trotter_steps",
    "trotter_dt",
    "evolution_time",
    "nominal_exponential_count",
    "ordering_build_time_seconds",
    "graph_vertices",
    "graph_edges",
    "number_of_colors",
    "graph_time_seconds",
    "coloring_time_seconds",
    "exact_sector_dimension",
    "exact_build_time_seconds",
    "exact_evolution_time_seconds",
    "trotter_runtime_seconds",
    "state_overlap_abs",
    "state_infidelity",
    "state_vector_2norm_error",
    "phase_aligned_state_2norm_error",
    "particle_number_leakage",
    "spin_sector_leakage",
    "state_infidelity_ratio_to_jw_raw",
    "coefficient_tolerance",
]

SUMMARY_FIELDNAMES = [
    "case_id",
    "bond_length",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "formula_order",
    "trotter_steps",
    "jw_raw_state_infidelity",
    "jw_coloring_state_infidelity",
    "fermionic_coloring_state_infidelity",
    "jw_coloring_ratio_to_jw_raw",
    "fermionic_coloring_ratio_to_jw_raw",
    "best_ordering",
]


# -----------------------------------------------------------------------------
# Numba-compiled Pauli-string state actions
# -----------------------------------------------------------------------------


@njit(cache=True)
def _parity_u64(value: np.uint64) -> int:
    parity = 0
    while value:
        parity ^= 1
        value &= value - np.uint64(1)
    return parity


@njit(cache=True)
def _apply_pauli_exponential_serial(
    state: np.ndarray,
    output: np.ndarray,
    flip_mask: np.uint64,
    sign_mask: np.uint64,
    y_mod_four: int,
    angle: float,
) -> None:
    cosine = math.cos(angle)
    sine = math.sin(angle)
    minus_i_sine = -1j * sine

    for output_index in range(state.size):
        source_index_u64 = np.uint64(output_index) ^ flip_mask
        source_index = int(source_index_u64)
        pauli_amplitude = state[source_index]

        if _parity_u64(source_index_u64 & sign_mask):
            pauli_amplitude = -pauli_amplitude

        if y_mod_four == 1:
            pauli_amplitude *= 1j
        elif y_mod_four == 2:
            pauli_amplitude = -pauli_amplitude
        elif y_mod_four == 3:
            pauli_amplitude *= -1j

        output[output_index] = (
            cosine * state[output_index]
            + minus_i_sine * pauli_amplitude
        )


@njit(cache=True, parallel=True)
def _apply_pauli_exponential_parallel(
    state: np.ndarray,
    output: np.ndarray,
    flip_mask: np.uint64,
    sign_mask: np.uint64,
    y_mod_four: int,
    angle: float,
) -> None:
    cosine = math.cos(angle)
    sine = math.sin(angle)
    minus_i_sine = -1j * sine

    for output_index in prange(state.size):
        source_index_u64 = np.uint64(output_index) ^ flip_mask
        source_index = int(source_index_u64)
        pauli_amplitude = state[source_index]

        if _parity_u64(source_index_u64 & sign_mask):
            pauli_amplitude = -pauli_amplitude

        if y_mod_four == 1:
            pauli_amplitude *= 1j
        elif y_mod_four == 2:
            pauli_amplitude = -pauli_amplitude
        elif y_mod_four == 3:
            pauli_amplitude *= -1j

        output[output_index] = (
            cosine * state[output_index]
            + minus_i_sine * pauli_amplitude
        )


def warm_up_numba() -> None:
    """Compile both kernels before benchmark timings begin."""
    state = np.zeros(4, dtype=np.complex128)
    output = np.empty_like(state)
    state[0] = 1.0
    _apply_pauli_exponential_serial(
        state,
        output,
        np.uint64(1),
        np.uint64(0),
        0,
        0.1,
    )
    _apply_pauli_exponential_parallel(
        state,
        output,
        np.uint64(1),
        np.uint64(0),
        0,
        0.1,
    )


def pauli_descriptor(
    pauli_key: PauliKey,
    n_qubits: int,
) -> tuple[np.uint64, np.uint64, int]:
    """Convert an OpenFermion Pauli key to bit masks for state action."""
    flip_mask = 0
    sign_mask = 0
    number_of_y = 0

    for qubit, pauli in pauli_key:
        if not 0 <= qubit < n_qubits:
            raise ValueError(f"Qubit index {qubit} is outside {n_qubits} qubits.")
        bit = 1 << (n_qubits - 1 - qubit)

        if pauli in {"X", "Y"}:
            flip_mask |= bit
        if pauli in {"Y", "Z"}:
            sign_mask |= bit
        if pauli == "Y":
            number_of_y += 1
        if pauli not in {"X", "Y", "Z"}:
            raise ValueError(f"Unsupported Pauli label {pauli!r}.")

    return (
        np.uint64(flip_mask),
        np.uint64(sign_mask),
        number_of_y % 4,
    )


def compile_ordered_terms(
    pauli_keys: Sequence[PauliKey],
    final_coefficients: dict[PauliKey, complex],
    n_qubits: int,
    tolerance: float,
) -> list[tuple[float, np.uint64, np.uint64, int]]:
    compiled = []
    for pauli_key in pauli_keys:
        coefficient = real_coefficient(
            final_coefficients[pauli_key],
            tolerance,
        )
        flip_mask, sign_mask, y_mod_four = pauli_descriptor(
            pauli_key,
            n_qubits,
        )
        compiled.append(
            (coefficient, flip_mask, sign_mask, y_mod_four)
        )
    return compiled


def apply_one_exponential(
    state: np.ndarray,
    work: np.ndarray,
    term: tuple[float, np.uint64, np.uint64, int],
    duration: float,
    parallel_threshold: int,
) -> tuple[np.ndarray, np.ndarray]:
    coefficient, flip_mask, sign_mask, y_mod_four = term
    angle = coefficient * duration

    if state.size >= parallel_threshold:
        _apply_pauli_exponential_parallel(
            state,
            work,
            flip_mask,
            sign_mask,
            y_mod_four,
            angle,
        )
    else:
        _apply_pauli_exponential_serial(
            state,
            work,
            flip_mask,
            sign_mask,
            y_mod_four,
            angle,
        )
    return work, state


def apply_second_order_slice(
    state: np.ndarray,
    work: np.ndarray,
    terms: Sequence[tuple[float, np.uint64, np.uint64, int]],
    duration: float,
    parallel_threshold: int,
) -> tuple[np.ndarray, np.ndarray]:
    if not terms:
        return state, work

    half_duration = 0.5 * duration
    for term in terms[:-1]:
        state, work = apply_one_exponential(
            state,
            work,
            term,
            half_duration,
            parallel_threshold,
        )

    state, work = apply_one_exponential(
        state,
        work,
        terms[-1],
        duration,
        parallel_threshold,
    )

    for term in reversed(terms[:-1]):
        state, work = apply_one_exponential(
            state,
            work,
            term,
            half_duration,
            parallel_threshold,
        )

    return state, work


def evolve_trotter_state(
    initial_state: np.ndarray,
    terms: Sequence[tuple[float, np.uint64, np.uint64, int]],
    formula_order: int,
    trotter_steps: int,
    evolution_time: float,
    parallel_threshold: int,
) -> tuple[np.ndarray, int]:
    """Apply first-, second-, or fourth-order product formulas."""
    state = initial_state.copy()
    work = np.empty_like(state)
    dt = evolution_time / trotter_steps
    number_of_terms = len(terms)

    for _ in range(trotter_steps):
        if formula_order == 1:
            for term in terms:
                state, work = apply_one_exponential(
                    state,
                    work,
                    term,
                    dt,
                    parallel_threshold,
                )
        elif formula_order == 2:
            state, work = apply_second_order_slice(
                state,
                work,
                terms,
                dt,
                parallel_threshold,
            )
        elif formula_order == 4:
            p = 1.0 / (4.0 - 4.0 ** (1.0 / 3.0))
            for scale in (p, p, 1.0 - 4.0 * p, p, p):
                state, work = apply_second_order_slice(
                    state,
                    work,
                    terms,
                    scale * dt,
                    parallel_threshold,
                )
        else:
            raise ValueError(f"Unsupported formula order {formula_order}.")

    if formula_order == 1:
        nominal_exponential_count = trotter_steps * number_of_terms
    elif formula_order == 2:
        nominal_exponential_count = trotter_steps * (
            2 * number_of_terms - 1
        )
    else:
        nominal_exponential_count = trotter_steps * (
            5 * (2 * number_of_terms - 1)
        )

    return state, nominal_exponential_count


# -----------------------------------------------------------------------------
# Exact fixed-particle and fixed-Sz reference
# -----------------------------------------------------------------------------


def remove_fermionic_identity(
    fermion_hamiltonian: FermionOperator,
    tolerance: float,
) -> FermionOperator:
    result = FermionOperator()
    for term_key, coefficient in fermion_hamiltonian.terms.items():
        if term_key and abs(coefficient) > tolerance:
            result += FermionOperator(term_key, coefficient)
    return clean_fermion_operator(result, tolerance)


def determinant_to_index(determinant: np.ndarray) -> int:
    weights = 1 << np.arange(determinant.size - 1, -1, -1)
    return int(np.dot(determinant.astype(np.int64), weights))


def number_sector_basis_indices(
    n_qubits: int,
    n_electrons: int,
) -> np.ndarray:
    """Match OpenFermion's excitation-ranked determinant ordering."""
    reference = np.array(
        [index < n_electrons for index in range(n_qubits)],
        dtype=bool,
    )
    occupied = np.where(reference)[0]
    unoccupied = np.where(~reference)[0]
    indices: list[int] = []

    maximum_order = min(len(occupied), len(unoccupied))
    for order in range(maximum_order + 1):
        for occupied_removed, unoccupied_added in itertools.product(
            itertools.combinations(occupied, order),
            itertools.combinations(unoccupied, order),
        ):
            determinant = reference.copy()
            determinant[list(occupied_removed)] = False
            determinant[list(unoccupied_added)] = True
            indices.append(determinant_to_index(determinant))

    return np.asarray(indices, dtype=np.int64)


def spin_sector_basis_indices(
    n_qubits: int,
    n_electrons: int,
) -> np.ndarray:
    """Match OpenFermion's spin-preserving excitation-ranked ordering."""
    if n_qubits % 2 != 0:
        raise ValueError("Spin-preserving reference requires an even qubit count.")

    reference = np.array(
        [index < n_electrons for index in range(n_qubits)],
        dtype=bool,
    )

    occupied_alpha = np.where(reference[::2])[0] * 2
    unoccupied_alpha = np.where(~reference[::2])[0] * 2
    occupied_beta = np.where(reference[1::2])[0] * 2 + 1
    unoccupied_beta = np.where(~reference[1::2])[0] * 2 + 1

    alpha_maximum = min(len(occupied_alpha), len(unoccupied_alpha))
    beta_maximum = min(len(occupied_beta), len(unoccupied_beta))
    indices: list[int] = []

    for total_order in range(alpha_maximum + beta_maximum + 1):
        for alpha_order in range(alpha_maximum + 1):
            beta_order = total_order - alpha_order
            if beta_order < 0 or beta_order > beta_maximum:
                continue

            products = itertools.product(
                itertools.combinations(occupied_alpha, alpha_order),
                itertools.combinations(unoccupied_alpha, alpha_order),
                itertools.combinations(occupied_beta, beta_order),
                itertools.combinations(unoccupied_beta, beta_order),
            )
            for (
                alpha_removed,
                alpha_added,
                beta_removed,
                beta_added,
            ) in products:
                determinant = reference.copy()
                determinant[list(alpha_removed)] = False
                determinant[list(alpha_added)] = True
                determinant[list(beta_removed)] = False
                determinant[list(beta_added)] = True
                indices.append(determinant_to_index(determinant))

    return np.asarray(indices, dtype=np.int64)


def build_hartree_fock_state(
    n_qubits: int,
    n_electrons: int,
) -> np.ndarray:
    basis_index = sum(
        1 << (n_qubits - 1 - qubit)
        for qubit in range(n_electrons)
    )
    state = np.zeros(2**n_qubits, dtype=np.complex128)
    state[basis_index] = 1.0
    return state


def exact_reference_state(
    fermion_hamiltonian: FermionOperator,
    n_qubits: int,
    n_electrons: int,
    evolution_time: float,
    tolerance: float,
    spin_preserving: bool,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Build and evolve the exact sparse Hamiltonian in a symmetry sector."""
    reference_determinant = np.array(
        [index < n_electrons for index in range(n_qubits)],
        dtype=bool,
    )
    fermion_without_identity = remove_fermionic_identity(
        fermion_hamiltonian,
        tolerance,
    )

    build_start = time.perf_counter()
    sparse_hamiltonian = get_number_preserving_sparse_operator(
        fermion_without_identity,
        num_qubits=n_qubits,
        num_electrons=n_electrons,
        spin_preserving=spin_preserving,
        reference_determinant=reference_determinant,
        excitation_level=None,
    )
    exact_build_time = time.perf_counter() - build_start

    if spin_preserving:
        basis_indices = spin_sector_basis_indices(n_qubits, n_electrons)
    else:
        basis_indices = number_sector_basis_indices(n_qubits, n_electrons)

    if sparse_hamiltonian.shape[0] != basis_indices.size:
        raise ValueError(
            "OpenFermion sector dimension does not match the reconstructed "
            f"basis order: {sparse_hamiltonian.shape[0]} != "
            f"{basis_indices.size}."
        )

    initial_sector_state = np.zeros(
        sparse_hamiltonian.shape[0],
        dtype=np.complex128,
    )
    initial_sector_state[0] = 1.0

    evolution_start = time.perf_counter()
    exact_sector_state = expm_multiply(
        (-1j * evolution_time) * sparse_hamiltonian,
        initial_sector_state,
    )
    exact_evolution_time = time.perf_counter() - evolution_start

    return (
        np.asarray(exact_sector_state, dtype=np.complex128),
        basis_indices,
        exact_build_time,
        exact_evolution_time,
    )


def compare_states(
    exact_sector_state: np.ndarray,
    exact_basis_indices: np.ndarray,
    approximate_full_state: np.ndarray,
    number_basis_indices: np.ndarray,
) -> dict[str, float]:
    approximate_in_exact_sector = approximate_full_state[exact_basis_indices]
    overlap = np.vdot(exact_sector_state, approximate_in_exact_sector)
    overlap_abs = min(1.0, max(0.0, float(abs(overlap))))
    fidelity = min(1.0, max(0.0, overlap_abs**2))

    number_weight = float(
        np.sum(np.abs(approximate_full_state[number_basis_indices]) ** 2).real
    )
    spin_weight = float(
        np.sum(np.abs(approximate_in_exact_sector) ** 2).real
    )
    number_weight = min(1.0, max(0.0, number_weight))
    spin_weight = min(1.0, max(0.0, spin_weight))

    return {
        "state_overlap_abs": overlap_abs,
        "state_infidelity": 1.0 - fidelity,
        "state_vector_2norm_error": math.sqrt(
            max(0.0, 2.0 - 2.0 * float(overlap.real))
        ),
        "phase_aligned_state_2norm_error": math.sqrt(
            max(0.0, 2.0 - 2.0 * overlap_abs)
        ),
        "particle_number_leakage": 1.0 - number_weight,
        "spin_sector_leakage": 1.0 - spin_weight,
    }


# -----------------------------------------------------------------------------
# CSV resume, summary, and plots
# -----------------------------------------------------------------------------


def blank_row() -> dict[str, Any]:
    return {field: "" for field in FIELDNAMES}


def load_resume_data(
    output_path: Path,
) -> tuple[
    set[tuple[str, str, int, int, float]],
    dict[tuple[str, int, int, float], float],
]:
    completed: set[tuple[str, str, int, int, float]] = set()
    raw_infidelities: dict[tuple[str, int, int, float], float] = {}

    if not output_path.exists():
        return completed, raw_infidelities

    with output_path.open("r", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row.get("status") != "success":
                continue
            try:
                key = (
                    row["case_id"],
                    row["ordering"],
                    int(row["formula_order"]),
                    int(row["trotter_steps"]),
                    float(row["evolution_time"]),
                )
                infidelity = float(row["state_infidelity"])
            except (KeyError, TypeError, ValueError):
                continue

            completed.add(key)
            if row["ordering"] == "jw_raw":
                raw_infidelities[(key[0], key[2], key[3], key[4])] = (
                    infidelity
                )

    return completed, raw_infidelities


def safe_ratio(numerator: float, denominator: float, tolerance: float) -> str | float:
    if denominator <= tolerance:
        return ""
    return numerator / denominator


def read_successful_rows(output_path: Path) -> list[dict[str, str]]:
    if not output_path.exists():
        return []
    with output_path.open("r", newline="", encoding="utf-8") as stream:
        return [
            row
            for row in csv.DictReader(stream)
            if row.get("status") == "success"
        ]


def write_comparison_summary(
    output_path: Path,
    summary_path: Path,
) -> None:
    rows = read_successful_rows(output_path)
    grouped: dict[
        tuple[str, str, int, int, int, int, int],
        dict[str, float],
    ] = defaultdict(dict)

    for row in rows:
        key = (
            row["case_id"],
            row["bond_length"],
            int(row["active_occupied"]),
            int(row["active_vacant"]),
            int(row["n_qubits"]),
            int(row["formula_order"]),
            int(row["trotter_steps"]),
        )
        grouped[key][row["ordering"]] = float(row["state_infidelity"])

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()

        for key in sorted(grouped, key=lambda item: (item[4], item[5], item[6])):
            (
                case_id,
                bond_length,
                active_occupied,
                active_vacant,
                n_qubits,
                formula_order,
                trotter_steps,
            ) = key
            errors = grouped[key]
            raw_error = errors.get("jw_raw", math.nan)
            jw_coloring_error = errors.get("jw_coloring", math.nan)
            fermionic_error = errors.get("fermionic_coloring", math.nan)
            finite_errors = {
                name: value
                for name, value in errors.items()
                if math.isfinite(value)
            }
            best_ordering = (
                min(finite_errors, key=finite_errors.get)
                if finite_errors
                else ""
            )

            writer.writerow(
                {
                    "case_id": case_id,
                    "bond_length": bond_length,
                    "active_occupied": active_occupied,
                    "active_vacant": active_vacant,
                    "n_qubits": n_qubits,
                    "formula_order": formula_order,
                    "trotter_steps": trotter_steps,
                    "jw_raw_state_infidelity": (
                        raw_error if math.isfinite(raw_error) else ""
                    ),
                    "jw_coloring_state_infidelity": (
                        jw_coloring_error
                        if math.isfinite(jw_coloring_error)
                        else ""
                    ),
                    "fermionic_coloring_state_infidelity": (
                        fermionic_error
                        if math.isfinite(fermionic_error)
                        else ""
                    ),
                    "jw_coloring_ratio_to_jw_raw": (
                        jw_coloring_error / raw_error
                        if math.isfinite(jw_coloring_error)
                        and math.isfinite(raw_error)
                        and raw_error > 0.0
                        else ""
                    ),
                    "fermionic_coloring_ratio_to_jw_raw": (
                        fermionic_error / raw_error
                        if math.isfinite(fermionic_error)
                        and math.isfinite(raw_error)
                        and raw_error > 0.0
                        else ""
                    ),
                    "best_ordering": best_ordering,
                }
            )


def write_comparison_plots(
    output_path: Path,
    plot_directory: Path,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is unavailable; skipping comparison plots.")
        return

    rows = read_successful_rows(output_path)
    if not rows:
        return

    plot_directory.mkdir(parents=True, exist_ok=True)
    combinations = sorted(
        {
            (int(row["formula_order"]), int(row["trotter_steps"]))
            for row in rows
        }
    )

    for formula_order, trotter_steps in combinations:
        figure, axis = plt.subplots(figsize=(8.5, 5.5))
        plotted = False

        for ordering in SUPPORTED_ORDERINGS:
            selected = [
                row
                for row in rows
                if row["ordering"] == ordering
                and int(row["formula_order"]) == formula_order
                and int(row["trotter_steps"]) == trotter_steps
            ]
            selected.sort(key=lambda row: int(row["n_qubits"]))
            if not selected:
                continue

            x_values = [int(row["n_qubits"]) for row in selected]
            y_values = [
                max(float(row["state_infidelity"]), np.finfo(float).tiny)
                for row in selected
            ]
            axis.plot(x_values, y_values, marker="o", label=ordering)
            plotted = True

        if not plotted:
            plt.close(figure)
            continue

        axis.set_yscale("log")
        axis.set_xlabel("Active-space qubits")
        axis.set_ylabel("State infidelity")
        axis.set_title(
            "B2/STO-6G Trotter ordering comparison\n"
            f"formula order {formula_order}, steps {trotter_steps}"
        )
        axis.grid(True, which="both", alpha=0.3)
        axis.legend()
        figure.tight_layout()

        plot_path = plot_directory / (
            f"b2_state_infidelity_order-{formula_order}_"
            f"steps-{trotter_steps}.png"
        )
        figure.savefig(plot_path, dpi=200)
        plt.close(figure)


# -----------------------------------------------------------------------------
# Benchmark driver
# -----------------------------------------------------------------------------


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare raw JW, JW coloring, and fermionic-color-induced JW "
            "orderings for all B2/STO-6G active spaces without dense matrices."
        )
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=Path("hamiltonian_generator/b2_active_space_library"),
        help="Root containing the B2 *.tensors.npz files.",
    )
    parser.add_argument(
        "--steps",
        type=int,
        nargs="+",
        default=[1],
        help=(
            "Trotter step counts. The default is 1 so all nine spaces, "
            "including 20 qubits, receive at least one comparison."
        ),
    )
    parser.add_argument(
        "--formula-orders",
        type=int,
        nargs="+",
        choices=[1, 2, 4],
        default=[1, 2],
        help="Product-formula orders to benchmark.",
    )
    parser.add_argument(
        "--orderings",
        nargs="+",
        choices=SUPPORTED_ORDERINGS,
        default=list(SUPPORTED_ORDERINGS),
        help="Ordering methods to compare.",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=1.0,
        dest="evolution_time",
        help="Total evolution time.",
    )
    parser.add_argument(
        "--min-qubits",
        type=int,
        default=4,
        help="Smallest active-space size to include.",
    )
    parser.add_argument(
        "--max-qubits",
        type=int,
        default=20,
        help="Largest active-space size to include.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/b2_active_space_trotter_results.csv"),
        help="Detailed result CSV.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("analysis/b2_active_space_comparison_summary.csv"),
        help="Wide comparison summary CSV.",
    )
    parser.add_argument(
        "--plot-directory",
        type=Path,
        default=Path("analysis/b2_active_space_figures"),
        help="Directory for infidelity comparison figures.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help="Coefficient and commutator compression tolerance.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Append to an existing CSV and skip completed combinations.",
    )
    parser.add_argument(
        "--no-spin-sector",
        action="store_true",
        help=(
            "Use only particle-number restriction for the exact reference. "
            "The default also fixes Sz to the Hartree-Fock sector."
        ),
    )
    parser.add_argument(
        "--parallel-threshold",
        type=int,
        default=2**16,
        help=(
            "Use the parallel Numba Pauli kernel at or above this state-vector "
            "dimension."
        ),
    )
    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    if not args.library.exists():
        raise FileNotFoundError(args.library)
    if any(step <= 0 for step in args.steps):
        raise ValueError("All Trotter step counts must be positive.")
    if args.evolution_time <= 0.0:
        raise ValueError("Evolution time must be positive.")
    if args.min_qubits < 1 or args.max_qubits < args.min_qubits:
        raise ValueError("Invalid qubit range.")
    if args.max_qubits > 20:
        raise ValueError("This B2 study contains at most 20 qubits.")
    if args.tolerance <= 0.0:
        raise ValueError("--tolerance must be positive.")
    if args.parallel_threshold <= 0:
        raise ValueError("--parallel-threshold must be positive.")


def select_b2_tensor_paths(
    library: Path,
    min_qubits: int,
    max_qubits: int,
) -> list[Path]:
    paths = discover_tensor_paths(library, case_pattern=None, limit=None)
    selected: list[Path] = []

    for path in paths:
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
        if not min_qubits <= n_qubits <= max_qubits:
            continue
        selected.append(path)

    selected.sort(
        key=lambda path: (
            int(path.name.split("_as-")[1].split("-")[0])
            + int(path.name.split("_as-")[1].split("-")[1].split(".")[0]),
            str(path),
        )
    )
    return selected


def add_metadata(row: dict[str, Any], metadata: Any) -> None:
    row.update(
        {
            "case_id": metadata.case_id,
            "tensor_path": str(metadata.tensor_path),
            "molecule": metadata.molecule,
            "bond_length": metadata.bond_length,
            "basis": metadata.basis,
            "active_occupied": metadata.active_occupied,
            "active_vacant": metadata.active_vacant,
            "n_qubits": metadata.n_qubits,
        }
    )


def benchmark_case(
    tensor_path: Path,
    requested_orderings: Sequence[str],
    formula_orders: Sequence[int],
    steps_list: Sequence[int],
    evolution_time: float,
    tolerance: float,
    spin_preserving: bool,
    parallel_threshold: int,
    completed: set[tuple[str, str, int, int, float]],
    raw_infidelities: dict[tuple[str, int, int, float], float],
    writer: csv.DictWriter,
    output_stream: Any,
) -> None:
    interaction, n_qubits = load_interaction_operator(tensor_path)
    metadata = parse_case_metadata(tensor_path, n_qubits)

    fermion_hamiltonian = clean_fermion_operator(
        get_fermion_operator(interaction),
        tolerance,
    )
    full_jw_hamiltonian = jordan_wigner(fermion_hamiltonian)
    full_jw_hamiltonian.compress(abs_tol=tolerance)

    orderings, hermitian_terms, final_coefficients = build_orderings(
        fermion_hamiltonian=fermion_hamiltonian,
        full_jw_hamiltonian=full_jw_hamiltonian,
        requested_orderings=requested_orderings,
        tolerance=tolerance,
    )

    if not hermitian_terms:
        hermitian_terms = build_hermitian_fermion_terms(
            fermion_hamiltonian,
            tolerance,
        )

    compiled_orderings = {
        name: compile_ordered_terms(
            orderings[name].pauli_keys,
            final_coefficients,
            n_qubits,
            tolerance,
        )
        for name in requested_orderings
    }

    print("  Building exact symmetry-sector reference...")
    (
        exact_sector_state,
        exact_basis_indices,
        exact_build_time,
        exact_evolution_time,
    ) = exact_reference_state(
        fermion_hamiltonian=fermion_hamiltonian,
        n_qubits=n_qubits,
        n_electrons=metadata.active_occupied,
        evolution_time=evolution_time,
        tolerance=tolerance,
        spin_preserving=spin_preserving,
    )

    number_basis_indices = number_sector_basis_indices(
        n_qubits,
        metadata.active_occupied,
    )
    initial_state = build_hartree_fock_state(
        n_qubits,
        metadata.active_occupied,
    )

    processing_order = ["jw_raw"] + [
        name for name in requested_orderings if name != "jw_raw"
    ]

    for formula_order in formula_orders:
        for trotter_steps in steps_list:
            raw_key = (
                metadata.case_id,
                formula_order,
                trotter_steps,
                evolution_time,
            )

            for ordering_name in processing_order:
                result_key = (
                    metadata.case_id,
                    ordering_name,
                    formula_order,
                    trotter_steps,
                    evolution_time,
                )
                if result_key in completed:
                    print(
                        "  SKIP completed: "
                        f"{ordering_name}, order={formula_order}, "
                        f"steps={trotter_steps}"
                    )
                    continue

                ordering = orderings[ordering_name]
                print(
                    "  RUN  "
                    f"{ordering_name:20s} "
                    f"order={formula_order} steps={trotter_steps}"
                )
                trotter_start = time.perf_counter()
                approximate_state, nominal_exponential_count = (
                    evolve_trotter_state(
                        initial_state=initial_state,
                        terms=compiled_orderings[ordering_name],
                        formula_order=formula_order,
                        trotter_steps=trotter_steps,
                        evolution_time=evolution_time,
                        parallel_threshold=parallel_threshold,
                    )
                )
                trotter_runtime = time.perf_counter() - trotter_start

                metrics = compare_states(
                    exact_sector_state=exact_sector_state,
                    exact_basis_indices=exact_basis_indices,
                    approximate_full_state=approximate_state,
                    number_basis_indices=number_basis_indices,
                )

                infidelity = metrics["state_infidelity"]
                if ordering_name == "jw_raw":
                    raw_infidelities[raw_key] = infidelity
                    infidelity_ratio: str | float = 1.0
                else:
                    if raw_key not in raw_infidelities:
                        raise RuntimeError(
                            "jw_raw must be processed before comparison orderings."
                        )
                    infidelity_ratio = safe_ratio(
                        infidelity,
                        raw_infidelities[raw_key],
                        tolerance,
                    )

                row = blank_row()
                add_metadata(row, metadata)
                row.update(
                    {
                        "status": "success",
                        "number_of_fermionic_terms": len(hermitian_terms),
                        "number_of_pauli_terms": len(final_coefficients),
                        "ordering": ordering_name,
                        "formula_order": formula_order,
                        "trotter_steps": trotter_steps,
                        "trotter_dt": evolution_time / trotter_steps,
                        "evolution_time": evolution_time,
                        "nominal_exponential_count": nominal_exponential_count,
                        "ordering_build_time_seconds": ordering.build_time_seconds,
                        "graph_vertices": ordering.graph_vertices,
                        "graph_edges": ordering.graph_edges,
                        "number_of_colors": ordering.number_of_colors,
                        "graph_time_seconds": ordering.graph_time_seconds,
                        "coloring_time_seconds": ordering.coloring_time_seconds,
                        "exact_sector_dimension": exact_sector_state.size,
                        "exact_build_time_seconds": exact_build_time,
                        "exact_evolution_time_seconds": exact_evolution_time,
                        "trotter_runtime_seconds": trotter_runtime,
                        **metrics,
                        "state_infidelity_ratio_to_jw_raw": infidelity_ratio,
                        "coefficient_tolerance": tolerance,
                    }
                )
                writer.writerow(row)
                output_stream.flush()
                completed.add(result_key)

                print(
                    "       "
                    f"infidelity={infidelity:.6e}  "
                    f"ratio={infidelity_ratio}  "
                    f"number_leakage={metrics['particle_number_leakage']:.3e}"
                )


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    warm_up_numba()

    requested_orderings = list(dict.fromkeys(args.orderings))
    if "jw_raw" not in requested_orderings:
        requested_orderings.insert(0, "jw_raw")
        print(
            "Added jw_raw automatically because it is required for baseline "
            "ratios."
        )

    tensor_paths = select_b2_tensor_paths(
        args.library,
        args.min_qubits,
        args.max_qubits,
    )
    if not tensor_paths:
        raise FileNotFoundError(
            f"No B-B/STO-6G tensor files found under {args.library}."
        )

    print(f"B2 tensor cases selected: {len(tensor_paths)}")
    print(f"Orderings: {requested_orderings}")
    print(f"Formula orders: {args.formula_orders}")
    print(f"Trotter steps: {args.steps}")
    print(f"Evolution time: {args.evolution_time}")
    print(f"Spin-preserving exact sector: {not args.no_spin_sector}")
    print(f"Detailed output: {args.output}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.resume:
        completed, raw_infidelities = load_resume_data(args.output)
        file_mode = "a"
        write_header = (
            not args.output.exists() or args.output.stat().st_size == 0
        )
        print(f"Completed rows loaded: {len(completed)}")
    else:
        completed = set()
        raw_infidelities = {}
        file_mode = "w"
        write_header = True

    successful_cases = 0
    failed_cases = 0

    with args.output.open(
        file_mode,
        newline="",
        encoding="utf-8",
    ) as output_stream:
        writer = csv.DictWriter(output_stream, fieldnames=FIELDNAMES)
        if write_header:
            writer.writeheader()
            output_stream.flush()

        for case_index, tensor_path in enumerate(tensor_paths, start=1):
            print()
            print("=" * 88)
            print(f"[{case_index}/{len(tensor_paths)}] {tensor_path}")
            try:
                benchmark_case(
                    tensor_path=tensor_path,
                    requested_orderings=requested_orderings,
                    formula_orders=args.formula_orders,
                    steps_list=args.steps,
                    evolution_time=args.evolution_time,
                    tolerance=args.tolerance,
                    spin_preserving=not args.no_spin_sector,
                    parallel_threshold=args.parallel_threshold,
                    completed=completed,
                    raw_infidelities=raw_infidelities,
                    writer=writer,
                    output_stream=output_stream,
                )
            except Exception as error:
                failed_cases += 1
                print(
                    f"FAILED {tensor_path}: "
                    f"{type(error).__name__}: {error}"
                )
                failure_row = blank_row()
                failure_row.update(
                    {
                        "status": "failed",
                        "error_message": f"{type(error).__name__}: {error}",
                        "case_id": tensor_path.name.removesuffix(
                            ".tensors.npz"
                        ),
                        "tensor_path": str(tensor_path),
                        "coefficient_tolerance": args.tolerance,
                    }
                )
                writer.writerow(failure_row)
                output_stream.flush()
            else:
                successful_cases += 1

    write_comparison_summary(args.output, args.summary)
    write_comparison_plots(args.output, args.plot_directory)

    print()
    print("=" * 88)
    print(f"Cases selected:   {len(tensor_paths)}")
    print(f"Cases successful: {successful_cases}")
    print(f"Cases failed:     {failed_cases}")
    print(f"Detailed results: {args.output}")
    print(f"Comparison table: {args.summary}")
    print(f"Figures:          {args.plot_directory}")
    print("=" * 88)


if __name__ == "__main__":
    main()
