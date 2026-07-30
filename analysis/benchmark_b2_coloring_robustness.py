#!/usr/bin/env python3
"""Test the robustness of B2 coloring-induced Trotter orderings.

This script is designed for the QHAT ``L-sweep`` branch. Place it in
``analysis/`` and run it from the repository root.

For the selected B2/STO-6G active spaces, the script keeps each greedy-coloring
partition fixed and changes only the schedule of its color groups and nodes.
It compares:

    jw_raw

    jw_current
    jw_reverse_groups
    jw_random_groups
    jw_random_within
    jw_random_both

    fermionic_current
    fermionic_reverse_groups
    fermionic_random_groups
    fermionic_random_within
    fermionic_random_both

The JW variants schedule Pauli-string color groups directly. The fermionic
variants schedule Hermitian fermionic-term color groups and then induce a
permutation of the final combined Jordan-Wigner Pauli strings using the same
first-occurrence rule as ``benchmark_L_sweep_trotter.py``.

Every schedule is evaluated with:

* the state-dependent leading BCH commutator norm on the Hartree-Fock state;
* first-order, matrix-free Trotter state infidelity against the exact
  particle-number- and Sz-preserving reference state.

The graph coloring is performed once per graph and case. Random samples change
only group order, within-group order, or both. This isolates scheduling from
partition construction.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import networkx as nx
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
        deterministic_coloring_order,
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
        deterministic_coloring_order,
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


PauliKey = tuple[tuple[int, str], ...]

RANDOM_SCHEDULES = (
    "random_groups",
    "random_within",
    "random_both",
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
    "number_of_fermionic_terms",
    "number_of_pauli_terms",
    "graph_level",
    "schedule",
    "sample_index",
    "random_seed",
    "graph_vertices",
    "graph_edges",
    "graph_density",
    "number_of_colors",
    "group_sizes",
    "group_order",
    "within_group_order_hash",
    "pauli_order_hash",
    "trotter_steps",
    "trotter_dt",
    "evolution_time",
    "nominal_exponential_count",
    "bch2_hf_state_norm",
    "leading_state_error_norm_estimate",
    "bch_squared_ratio_to_raw",
    "schedule_build_time_seconds",
    "bch_evaluation_time_seconds",
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
    "state_infidelity_ratio_to_raw_jw",
    "state_infidelity_ratio_to_current_coloring",
    "coefficient_tolerance",
]

SUMMARY_FIELDNAMES = [
    "case_id",
    "n_qubits",
    "graph_level",
    "schedule",
    "number_of_samples",
    "minimum_bch2_hf_state_norm",
    "median_bch2_hf_state_norm",
    "maximum_bch2_hf_state_norm",
    "mean_bch2_hf_state_norm",
    "std_bch2_hf_state_norm",
    "minimum_state_infidelity",
    "median_state_infidelity",
    "maximum_state_infidelity",
    "mean_state_infidelity",
    "std_state_infidelity",
    "fraction_beating_raw_jw",
    "fraction_beating_current_coloring",
    "pearson_correlation_log_bch2_squared_vs_log_infidelity",
]


@dataclass(frozen=True)
class ScheduleSpec:
    """One fixed or random schedule before Trotter evolution."""

    graph_level: str
    schedule: str
    sample_index: int
    random_seed: int
    node_order: tuple[int, ...]
    group_order: tuple[int, ...]
    grouped_nodes: tuple[tuple[int, ...], ...]


@dataclass
class HFCommutatorEvaluator:
    """Fast state-dependent BCH evaluator for permutations of one Pauli set.

    For every noncommuting unordered Pauli pair (i, j), this class precomputes
    the amplitude of [H_i, H_j] acting on the Hartree-Fock basis state. For a
    new ordering, only the sign of that pair changes according to whether i or
    j appears first. Equal output basis states are combined with NumPy bincount.
    """

    number_of_pauli_terms: int
    left_indices: np.ndarray
    right_indices: np.ndarray
    target_bins: np.ndarray
    pair_amplitudes: np.ndarray
    number_of_target_bins: int

    @classmethod
    def build(
        cls,
        pauli_keys: Sequence[PauliKey],
        coefficients: dict[PauliKey, complex],
        pauli_graph: nx.Graph,
        n_qubits: int,
        n_electrons: int,
        tolerance: float,
    ) -> "HFCommutatorEvaluator":
        real_coefficients = np.asarray(
            [real_coefficient(coefficients[key], tolerance) for key in pauli_keys],
            dtype=np.float64,
        )
        hf_index = hartree_fock_basis_index(n_qubits, n_electrons)

        left_values: list[int] = []
        right_values: list[int] = []
        targets: list[int] = []
        amplitudes: list[complex] = []

        for first, second in pauli_graph.edges:
            left = min(int(first), int(second))
            right = max(int(first), int(second))
            product_phase, product_key = multiply_pauli_keys(
                pauli_keys[left],
                pauli_keys[right],
            )
            target, pauli_phase = apply_pauli_to_basis_index(
                product_key,
                hf_index,
                n_qubits,
            )
            amplitude = (
                2.0
                * real_coefficients[left]
                * real_coefficients[right]
                * product_phase
                * pauli_phase
            )
            if abs(amplitude) <= tolerance:
                continue
            left_values.append(left)
            right_values.append(right)
            targets.append(target)
            amplitudes.append(amplitude)

        if not amplitudes:
            return cls(
                number_of_pauli_terms=len(pauli_keys),
                left_indices=np.empty(0, dtype=np.int32),
                right_indices=np.empty(0, dtype=np.int32),
                target_bins=np.empty(0, dtype=np.int32),
                pair_amplitudes=np.empty(0, dtype=np.complex128),
                number_of_target_bins=0,
            )

        _, target_bins = np.unique(
            np.asarray(targets, dtype=np.int64),
            return_inverse=True,
        )
        return cls(
            number_of_pauli_terms=len(pauli_keys),
            left_indices=np.asarray(left_values, dtype=np.int32),
            right_indices=np.asarray(right_values, dtype=np.int32),
            target_bins=np.asarray(target_bins, dtype=np.int32),
            pair_amplitudes=np.asarray(amplitudes, dtype=np.complex128),
            number_of_target_bins=int(np.max(target_bins)) + 1,
        )

    def evaluate(self, pauli_order_indices: Sequence[int]) -> float:
        if len(pauli_order_indices) != self.number_of_pauli_terms:
            raise ValueError(
                "Pauli ordering length does not match the evaluator: "
                f"{len(pauli_order_indices)} != {self.number_of_pauli_terms}."
            )
        if self.pair_amplitudes.size == 0:
            return 0.0

        positions = np.empty(self.number_of_pauli_terms, dtype=np.int32)
        positions[np.asarray(pauli_order_indices, dtype=np.int32)] = np.arange(
            self.number_of_pauli_terms,
            dtype=np.int32,
        )
        signs = np.where(
            positions[self.left_indices] < positions[self.right_indices],
            1.0,
            -1.0,
        )
        signed_amplitudes = self.pair_amplitudes * signs
        real_parts = np.bincount(
            self.target_bins,
            weights=signed_amplitudes.real,
            minlength=self.number_of_target_bins,
        )
        imaginary_parts = np.bincount(
            self.target_bins,
            weights=signed_amplitudes.imag,
            minlength=self.number_of_target_bins,
        )
        return float(
            math.sqrt(
                float(
                    np.dot(real_parts, real_parts)
                    + np.dot(imaginary_parts, imaginary_parts)
                )
            )
        )


_SINGLE_QUBIT_PRODUCT: dict[tuple[str, str], tuple[complex, str | None]] = {
    ("X", "X"): (1.0 + 0.0j, None),
    ("Y", "Y"): (1.0 + 0.0j, None),
    ("Z", "Z"): (1.0 + 0.0j, None),
    ("X", "Y"): (1.0j, "Z"),
    ("Y", "X"): (-1.0j, "Z"),
    ("Y", "Z"): (1.0j, "X"),
    ("Z", "Y"): (-1.0j, "X"),
    ("Z", "X"): (1.0j, "Y"),
    ("X", "Z"): (-1.0j, "Y"),
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Test current, reversed, and randomized JW- and fermionic-color "
            "schedules for B2/STO-6G active spaces."
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
        "--samples",
        type=int,
        default=100,
        help=(
            "Random samples for each random schedule family and graph level. "
            "Default: 100. Use a smaller value for a quick test."
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=125,
        help="Base random seed. Default: 125.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/b2_coloring_robustness_results.csv"),
        help="Detailed result CSV.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("analysis/b2_coloring_robustness_summary.csv"),
        help="Aggregated result CSV.",
    )
    parser.add_argument(
        "--plot-directory",
        type=Path,
        default=Path("analysis/b2_coloring_robustness_figures"),
        help="Directory for BCH-versus-infidelity plots.",
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
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Append to the CSV and skip completed schedules.",
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
    if args.samples < 0:
        raise ValueError("--samples cannot be negative.")
    if args.tolerance <= 0.0:
        raise ValueError("--tolerance must be positive.")
    if args.parallel_threshold <= 0:
        raise ValueError("--parallel-threshold must be positive.")


def select_cases(
    library: Path,
    requested_qubits: Sequence[int],
    bond_length: str,
) -> list[Path]:
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


def deterministic_color_groups(graph: nx.Graph) -> list[list[int]]:
    """Reproduce QHAT's largest-first coloring and deterministic group order."""
    coloring = nx.greedy_color(graph, strategy="largest_first")
    if not coloring:
        return []
    return [
        [int(node) for node in graph.nodes if coloring[node] == color]
        for color in sorted(set(coloring.values()))
    ]


def graph_density(graph: nx.Graph) -> float:
    vertices = graph.number_of_nodes()
    if vertices < 2:
        return 0.0
    return 2.0 * graph.number_of_edges() / (vertices * (vertices - 1))


def flatten(groups: Sequence[Sequence[int]]) -> list[int]:
    return [node for group in groups for node in group]


def copy_and_shuffle(values: Sequence[int], rng: np.random.Generator) -> list[int]:
    result = list(values)
    rng.shuffle(result)
    return result


def schedule_seed(
    base_seed: int,
    n_qubits: int,
    graph_level: str,
    schedule: str,
    sample_index: int,
) -> int:
    graph_code = {"jw": 1, "fermionic": 2}[graph_level]
    schedule_code = {
        "random_groups": 1,
        "random_within": 2,
        "random_both": 3,
    }[schedule]
    sequence = np.random.SeedSequence(
        [base_seed, n_qubits, graph_code, schedule_code, sample_index]
    )
    return int(sequence.generate_state(1, dtype=np.uint32)[0])


def make_schedule_spec(
    graph_level: str,
    schedule: str,
    sample_index: int,
    random_seed: int,
    base_groups: Sequence[Sequence[int]],
) -> ScheduleSpec:
    group_order = list(range(len(base_groups)))
    grouped_nodes = [list(group) for group in base_groups]

    if schedule == "current":
        pass
    elif schedule == "reverse_groups":
        group_order.reverse()
        grouped_nodes = [grouped_nodes[index] for index in group_order]
    elif schedule in RANDOM_SCHEDULES:
        rng = np.random.default_rng(random_seed)
        if schedule in {"random_groups", "random_both"}:
            rng.shuffle(group_order)
        if schedule in {"random_within", "random_both"}:
            grouped_nodes = [
                copy_and_shuffle(group, rng)
                for group in grouped_nodes
            ]
        grouped_nodes = [grouped_nodes[index] for index in group_order]
    else:
        raise ValueError(f"Unsupported schedule {schedule!r}.")

    return ScheduleSpec(
        graph_level=graph_level,
        schedule=schedule,
        sample_index=sample_index,
        random_seed=random_seed,
        node_order=tuple(flatten(grouped_nodes)),
        group_order=tuple(group_order),
        grouped_nodes=tuple(tuple(group) for group in grouped_nodes),
    )


def build_schedule_specs(
    graph_level: str,
    groups: Sequence[Sequence[int]],
    n_qubits: int,
    samples: int,
    base_seed: int,
) -> list[ScheduleSpec]:
    specs = [
        make_schedule_spec(
            graph_level,
            "current",
            -1,
            -1,
            groups,
        ),
        make_schedule_spec(
            graph_level,
            "reverse_groups",
            -1,
            -1,
            groups,
        ),
    ]
    for schedule in RANDOM_SCHEDULES:
        for sample_index in range(samples):
            seed = schedule_seed(
                base_seed,
                n_qubits,
                graph_level,
                schedule,
                sample_index,
            )
            specs.append(
                make_schedule_spec(
                    graph_level,
                    schedule,
                    sample_index,
                    seed,
                    groups,
                )
            )
    return specs


def multiply_pauli_keys(
    left: PauliKey,
    right: PauliKey,
) -> tuple[complex, PauliKey]:
    left_map = dict(left)
    right_map = dict(right)
    phase = 1.0 + 0.0j
    product: list[tuple[int, str]] = []

    for qubit in sorted(left_map.keys() | right_map.keys()):
        left_label = left_map.get(qubit)
        right_label = right_map.get(qubit)
        if left_label is None:
            if right_label is None:
                raise RuntimeError("Invalid empty local Pauli product.")
            product.append((qubit, right_label))
        elif right_label is None:
            product.append((qubit, left_label))
        else:
            local_phase, local_label = _SINGLE_QUBIT_PRODUCT[
                (left_label, right_label)
            ]
            phase *= local_phase
            if local_label is not None:
                product.append((qubit, local_label))

    return phase, tuple(product)


def apply_pauli_to_basis_index(
    pauli_key: PauliKey,
    basis_index: int,
    n_qubits: int,
) -> tuple[int, complex]:
    target = basis_index
    phase = 1.0 + 0.0j

    for qubit, label in pauli_key:
        bit_position = n_qubits - 1 - qubit
        bit_mask = 1 << bit_position
        occupied = bool(target & bit_mask)

        if label == "X":
            target ^= bit_mask
        elif label == "Y":
            phase *= -1.0j if occupied else 1.0j
            target ^= bit_mask
        elif label == "Z":
            if occupied:
                phase = -phase
        else:
            raise ValueError(f"Unsupported Pauli label {label!r}.")

    return target, phase


def hartree_fock_basis_index(n_qubits: int, n_electrons: int) -> int:
    return sum(
        1 << (n_qubits - 1 - qubit)
        for qubit in range(n_electrons)
    )


def precompute_fermion_to_pauli_indices(
    hermitian_terms: Sequence[Any],
    final_coefficients: dict[PauliKey, complex],
    raw_index_by_key: dict[PauliKey, int],
    tolerance: float,
) -> list[tuple[int, ...]]:
    mapped_indices: list[tuple[int, ...]] = []
    for fermionic_term in hermitian_terms:
        mapped = jordan_wigner(fermionic_term.operator)
        mapped.compress(abs_tol=tolerance)
        indices: list[int] = []
        for pauli_key, coefficient in mapped.terms.items():
            if pauli_key == () or abs(coefficient) <= tolerance:
                continue
            if pauli_key not in final_coefficients:
                continue
            indices.append(raw_index_by_key[pauli_key])
        mapped_indices.append(tuple(indices))
    return mapped_indices


def induced_pauli_order_indices(
    fermionic_node_order: Sequence[int],
    fermion_to_pauli_indices: Sequence[Sequence[int]],
    number_of_pauli_terms: int,
) -> list[int]:
    result: list[int] = []
    seen = np.zeros(number_of_pauli_terms, dtype=np.bool_)

    for node in fermionic_node_order:
        for pauli_index in fermion_to_pauli_indices[node]:
            if not seen[pauli_index]:
                seen[pauli_index] = True
                result.append(int(pauli_index))

    for pauli_index in range(number_of_pauli_terms):
        if not seen[pauli_index]:
            seen[pauli_index] = True
            result.append(pauli_index)

    if len(result) != number_of_pauli_terms or len(set(result)) != len(result):
        raise RuntimeError("Induced fermionic Pauli order is not a permutation.")
    return result


def hash_nested_integers(values: Sequence[Sequence[int]]) -> str:
    payload = json.dumps([list(group) for group in values], separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def hash_integer_order(values: Sequence[int]) -> str:
    payload = ",".join(str(value) for value in values)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def safe_ratio(numerator: float, denominator: float) -> str | float:
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        return ""
    if denominator <= 0.0:
        return ""
    return numerator / denominator


def blank_row() -> dict[str, Any]:
    return {field: "" for field in FIELDNAMES}


def result_key(
    case_id: str,
    graph_level: str,
    schedule: str,
    sample_index: int,
    random_seed: int,
    trotter_steps: int,
    evolution_time: float,
) -> tuple[str, str, str, int, int, int, float]:
    return (
        case_id,
        graph_level,
        schedule,
        sample_index,
        random_seed,
        trotter_steps,
        evolution_time,
    )


def load_resume_data(
    output_path: Path,
) -> tuple[
    set[tuple[str, str, str, int, int, int, float]],
    dict[str, float],
    dict[tuple[str, str], float],
    dict[str, float],
]:
    completed: set[tuple[str, str, str, int, int, int, float]] = set()
    raw_infidelities: dict[str, float] = {}
    current_infidelities: dict[tuple[str, str], float] = {}
    raw_bch_norms: dict[str, float] = {}

    if not output_path.exists():
        return completed, raw_infidelities, current_infidelities, raw_bch_norms

    with output_path.open("r", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row.get("status") != "success":
                continue
            try:
                key = result_key(
                    row["case_id"],
                    row["graph_level"],
                    row["schedule"],
                    int(row["sample_index"]),
                    int(row["random_seed"]),
                    int(row["trotter_steps"]),
                    float(row["evolution_time"]),
                )
                infidelity = float(row["state_infidelity"])
                bch_norm = float(row["bch2_hf_state_norm"])
            except (KeyError, TypeError, ValueError):
                continue

            completed.add(key)
            if row["schedule"] == "jw_raw":
                raw_infidelities[row["case_id"]] = infidelity
                raw_bch_norms[row["case_id"]] = bch_norm
            elif row["schedule"] == "current":
                current_infidelities[(row["case_id"], row["graph_level"])] = (
                    infidelity
                )

    return completed, raw_infidelities, current_infidelities, raw_bch_norms


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


def build_case_schedules(
    n_qubits: int,
    samples: int,
    base_seed: int,
    raw_pauli_keys: Sequence[PauliKey],
    final_coefficients: dict[PauliKey, complex],
    hermitian_terms: Sequence[Any],
    jw_graph: nx.Graph,
    fermionic_graph: nx.Graph,
    tolerance: float,
) -> tuple[
    list[tuple[ScheduleSpec | None, list[int]]],
    dict[str, list[list[int]]],
]:
    raw_index_by_key = {
        key: index for index, key in enumerate(raw_pauli_keys)
    }
    fermion_to_pauli_indices = precompute_fermion_to_pauli_indices(
        hermitian_terms,
        final_coefficients,
        raw_index_by_key,
        tolerance,
    )

    jw_groups = deterministic_color_groups(jw_graph)
    fermionic_groups = deterministic_color_groups(fermionic_graph)
    groups_by_level = {
        "jw": jw_groups,
        "fermionic": fermionic_groups,
    }

    schedules: list[tuple[ScheduleSpec | None, list[int]]] = [
        (None, list(range(len(raw_pauli_keys))))
    ]

    # Current schedules are evaluated first, then reversed schedules, then
    # random families. This makes raw/current baselines available immediately.
    all_specs = {
        level: build_schedule_specs(
            level,
            groups,
            n_qubits,
            samples,
            base_seed,
        )
        for level, groups in groups_by_level.items()
    }

    ordered_specs: list[ScheduleSpec] = []
    for schedule_name in ("current", "reverse_groups"):
        for level in ("jw", "fermionic"):
            ordered_specs.extend(
                spec
                for spec in all_specs[level]
                if spec.schedule == schedule_name
            )
    for schedule_name in RANDOM_SCHEDULES:
        for sample_index in range(samples):
            for level in ("jw", "fermionic"):
                ordered_specs.extend(
                    spec
                    for spec in all_specs[level]
                    if spec.schedule == schedule_name
                    and spec.sample_index == sample_index
                )

    for spec in ordered_specs:
        if spec.graph_level == "jw":
            pauli_order_indices = list(spec.node_order)
        else:
            pauli_order_indices = induced_pauli_order_indices(
                spec.node_order,
                fermion_to_pauli_indices,
                len(raw_pauli_keys),
            )
        schedules.append((spec, pauli_order_indices))

    return schedules, groups_by_level


def benchmark_case(
    tensor_path: Path,
    args: argparse.Namespace,
    completed: set[tuple[str, str, str, int, int, int, float]],
    raw_infidelities: dict[str, float],
    current_infidelities: dict[tuple[str, str], float],
    raw_bch_norms: dict[str, float],
    writer: csv.DictWriter,
    output_stream: Any,
) -> None:
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
    hermitian_terms = build_hermitian_fermion_terms(
        fermion_hamiltonian,
        args.tolerance,
    )

    print("  Building fixed graph partitions...")
    jw_graph, _ = build_pauli_noncommutation_graph(raw_pauli_keys)
    fermionic_graph, _ = build_fermionic_noncommutation_graph(
        hermitian_terms,
        args.tolerance,
    )

    schedules, groups_by_level = build_case_schedules(
        n_qubits=n_qubits,
        samples=args.samples,
        base_seed=args.seed,
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        hermitian_terms=hermitian_terms,
        jw_graph=jw_graph,
        fermionic_graph=fermionic_graph,
        tolerance=args.tolerance,
    )

    # Confirm that the reconstructed current node schedules exactly match the
    # branch's deterministic coloring rule before randomization begins.
    current_jw_spec, current_jw_indices = schedules[1]
    current_fermion_spec, current_fermion_indices = schedules[2]
    if current_jw_spec is None or current_fermion_spec is None:
        raise RuntimeError("Current schedule reconstruction failed.")

    branch_jw_nodes, _, _ = deterministic_coloring_order(jw_graph)
    branch_fermion_nodes, _, _ = deterministic_coloring_order(fermionic_graph)
    if list(current_jw_spec.node_order) != branch_jw_nodes:
        raise RuntimeError("Reconstructed JW current schedule differs from QHAT.")
    if list(current_fermion_spec.node_order) != branch_fermion_nodes:
        raise RuntimeError(
            "Reconstructed fermionic current schedule differs from QHAT."
        )

    validate_pauli_order(
        "jw_current",
        [raw_pauli_keys[index] for index in current_jw_indices],
        raw_pauli_keys,
    )
    validate_pauli_order(
        "fermionic_current",
        [raw_pauli_keys[index] for index in current_fermion_indices],
        raw_pauli_keys,
    )

    print("  Precomputing state-dependent commutator actions...")
    bch_evaluator = HFCommutatorEvaluator.build(
        pauli_keys=raw_pauli_keys,
        coefficients=final_coefficients,
        pauli_graph=jw_graph,
        n_qubits=n_qubits,
        n_electrons=metadata.active_occupied,
        tolerance=args.tolerance,
    )

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

    total = len(schedules)
    for schedule_index, (spec, pauli_order_indices) in enumerate(
        schedules,
        start=1,
    ):
        if spec is None:
            graph_level = "baseline"
            schedule_name = "jw_raw"
            sample_index = -1
            random_seed_value = -1
            group_order: Sequence[int] = []
            grouped_nodes: Sequence[Sequence[int]] = []
            graph = None
            groups: Sequence[Sequence[int]] = []
        else:
            graph_level = spec.graph_level
            schedule_name = spec.schedule
            sample_index = spec.sample_index
            random_seed_value = spec.random_seed
            group_order = spec.group_order
            grouped_nodes = spec.grouped_nodes
            graph = graph_by_level[graph_level]
            groups = groups_by_level[graph_level]

        key = result_key(
            metadata.case_id,
            graph_level,
            schedule_name,
            sample_index,
            random_seed_value,
            args.steps,
            args.evolution_time,
        )
        if key in completed:
            print(
                f"  [{schedule_index}/{total}] SKIP "
                f"{graph_level}:{schedule_name}:{sample_index}"
            )
            continue

        print(
            f"  [{schedule_index}/{total}] RUN  "
            f"{graph_level}:{schedule_name}:{sample_index}"
        )
        row = blank_row()
        add_metadata(row, metadata)
        row.update(
            {
                "status": "success",
                "number_of_fermionic_terms": len(hermitian_terms),
                "number_of_pauli_terms": len(raw_pauli_keys),
                "graph_level": graph_level,
                "schedule": schedule_name,
                "sample_index": sample_index,
                "random_seed": random_seed_value,
                "graph_vertices": (
                    graph.number_of_nodes() if graph is not None else ""
                ),
                "graph_edges": (
                    graph.number_of_edges() if graph is not None else ""
                ),
                "graph_density": (
                    graph_density(graph) if graph is not None else ""
                ),
                "number_of_colors": len(groups) if graph is not None else "",
                "group_sizes": json.dumps(
                    [len(group) for group in groups],
                    separators=(",", ":"),
                ),
                "group_order": json.dumps(
                    list(group_order),
                    separators=(",", ":"),
                ),
                "within_group_order_hash": hash_nested_integers(grouped_nodes),
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
            schedule_build_start = time.perf_counter()
            ordered_terms = [
                compiled_raw_terms[index]
                for index in pauli_order_indices
            ]
            row["schedule_build_time_seconds"] = (
                time.perf_counter() - schedule_build_start
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

            if schedule_name == "jw_raw":
                raw_bch_norms[metadata.case_id] = bch_norm
                row["bch_squared_ratio_to_raw"] = 1.0
            else:
                raw_bch = raw_bch_norms.get(metadata.case_id)
                if raw_bch is None:
                    raise RuntimeError("Raw JW BCH baseline was not processed first.")
                ratio = safe_ratio(bch_norm, raw_bch)
                row["bch_squared_ratio_to_raw"] = (
                    ratio**2 if isinstance(ratio, float) else ""
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
            row.update(metrics)
            infidelity = metrics["state_infidelity"]

            if schedule_name == "jw_raw":
                raw_infidelities[metadata.case_id] = infidelity
                row["state_infidelity_ratio_to_raw_jw"] = 1.0
            else:
                raw_error = raw_infidelities.get(metadata.case_id)
                if raw_error is None:
                    raise RuntimeError("Raw JW error baseline was not processed first.")
                row["state_infidelity_ratio_to_raw_jw"] = safe_ratio(
                    infidelity,
                    raw_error,
                )

            if schedule_name == "current":
                current_infidelities[(metadata.case_id, graph_level)] = infidelity
                row["state_infidelity_ratio_to_current_coloring"] = 1.0
            elif graph_level in {"jw", "fermionic"}:
                current_error = current_infidelities.get(
                    (metadata.case_id, graph_level)
                )
                if current_error is None:
                    raise RuntimeError(
                        f"Current {graph_level} coloring was not processed first."
                    )
                row["state_infidelity_ratio_to_current_coloring"] = safe_ratio(
                    infidelity,
                    current_error,
                )

            writer.writerow(row)
            output_stream.flush()
            completed.add(key)
            print(
                "       "
                f"BCH_HF={bch_norm:.6e}  "
                f"infidelity={infidelity:.6e}  "
                f"ratio/raw={row['state_infidelity_ratio_to_raw_jw']}"
            )
        except Exception as error:
            row["status"] = "failed"
            row["error_message"] = f"{type(error).__name__}: {error}"
            writer.writerow(row)
            output_stream.flush()
            print(f"       FAILED: {type(error).__name__}: {error}")


def read_successful_rows(output_path: Path) -> list[dict[str, str]]:
    if not output_path.exists():
        return []
    with output_path.open("r", newline="", encoding="utf-8") as stream:
        return [
            row
            for row in csv.DictReader(stream)
            if row.get("status") == "success"
        ]


def finite_float(row: dict[str, str], field: str) -> float | None:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def population_statistics(values: Sequence[float]) -> dict[str, float | str]:
    if not values:
        return {
            "minimum": "",
            "median": "",
            "maximum": "",
            "mean": "",
            "std": "",
        }
    array = np.asarray(values, dtype=np.float64)
    return {
        "minimum": float(np.min(array)),
        "median": float(np.median(array)),
        "maximum": float(np.max(array)),
        "mean": float(np.mean(array)),
        "std": float(np.std(array)),
    }


def log_pearson_correlation(
    bch_values: Sequence[float],
    infidelities: Sequence[float],
) -> str | float:
    paired = [
        (bch, infidelity)
        for bch, infidelity in zip(bch_values, infidelities)
        if bch > 0.0 and infidelity > 0.0
    ]
    if len(paired) < 3:
        return ""
    x_values = np.log10([bch**2 for bch, _ in paired])
    y_values = np.log10([infidelity for _, infidelity in paired])
    if float(np.std(x_values)) == 0.0 or float(np.std(y_values)) == 0.0:
        return ""
    return float(np.corrcoef(x_values, y_values)[0, 1])


def write_summary(output_path: Path, summary_path: Path) -> None:
    rows = read_successful_rows(output_path)
    if not rows:
        return

    raw_by_case = {
        row["case_id"]: float(row["state_infidelity"])
        for row in rows
        if row["schedule"] == "jw_raw"
    }
    current_by_case_level = {
        (row["case_id"], row["graph_level"]): float(row["state_infidelity"])
        for row in rows
        if row["schedule"] == "current"
    }

    grouped: dict[tuple[str, int, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (
            row["case_id"],
            int(row["n_qubits"]),
            row["graph_level"],
            row["schedule"],
        )
        grouped[key].append(row)

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()

        for key in sorted(grouped, key=lambda item: (item[1], item[2], item[3])):
            case_id, n_qubits, graph_level, schedule = key
            group_rows = grouped[key]
            bch_values = [
                value
                for row in group_rows
                if (value := finite_float(row, "bch2_hf_state_norm")) is not None
            ]
            infidelities = [
                value
                for row in group_rows
                if (value := finite_float(row, "state_infidelity")) is not None
            ]
            bch_stats = population_statistics(bch_values)
            error_stats = population_statistics(infidelities)
            raw_error = raw_by_case.get(case_id)
            current_error = current_by_case_level.get((case_id, graph_level))

            writer.writerow(
                {
                    "case_id": case_id,
                    "n_qubits": n_qubits,
                    "graph_level": graph_level,
                    "schedule": schedule,
                    "number_of_samples": len(group_rows),
                    "minimum_bch2_hf_state_norm": bch_stats["minimum"],
                    "median_bch2_hf_state_norm": bch_stats["median"],
                    "maximum_bch2_hf_state_norm": bch_stats["maximum"],
                    "mean_bch2_hf_state_norm": bch_stats["mean"],
                    "std_bch2_hf_state_norm": bch_stats["std"],
                    "minimum_state_infidelity": error_stats["minimum"],
                    "median_state_infidelity": error_stats["median"],
                    "maximum_state_infidelity": error_stats["maximum"],
                    "mean_state_infidelity": error_stats["mean"],
                    "std_state_infidelity": error_stats["std"],
                    "fraction_beating_raw_jw": (
                        float(np.mean(np.asarray(infidelities) < raw_error))
                        if infidelities and raw_error is not None
                        else ""
                    ),
                    "fraction_beating_current_coloring": (
                        float(np.mean(np.asarray(infidelities) < current_error))
                        if infidelities and current_error is not None
                        else ""
                    ),
                    "pearson_correlation_log_bch2_squared_vs_log_infidelity": (
                        log_pearson_correlation(bch_values, infidelities)
                        if len(bch_values) == len(infidelities)
                        else ""
                    ),
                }
            )


def write_plots(output_path: Path, plot_directory: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is unavailable; skipping plots.")
        return

    rows = read_successful_rows(output_path)
    if not rows:
        return
    plot_directory.mkdir(parents=True, exist_ok=True)

    qubit_counts = sorted({int(row["n_qubits"]) for row in rows})
    for n_qubits in qubit_counts:
        selected = [row for row in rows if int(row["n_qubits"]) == n_qubits]
        figure, axis = plt.subplots(figsize=(8.5, 5.5))
        plotted = False

        labels = sorted(
            {
                f"{row['graph_level']}:{row['schedule']}"
                for row in selected
                if float(row["bch2_hf_state_norm"]) > 0.0
                and float(row["state_infidelity"]) > 0.0
            }
        )
        for label in labels:
            x_values: list[float] = []
            y_values: list[float] = []
            for row in selected:
                row_label = f"{row['graph_level']}:{row['schedule']}"
                if row_label != label:
                    continue
                bch = float(row["bch2_hf_state_norm"])
                error = float(row["state_infidelity"])
                if bch <= 0.0 or error <= 0.0:
                    continue
                x_values.append(bch**2)
                y_values.append(error)
            if x_values:
                axis.scatter(x_values, y_values, label=label, alpha=0.75)
                plotted = True

        if not plotted:
            plt.close(figure)
            continue
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel(r"$\|C_2|\psi_{HF}\rangle\|^2$")
        axis.set_ylabel("State infidelity")
        axis.set_title(
            f"B2/STO-6G {n_qubits}-qubit coloring robustness\n"
            "first order, fixed time and Trotter steps"
        )
        axis.grid(True, which="both", alpha=0.3)
        axis.legend(fontsize="small", ncol=2)
        figure.tight_layout()
        figure.savefig(
            plot_directory / f"b2_{n_qubits}q_bch_vs_infidelity.png",
            dpi=200,
        )
        plt.close(figure)


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    warm_up_numba()

    tensor_paths = select_cases(
        args.library,
        args.qubits,
        args.bond_length,
    )
    schedules_per_case = 1 + 2 * (2 + 3 * args.samples)

    print(f"B2 cases selected: {len(tensor_paths)}")
    print(f"Active-space qubits: {sorted(args.qubits)}")
    print(f"Random samples per family: {args.samples}")
    print(f"Schedules per case: {schedules_per_case}")
    print(f"First-order Trotter steps: {args.steps}")
    print(f"Evolution time: {args.evolution_time}")
    print(f"Output: {args.output}")
    print(
        "Note: 18-qubit random schedules can be expensive. For a smoke test, "
        "run first with --samples 2 --qubits 12."
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.resume:
        (
            completed,
            raw_infidelities,
            current_infidelities,
            raw_bch_norms,
        ) = load_resume_data(args.output)
        file_mode = "a"
        write_header = (
            not args.output.exists() or args.output.stat().st_size == 0
        )
        print(f"Completed schedules loaded: {len(completed)}")
    else:
        completed = set()
        raw_infidelities = {}
        current_infidelities = {}
        raw_bch_norms = {}
        file_mode = "w"
        write_header = True

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
            print("=" * 96)
            print(f"[{case_index}/{len(tensor_paths)}] {tensor_path}")
            try:
                benchmark_case(
                    tensor_path=tensor_path,
                    args=args,
                    completed=completed,
                    raw_infidelities=raw_infidelities,
                    current_infidelities=current_infidelities,
                    raw_bch_norms=raw_bch_norms,
                    writer=writer,
                    output_stream=output_stream,
                )
            except Exception as error:
                print(
                    f"CASE FAILED {tensor_path}: "
                    f"{type(error).__name__}: {error}"
                )
                failure = blank_row()
                failure.update(
                    {
                        "status": "failed",
                        "error_message": f"{type(error).__name__}: {error}",
                        "case_id": tensor_path.name.removesuffix(
                            ".tensors.npz"
                        ),
                        "tensor_path": str(tensor_path),
                        "trotter_steps": args.steps,
                        "evolution_time": args.evolution_time,
                        "coefficient_tolerance": args.tolerance,
                    }
                )
                writer.writerow(failure)
                output_stream.flush()

    write_summary(args.output, args.summary)
    write_plots(args.output, args.plot_directory)

    print()
    print("=" * 96)
    print(f"Detailed results: {args.output}")
    print(f"Summary:          {args.summary}")
    print(f"Figures:          {args.plot_directory}")
    print("=" * 96)


if __name__ == "__main__":
    main()

