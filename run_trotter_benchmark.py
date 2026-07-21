#!/usr/bin/env python3
"""Run focused fermionic-vs-JW Trotter ordering benchmarks.

The benchmark uses three representative QHAT Hamiltonians and compares two
decompositions of the same Hamiltonian:

* Hermitian fermionic terms reconstructed from QHAT tensors, then mapped to JW.
* Individual JW Pauli strings read from QHAT's ``.dat`` output.

For each decomposition it records the noncommutation graph, greedy coloring,
color-block commutator norms, colored ordering error, random vertex-ordering
errors, and particle-number leakage from the active-space Hartree--Fock state.
Only first- and second-order product formulas are implemented intentionally.
"""

from __future__ import annotations

import argparse
import json
import zlib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import networkx as nx
import numpy as np
import pandas as pd
from scipy.linalg import expm, norm


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}


@dataclass(frozen=True)
class CaseSpec:
    case_id: str
    molecule: str
    bond_length: float
    basis: str
    active_occupied: int
    active_vacant: int
    dat_path: str
    tensor_path: str
    initial_state_path: str
    selection_reason: str


DEFAULT_CASES = (
    CaseSpec(
        case_id="similar_4q",
        molecule="Be-Be",
        bond_length=2.244,
        basis="sto-6g",
        active_occupied=2,
        active_vacant=2,
        dat_path=(
            "demo_L_sweep_library/Be-Be/2.24/sto-6g/"
            "Be-Be_2.24_sto-6g_as-002-002_jw.dat"
        ),
        tensor_path=(
            "demo_L_sweep_library/Be-Be/2.24/sto-6g/"
            "Be-Be_2.24_sto-6g_as-002-002.tensors.npz"
        ),
        initial_state_path=(
            "demo_L_sweep_library/Be-Be/2.24/sto-6g/"
            "Be-Be_2.24_sto-6g_as-002-002_jw.npy"
        ),
        selection_reason="4-qubit case where both colored methods have similar error",
    ),
    CaseSpec(
        case_id="fermionic_first_order_6q",
        molecule="B-B",
        bond_length=2.38,
        basis="sto-6g",
        active_occupied=2,
        active_vacant=4,
        dat_path=(
            "demo_L_sweep_library/B-B/2.38/sto-6g/"
            "B-B_2.38_sto-6g_as-002-004_jw.dat"
        ),
        tensor_path=(
            "demo_L_sweep_library/B-B/2.38/sto-6g/"
            "B-B_2.38_sto-6g_as-002-004.tensors.npz"
        ),
        initial_state_path=(
            "demo_L_sweep_library/B-B/2.38/sto-6g/"
            "B-B_2.38_sto-6g_as-002-004_jw.npy"
        ),
        selection_reason="6-qubit case where fermionic coloring is better at first order",
    ),
    CaseSpec(
        case_id="jw_second_order_6q",
        molecule="B-B",
        bond_length=2.38,
        basis="hgbs-5",
        active_occupied=2,
        active_vacant=4,
        dat_path=(
            "demo_L_sweep_library/B-B/2.38/hgbs-5/"
            "B-B_2.38_hgbs-5_as-002-004_jw.dat"
        ),
        tensor_path=(
            "demo_L_sweep_library/B-B/2.38/hgbs-5/"
            "B-B_2.38_hgbs-5_as-002-004.tensors.npz"
        ),
        initial_state_path=(
            "demo_L_sweep_library/B-B/2.38/hgbs-5/"
            "B-B_2.38_hgbs-5_as-002-004_jw.npy"
        ),
        selection_reason="6-qubit case where direct JW coloring is much better at second order",
    ),
)


def parse_qhat_dat(path: Path) -> tuple[dict[str, str], list[tuple[str, complex]]]:
    metadata: dict[str, str] = {}
    terms: list[tuple[str, complex]] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                if "=" in line:
                    key, value = line[1:].split("=", 1)
                    metadata[key.strip()] = value.strip()
                continue
            coefficient, pauli = line.split()
            terms.append((pauli, complex(coefficient)))
    if not terms:
        raise ValueError(f"No Pauli terms found in {path}")
    return metadata, terms


def kron_all(items: Iterable[np.ndarray]) -> np.ndarray:
    result = np.array([[1]], dtype=complex)
    for item in items:
        result = np.kron(result, item)
    return result


def pauli_matrix(word: str) -> np.ndarray:
    return kron_all(PAULI[letter] for letter in word)


def paulis_commute(left: str, right: str) -> bool:
    mismatches = sum(
        a != "I" and b != "I" and a != b for a, b in zip(left, right)
    )
    return mismatches % 2 == 0


def ladder_matrices(n_qubits: int) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Return JW annihilation and creation matrices in QHAT string order."""
    annihilation: list[np.ndarray] = []
    creation: list[np.ndarray] = []
    for target in range(n_qubits):
        factors_a: list[np.ndarray] = []
        factors_c: list[np.ndarray] = []
        for qubit in range(n_qubits):
            if qubit < target:
                factors_a.append(Z)
                factors_c.append(Z)
            elif qubit == target:
                factors_a.append((X + 1j * Y) / 2)
                factors_c.append((X - 1j * Y) / 2)
            else:
                factors_a.append(I2)
                factors_c.append(I2)
        annihilation.append(kron_all(factors_a))
        creation.append(kron_all(factors_c))
    return annihilation, creation


def monomial_matrix(
    key: tuple[tuple[int, int], ...],
    annihilation: Sequence[np.ndarray],
    creation: Sequence[np.ndarray],
) -> np.ndarray:
    result = np.eye(annihilation[0].shape[0], dtype=complex)
    for mode, action in key:
        result = result @ (creation[mode] if action else annihilation[mode])
    return result


def canonicalize_tensor_monomial(
    key: tuple[tuple[int, int], ...],
) -> tuple[tuple[tuple[int, int], ...] | None, int]:
    """Canonicalize an already normal-ordered QHAT tensor monomial.

    OpenFermion's ``normal_ordered`` convention places creation and annihilation
    operators in descending mode order. QHAT tensor monomials already have every
    creation operator before every annihilation operator, so only antisymmetric
    sorting and combining are required; no contraction terms can be introduced.
    Repeated equal-action operators make the monomial identically zero.
    """
    actions = [action for _, action in key]
    creation_count = sum(action == 1 for action in actions)
    if actions != [1] * creation_count + [0] * (len(actions) - creation_count):
        raise ValueError(f"Tensor monomial is not normal ordered: {key}")

    creators = [mode for mode, action in key if action == 1]
    annihilators = [mode for mode, action in key if action == 0]
    if len(set(creators)) != len(creators) or len(set(annihilators)) != len(
        annihilators
    ):
        return None, 0

    descending_swaps = sum(
        values[left] < values[right]
        for values in (creators, annihilators)
        for left in range(len(values))
        for right in range(left + 1, len(values))
    )
    canonical = tuple((mode, 1) for mode in sorted(creators, reverse=True)) + tuple(
        (mode, 0) for mode in sorted(annihilators, reverse=True)
    )
    return canonical, -1 if descending_swaps % 2 else 1


def dagger_key(key: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    raw_dagger = tuple((mode, 1 - action) for mode, action in reversed(key))
    canonical, sign = canonicalize_tensor_monomial(raw_dagger)
    if canonical is None or sign != 1:
        raise ValueError(f"Unexpected dagger canonicalization for {key}")
    return canonical


def tensor_monomials(
    tensor_path: Path, tolerance: float = 1e-12
) -> tuple[dict[tuple[tuple[int, int], ...], complex], int]:
    with np.load(tensor_path) as data:
        constant = complex(data["constant"]) if "constant" in data.files else 0.0
        one_body = data["one_body"]
        two_body = data["two_body"]

    coefficients: dict[tuple[tuple[int, int], ...], complex] = {}

    def accumulate(key: tuple[tuple[int, int], ...], coefficient: complex) -> None:
        canonical, sign = canonicalize_tensor_monomial(key)
        if canonical is not None:
            coefficients[canonical] = coefficients.get(canonical, 0.0) + (
                sign * coefficient
            )

    accumulate((), constant)
    for p, q in np.argwhere(np.abs(one_body) > tolerance):
        accumulate(((int(p), 1), (int(q), 0)), complex(one_body[p, q]))
    for p, q, r, s in np.argwhere(np.abs(two_body) > tolerance):
        accumulate(
            ((int(p), 1), (int(q), 1), (int(r), 0), (int(s), 0)),
            complex(two_body[p, q, r, s]),
        )
    combined = {
        key: coefficient
        for key, coefficient in coefficients.items()
        if abs(coefficient) > tolerance
    }
    return combined, one_body.shape[0]


def format_monomial(key: tuple[tuple[int, int], ...]) -> str:
    if not key:
        return "I"
    return " ".join(f"a{mode}{'^' if action else ''}" for mode, action in key)


def hermitian_vertices_from_tensors(
    tensor_path: Path, tolerance: float = 1e-12
) -> tuple[list[np.ndarray], list[str], int]:
    coefficients, n_qubits = tensor_monomials(tensor_path, tolerance)
    annihilation, creation = ladder_matrices(n_qubits)
    identity = np.eye(2**n_qubits, dtype=complex)
    used: set[tuple[tuple[int, int], ...]] = set()
    vertices: list[np.ndarray] = []
    labels: list[str] = []

    for key, coefficient in coefficients.items():
        if key in used:
            continue
        partner = dagger_key(key)
        matrix = (
            monomial_matrix(key, annihilation, creation) if key else identity.copy()
        )
        if partner == key:
            vertex = coefficient * matrix
            label = f"{coefficient.real:+.8g} {format_monomial(key)}"
            used.add(key)
        else:
            partner_coefficient = coefficients.get(partner, 0.0)
            partner_matrix = monomial_matrix(partner, annihilation, creation)
            vertex = coefficient * matrix + partner_coefficient * partner_matrix
            label = (
                f"{coefficient.real:+.8g} {format_monomial(key)} + "
                f"{partner_coefficient.real:+.8g} {format_monomial(partner)}"
            )
            used.update((key, partner))
        if norm(vertex, 2) > tolerance:
            if norm(vertex - vertex.conj().T, 2) > 1e-10:
                raise ValueError(f"Non-Hermitian reconstructed vertex: {label}")
            vertices.append(vertex)
            labels.append(label)
    return vertices, labels, n_qubits


def matrix_noncommutation_graph(
    vertices: Sequence[np.ndarray], tolerance: float = 1e-10
) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(vertices)))
    for left in range(len(vertices)):
        for right in range(left + 1, len(vertices)):
            commutator = vertices[left] @ vertices[right] - vertices[right] @ vertices[left]
            denominator = 2 * norm(vertices[left], 2) * norm(vertices[right], 2)
            normalized = norm(commutator, 2) / denominator if denominator else 0.0
            if normalized > tolerance:
                graph.add_edge(left, right)
    return graph


def pauli_noncommutation_graph(words: Sequence[str]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(words)))
    for left in range(len(words)):
        for right in range(left + 1, len(words)):
            if not paulis_commute(words[left], words[right]):
                graph.add_edge(left, right)
    return graph


def greedy_groups(graph: nx.Graph) -> tuple[dict[int, int], dict[int, list[int]]]:
    # Explicit tie-breaking keeps results stable across NetworkX versions.
    node_order = sorted(graph.nodes, key=lambda node: (-graph.degree[node], node))
    coloring: dict[int, int] = {}
    for node in node_order:
        neighbor_colors = {
            coloring[neighbor] for neighbor in graph.neighbors(node) if neighbor in coloring
        }
        color = 0
        while color in neighbor_colors:
            color += 1
        coloring[node] = color
    groups: dict[int, list[int]] = {}
    for node, color in coloring.items():
        groups.setdefault(color, []).append(node)
    for nodes in groups.values():
        nodes.sort()
    return coloring, dict(sorted(groups.items()))


def colored_vertex_order(groups: dict[int, list[int]]) -> list[int]:
    return [node for color in sorted(groups) for node in groups[color]]


def block_commutators(
    vertices: Sequence[np.ndarray], groups: dict[int, list[int]]
) -> tuple[list[np.ndarray], list[dict[str, float | int]], dict[str, float]]:
    zero = np.zeros_like(vertices[0])
    blocks = [sum((vertices[node] for node in groups[color]), zero.copy()) for color in groups]
    rows: list[dict[str, float | int]] = []
    ordered_sum = np.zeros_like(vertices[0])
    for left in range(len(blocks)):
        for right in range(left + 1, len(blocks)):
            commutator = blocks[left] @ blocks[right] - blocks[right] @ blocks[left]
            ordered_sum += commutator
            left_norm = float(norm(blocks[left], 2))
            right_norm = float(norm(blocks[right], 2))
            spectral = float(norm(commutator, 2))
            denominator = 2 * left_norm * right_norm
            left_nested = blocks[left] @ commutator - commutator @ blocks[left]
            right_nested = blocks[right] @ commutator - commutator @ blocks[right]
            rows.append(
                {
                    "left_color": left,
                    "right_color": right,
                    "left_block_spectral_norm": left_norm,
                    "right_block_spectral_norm": right_norm,
                    "commutator_spectral_norm": spectral,
                    "commutator_frobenius_norm": float(norm(commutator, "fro")),
                    "normalized_commutator_norm": spectral / denominator
                    if denominator
                    else 0.0,
                    "left_nested_commutator_norm": float(norm(left_nested, 2)),
                    "right_nested_commutator_norm": float(norm(right_nested, 2)),
                }
            )
    aggregates = {
        "ordered_commutator_sum_norm": float(norm(ordered_sum, 2)),
        "block_commutator_rss": float(
            np.sqrt(sum(float(row["commutator_spectral_norm"]) ** 2 for row in rows))
        ),
        "nested_commutator_norm_sum": sum(
            float(row["left_nested_commutator_norm"])
            + float(row["right_nested_commutator_norm"])
            for row in rows
        ),
    }
    return blocks, rows, aggregates


def product_formula_unitary(
    vertices: Sequence[np.ndarray],
    vertex_order: Sequence[int],
    total_time: float,
    nsteps: int,
    order: int,
    exponential_cache: dict[tuple[int, float], np.ndarray] | None = None,
) -> np.ndarray:
    if order not in (1, 2):
        raise ValueError("Only first- and second-order product formulas are supported")
    if nsteps < 1:
        raise ValueError("nsteps must be positive")

    cache = exponential_cache if exponential_cache is not None else {}
    identity = np.eye(vertices[0].shape[0], dtype=complex)
    dt = total_time / nsteps
    if order == 1:
        schedule = [(index, dt) for index in vertex_order]
    else:
        schedule = (
            [(index, dt / 2) for index in vertex_order]
            + [(index, dt / 2) for index in reversed(vertex_order)]
        )

    step = identity.copy()
    for index, duration in schedule:
        cache_key = (index, duration)
        if cache_key not in cache:
            cache[cache_key] = expm(-1j * duration * vertices[index])
        step = cache[cache_key] @ step
    return np.linalg.matrix_power(step, nsteps)


def hartree_fock_state(n_qubits: int, active_electrons: int) -> np.ndarray:
    if not 0 <= active_electrons <= n_qubits:
        raise ValueError("Active electron count must be between zero and n_qubits")
    index = sum(1 << (n_qubits - 1 - qubit) for qubit in range(active_electrons))
    state = np.zeros(2**n_qubits, dtype=complex)
    state[index] = 1.0
    return state


def particle_number_leakage(state: np.ndarray, n_qubits: int, particles: int) -> float:
    sector_probability = sum(
        abs(amplitude) ** 2
        for index, amplitude in enumerate(state)
        if index.bit_count() == particles
    )
    leakage = 1.0 - float(sector_probability.real)
    return float(np.clip(leakage, 0.0, 1.0))


def particle_number_moments(state: np.ndarray) -> tuple[float, float]:
    probabilities = np.abs(state) ** 2
    numbers = np.fromiter(
        (index.bit_count() for index in range(state.size)), dtype=float, count=state.size
    )
    expectation = float(np.dot(probabilities, numbers).real)
    second_moment = float(np.dot(probabilities, numbers**2).real)
    return expectation, max(0.0, second_moment - expectation**2)


def is_scalar_identity(matrix: np.ndarray, tolerance: float = 1e-10) -> bool:
    dimension = matrix.shape[0]
    scalar = np.trace(matrix) / dimension
    return norm(matrix - scalar * np.eye(dimension), 2) <= tolerance


def resolve_cases(repo_root: Path, manifest: Path | None) -> list[CaseSpec]:
    if manifest is None:
        cases = list(DEFAULT_CASES)
    else:
        raw_cases = json.loads(manifest.read_text())
        cases = [CaseSpec(**item) for item in raw_cases]
    missing = [
        repo_root / relative
        for case in cases
        for relative in (case.dat_path, case.tensor_path, case.initial_state_path)
        if not (repo_root / relative).exists()
    ]
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(
            "Missing generated L-sweep benchmark inputs:\n"
            f"{formatted}\n"
            "Generate the molecular library first; see "
            "QHAT_L_SWEEP_ANALYSIS_REPRODUCTION.md."
        )
    return cases


def json_list(values: Iterable[object]) -> str:
    return json.dumps(list(values), separators=(",", ":"))


def run_case(
    case: CaseSpec,
    repo_root: Path,
    total_time: float,
    steps: Sequence[int],
    random_orderings: int,
    root_seed: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    dat_path = repo_root / case.dat_path
    tensor_path = repo_root / case.tensor_path
    initial_state_path = repo_root / case.initial_state_path
    metadata, pauli_terms = parse_qhat_dat(dat_path)
    fermion_vertices, fermion_labels, n_qubits = hermitian_vertices_from_tensors(
        tensor_path
    )
    pauli_words = [word for word, _ in pauli_terms]
    pauli_vertices = [coefficient * pauli_matrix(word) for word, coefficient in pauli_terms]

    if n_qubits != case.active_occupied + case.active_vacant:
        raise ValueError(f"Active-space size mismatch for {case.case_id}")
    metadata_qubits = int(metadata.get("number of qubits", n_qubits))
    if metadata_qubits != n_qubits:
        raise ValueError(f"DAT/tensor qubit mismatch for {case.case_id}")
    metadata_bond_length = float(
        metadata.get("interatomic separation (angstroms)", case.bond_length)
    )
    if not np.isclose(metadata_bond_length, case.bond_length, atol=1e-12):
        raise ValueError(f"DAT/manifest bond-length mismatch for {case.case_id}")

    zero = np.zeros_like(fermion_vertices[0])
    hamiltonian_tensor = sum(fermion_vertices, zero.copy())
    hamiltonian_dat = sum(pauli_vertices, zero.copy())
    reconstruction_error = float(norm(hamiltonian_tensor - hamiltonian_dat, 2))
    if reconstruction_error > 1e-8:
        raise ValueError(
            f"Tensor-to-JW reconstruction failed for {case.case_id}: "
            f"{reconstruction_error:.3e}"
        )

    exact_unitary = expm(-1j * total_time * hamiltonian_dat)
    initial_state = np.load(initial_state_path, allow_pickle=False)
    if initial_state.shape != (2**n_qubits,):
        raise ValueError(f"Initial-state shape mismatch for {case.case_id}")
    if not np.isclose(norm(initial_state), 1.0):
        raise ValueError(f"Initial state is not normalized for {case.case_id}")
    exact_state = exact_unitary @ initial_state
    exact_leakage = particle_number_leakage(
        exact_state, n_qubits, case.active_occupied
    )
    if exact_leakage > 1e-10:
        raise ValueError(
            f"Exact evolution leaks particle number for {case.case_id}: {exact_leakage:.3e}"
        )
    exact_number_expectation, exact_number_variance = particle_number_moments(exact_state)

    decompositions = (
        (
            "fermionic_commutation_then_JW",
            fermion_vertices,
            fermion_labels,
            matrix_noncommutation_graph(fermion_vertices),
        ),
        (
            "JW_Pauli_string_commutation",
            pauli_vertices,
            pauli_words,
            pauli_noncommutation_graph(pauli_words),
        ),
    )

    result_rows: list[dict[str, object]] = []
    graph_rows: list[dict[str, object]] = []
    commutator_rows: list[dict[str, object]] = []

    for scheme, vertices, labels, graph in decompositions:
        coloring, groups = greedy_groups(graph)
        default_order = colored_vertex_order(groups)
        blocks, block_rows, block_aggregates = block_commutators(vertices, groups)
        block_norms = [float(norm(block, 2)) for block in blocks]
        graph_rows.append(
            {
                **asdict(case),
                "scheme": scheme,
                "qubits": n_qubits,
                "vertices": len(vertices),
                "noncommuting_edges": graph.number_of_edges(),
                "colors": len(groups),
                "coloring_json": json.dumps(coloring, sort_keys=True),
                "color_groups_json": json.dumps(groups, sort_keys=True),
                "color_order_json": json_list(groups.keys()),
                "colored_vertex_order_json": json_list(default_order),
                "vertex_labels_json": json.dumps(labels),
                "block_spectral_norms_json": json_list(block_norms),
                "block_commutator_norm_sum": sum(
                    float(row["commutator_spectral_norm"]) for row in block_rows
                ),
                "block_commutator_norm_max": max(
                    (float(row["commutator_spectral_norm"]) for row in block_rows),
                    default=0.0,
                ),
                **block_aggregates,
                "tensor_to_jw_reconstruction_error": reconstruction_error,
                "exact_particle_number_leakage": exact_leakage,
                "exact_particle_number_expectation": exact_number_expectation,
                "exact_particle_number_variance": exact_number_variance,
            }
        )
        for row in block_rows:
            commutator_rows.append({**asdict(case), "scheme": scheme, **row})

        identity_indices = [
            index for index, vertex in enumerate(vertices) if is_scalar_identity(vertex)
        ]
        randomizable_indices = [
            index for index in range(len(vertices)) if index not in identity_indices
        ]
        scheme_seed = np.random.SeedSequence(
            [
                root_seed,
                zlib.crc32(case.case_id.encode()),
                zlib.crc32(scheme.encode()),
            ]
        )
        rng = np.random.default_rng(scheme_seed)
        # Sampling with replacement avoids an infinite uniqueness loop for a custom
        # manifest containing very few vertices. Duplicates are vanishingly unlikely
        # for the default cases and are valid independent random trials if they occur.
        random_permutations = [
            identity_indices + rng.permutation(randomizable_indices).tolist()
            for _ in range(random_orderings)
        ]
        orderings = [("colored", 0, default_order)] + [
            ("random", sample + 1, permutation)
            for sample, permutation in enumerate(random_permutations)
        ]
        for formula_order in (1, 2):
            for nsteps in steps:
                exponential_cache: dict[tuple[int, float], np.ndarray] = {}
                for ordering_kind, ordering_id, vertex_order in orderings:
                    approximate = product_formula_unitary(
                        vertices,
                        vertex_order,
                        total_time,
                        nsteps,
                        formula_order,
                        exponential_cache,
                    )
                    approximate_state = approximate @ initial_state
                    particle_expectation, particle_variance = particle_number_moments(
                        approximate_state
                    )
                    overlap = np.vdot(exact_state, approximate_state)
                    result_rows.append(
                        {
                            **asdict(case),
                            "scheme": scheme,
                            "qubits": n_qubits,
                            "formula_order": formula_order,
                            "steps": nsteps,
                            "total_time": total_time,
                            "ordering_kind": ordering_kind,
                            "ordering_id": ordering_id,
                            "vertex_order_json": json_list(vertex_order),
                            "spectral_error": float(norm(approximate - exact_unitary, 2)),
                            "state_error": float(norm(approximate_state - exact_state)),
                            "state_infidelity": max(0.0, 1.0 - float(abs(overlap) ** 2)),
                            "particle_number_leakage": particle_number_leakage(
                                approximate_state,
                                n_qubits,
                                case.active_occupied,
                            ),
                            "particle_number_expectation": particle_expectation,
                            "particle_number_variance": particle_variance,
                        }
                    )
    return result_rows, graph_rows, commutator_rows


def parse_steps(value: str) -> tuple[int, ...]:
    steps = tuple(int(token.strip()) for token in value.split(",") if token.strip())
    if not steps or any(step < 1 for step in steps):
        raise argparse.ArgumentTypeError("Steps must be a comma-separated list of positives")
    if len(set(steps)) != len(steps):
        raise argparse.ArgumentTypeError("Steps must not contain duplicates")
    return steps


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--manifest", type=Path, help="Optional JSON list of CaseSpec objects")
    parser.add_argument("--output-dir", type=Path, default=Path("trotter_benchmark_results"))
    parser.add_argument("--random-orderings", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260721)
    parser.add_argument("--time", type=float, default=1.0, dest="total_time")
    parser.add_argument("--steps", type=parse_steps, default=(1, 2, 5, 10))
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.random_orderings < 1:
        raise ValueError("--random-orderings must be positive")
    if args.total_time <= 0:
        raise ValueError("--time must be positive")

    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir
    if not output_dir.is_absolute():
        output_dir = repo_root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = resolve_cases(repo_root, args.manifest)
    all_results: list[dict[str, object]] = []
    all_graphs: list[dict[str, object]] = []
    all_commutators: list[dict[str, object]] = []
    for index, case in enumerate(cases, start=1):
        print(f"[{index}/{len(cases)}] {case.case_id}: {case.selection_reason}", flush=True)
        results, graphs, commutators = run_case(
            case,
            repo_root,
            args.total_time,
            args.steps,
            args.random_orderings,
            args.seed,
        )
        all_results.extend(results)
        all_graphs.extend(graphs)
        all_commutators.extend(commutators)
        # Checkpoint each completed case so an interrupted long run retains work.
        pd.DataFrame(all_results).to_csv(output_dir / "ordering_trials.csv", index=False)
        pd.DataFrame(all_graphs).to_csv(output_dir / "graph_diagnostics.csv", index=False)
        pd.DataFrame(all_commutators).to_csv(
            output_dir / "block_commutators.csv", index=False
        )
    (output_dir / "benchmark_config.json").write_text(
        json.dumps(
            {
                "cases": [asdict(case) for case in cases],
                "random_orderings": args.random_orderings,
                "seed": args.seed,
                "total_time": args.total_time,
                "steps": list(args.steps),
                "implemented_orders": [1, 2],
                "randomization_unit": "individual nonidentity decomposition vertices",
            },
            indent=2,
        )
        + "\n"
    )
    print(f"Wrote benchmark outputs to {output_dir}")


if __name__ == "__main__":
    main()
