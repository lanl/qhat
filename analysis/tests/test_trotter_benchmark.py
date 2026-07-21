"""Validation tests for the focused fermionic-vs-JW benchmark workflow."""

from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm, norm

from run_trotter_benchmark import (
    DEFAULT_CASES,
    canonicalize_tensor_monomial,
    greedy_groups,
    hermitian_vertices_from_tensors,
    matrix_noncommutation_graph,
    parse_qhat_dat,
    pauli_matrix,
    pauli_noncommutation_graph,
    paulis_commute,
    product_formula_unitary,
)


REPO_ROOT = Path(__file__).resolve().parents[2]

EXPECTED_GRAPHS = {
    "similar_4q": ((13, 12, 2), (15, 16, 2)),
    "fermionic_first_order_6q": ((28, 90, 4), (34, 168, 3)),
    "jw_second_order_6q": ((28, 90, 4), (34, 168, 3)),
}

EXPECTED_STEP_10 = {
    "similar_4q": {
        "fermion": (9.137841166e-4, 3.225636544e-6),
        "pauli": (9.137841166e-4, 3.298905728e-6),
    },
    "fermionic_first_order_6q": {
        "fermion": (2.353589861e-3, 4.151929638e-5),
        "pauli": (2.892363911e-3, 5.138663716e-5),
    },
    "jw_second_order_6q": {
        "fermion": (8.786419899e-4, 1.239783101e-5),
        "pauli": (3.633367543e-4, 1.239658190e-6),
    },
}

BENCHMARK_INPUTS_AVAILABLE = all(
    (REPO_ROOT / relative).exists()
    for case in DEFAULT_CASES
    for relative in (case.dat_path, case.tensor_path, case.initial_state_path)
)


def test_pauli_commutation_parity():
    assert not paulis_commute("XI", "ZI")
    assert paulis_commute("XI", "IZ")
    assert paulis_commute("XX", "ZZ")


def test_tensor_monomial_canonicalization_matches_normal_order_convention():
    canonical, sign = canonicalize_tensor_monomial(
        ((0, 1), (1, 1), (2, 0), (3, 0))
    )
    assert canonical == ((1, 1), (0, 1), (3, 0), (2, 0))
    assert sign == 1

    canonical, sign = canonicalize_tensor_monomial(
        ((1, 1), (0, 1), (2, 0), (3, 0))
    )
    assert canonical == ((1, 1), (0, 1), (3, 0), (2, 0))
    assert sign == -1

    assert canonicalize_tensor_monomial(((0, 1), (0, 1))) == (None, 0)


@pytest.mark.parametrize("case", DEFAULT_CASES, ids=lambda case: case.case_id)
@pytest.mark.skipif(
    not BENCHMARK_INPUTS_AVAILABLE,
    reason="generated L-sweep inputs are intentionally not committed",
)
def test_selected_case_reconstruction_graphs_and_colored_errors(case):
    metadata, terms = parse_qhat_dat(REPO_ROOT / case.dat_path)
    fermion_vertices, _, n_qubits = hermitian_vertices_from_tensors(
        REPO_ROOT / case.tensor_path
    )
    pauli_words = [word for word, _ in terms]
    pauli_vertices = [coefficient * pauli_matrix(word) for word, coefficient in terms]
    assert int(metadata["number of qubits"]) == n_qubits

    zero = np.zeros_like(fermion_vertices[0])
    hamiltonian_f = sum(fermion_vertices, zero.copy())
    hamiltonian_p = sum(pauli_vertices, zero.copy())
    assert norm(hamiltonian_f - hamiltonian_p, 2) < 1e-8

    fermion_graph = matrix_noncommutation_graph(fermion_vertices)
    pauli_graph = pauli_noncommutation_graph(pauli_words)
    _, fermion_groups = greedy_groups(fermion_graph)
    _, pauli_groups = greedy_groups(pauli_graph)
    expected_f, expected_p = EXPECTED_GRAPHS[case.case_id]
    assert (len(fermion_vertices), fermion_graph.number_of_edges(), len(fermion_groups)) == expected_f
    assert (len(pauli_vertices), pauli_graph.number_of_edges(), len(pauli_groups)) == expected_p

    exact = expm(-1j * hamiltonian_p)
    for name, vertices, groups in (
        ("fermion", fermion_vertices, fermion_groups),
        ("pauli", pauli_vertices, pauli_groups),
    ):
        vertex_order = [node for color in groups for node in groups[color]]
        for formula_order, expected in zip(
            (1, 2), EXPECTED_STEP_10[case.case_id][name]
        ):
            approximate = product_formula_unitary(
                vertices, vertex_order, total_time=1.0, nsteps=10, order=formula_order
            )
            assert norm(approximate - exact, 2) == pytest.approx(expected, rel=2e-7)
