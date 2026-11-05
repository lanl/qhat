"""
Tests for ExactStatePreparation.
"""

from collections import Counter

import numpy as np
import pytest

from qualtran.bloqs.basic_gates import CNOT, GlobalPhase, Ry, Rz
from qualtran.cirq_interop.t_complexity_protocol import t_complexity
from qualtran.resource_counting import get_cost_value, QubitCount, SympySymbolAllocator

from qhat.common.exact_state_preparation import ExactStatePreparation


def random_state(n, complex_, seed):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=2**n) + (1j * rng.normal(size=2**n) if complex_ else 0)
    return a / np.linalg.norm(a)


def gate_by_gate(bloq):
    """The circuit's matrix, contracted from its decomposition (not the dense override)."""
    return bloq.decompose_bloq().tensor_contract()


def gate_counts(bloq):
    """Number of Ry, Rz, CNOT and GlobalPhase gates in the decomposition (no Split/Join)."""
    counts = Counter()
    for binst in bloq.decompose_bloq().bloq_instances:
        if isinstance(binst.bloq, (Ry, Rz, CNOT, GlobalPhase)):
            counts[type(binst.bloq).__name__] += 1
    return counts


STATES = {
    "real": random_state(3, False, 1),
    "complex": random_state(3, True, 2),
    # zeros, an all-zero prefix, and negative signs
    "sparse": np.array([0, 0.6, 0, 0, -0.48, 0, 0, 0.64]),
    "basis state": np.eye(8)[5],
    "complex with zeros": np.array([0, 1j, 0, 0, 0, 0, -1, 0]) / np.sqrt(2),
}


@pytest.mark.parametrize("name", STATES)
@pytest.mark.parametrize("contract", [ExactStatePreparation.tensor_contract, gate_by_gate],
                         ids=["dense", "gates"])
def test_prepares_amplitudes_exactly(name, contract):
    a = STATES[name]
    np.testing.assert_allclose(contract(ExactStatePreparation(a))[:, 0], a, atol=1e-14)


@pytest.mark.parametrize("n", [1, 2, 4])
@pytest.mark.parametrize("complex_", [False, True])
def test_dense_matrix_is_the_circuit(n, complex_):
    """The whole unitary, not only its first column, matches the gate-by-gate circuit."""
    bloq = ExactStatePreparation(random_state(n, complex_, n))
    u = bloq.tensor_contract()
    np.testing.assert_allclose(u, gate_by_gate(bloq), atol=1e-13)
    np.testing.assert_allclose(u.conj().T @ u, np.eye(2**n), atol=1e-13)


@pytest.mark.parametrize("n", [1, 2, 3, 4])
@pytest.mark.parametrize("complex_", [False, True])
def test_call_graph_counts_match_decomposition(n, complex_):
    bloq = ExactStatePreparation(random_state(n, complex_, 10 + n))
    from_graph = Counter()
    for callee, count in bloq.build_call_graph(SympySymbolAllocator()):
        from_graph[type(callee).__name__] += count
    assert from_graph == gate_counts(bloq)
    expected = {"Ry": 2**n - 1, "CNOT": 2**n - 2}
    if complex_:
        expected = {"Ry": 2**n - 1, "Rz": 2**n - 1, "CNOT": 2 * (2**n - 2), "GlobalPhase": 1}
    assert +from_graph == Counter(expected)


def test_real_amplitudes_use_no_rz():
    """Signs of real amplitudes come from the last Ry level, not from Rz rotations."""
    counts = gate_counts(ExactStatePreparation(STATES["sparse"]))
    assert counts["Rz"] == 0 and counts["GlobalPhase"] == 0


def test_rotation_count_does_not_depend_on_angles():
    """Zero angles still count as rotations, so the cost is that of a generic state."""
    basis = t_complexity(ExactStatePreparation(np.eye(8)[0]))
    generic = t_complexity(ExactStatePreparation(random_state(3, False, 5)))
    assert basis == generic
    assert (generic.rotations, generic.clifford) == (7, 6)


def test_qubit_count_has_no_scratch():
    assert get_cost_value(ExactStatePreparation(random_state(3, True, 7)), QubitCount()) == 3


@pytest.mark.parametrize("amplitudes, match", [
    ([1.0], "2\\^n"), ([0.6, 0.0, 0.8], "2\\^n"), ([1.0, 1.0], "unit norm"),
])
def test_invalid_amplitudes_raise(amplitudes, match):
    with pytest.raises(ValueError, match=match):
        ExactStatePreparation(amplitudes)
