"""Cirq interop of the three evolution bloqs: base, adjoint, powers and controlled versions.

Three routes to Cirq are each compared with an explicitly built matrix (the controlled unitary is
U on the matching control block and the identity elsewhere), exactly, not up to a phase:
`cirq.unitary(BloqAsCirqGate(bloq))`, `bloq.as_cirq_op(...)` in a circuit, and the decomposition
converted with `to_cirq_circuit`. The terms include Y and terms that don't commute, because
Qualtran 0.4 converts `SGate(is_adjoint=True)` to a plain `cirq.S` (see `PauliStringEvolution`).
"""

import cirq
import numpy as np
import pytest
from qualtran import CtrlSpec
from qualtran._infra.gate_with_registers import get_named_qubits
from qualtran.cirq_interop import BloqAsCirqGate
from scipy.linalg import expm

from qhat.common.commuting_pauli_string_evolution import CommutingPauliStringEvolution
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.pauli_utils import pauli_string_to_matrix
from qhat.common.trotter_flattened import Trotterization

N = 3
TIME = 0.43
TERMS = (("XYZ", 0.7), ("IZY", -0.45), ("YXX", 0.31))  # Y in every term; they don't commute
COMMUTING_TERMS = (("XZI", 0.7), ("IZY", -0.3), ("XZY", 0.2))


def _matrix(pauli_string, coefficient, time=TIME):
    return expm(-1j * coefficient * time * pauli_string_to_matrix(pauli_string))


def _trotter_matrix(trotter):
    """Product of the term exponentials in the expanded sequence (later terms on the left)."""
    u = np.eye(2 ** N, dtype=complex)
    for idx, coeff in trotter.expanded_sequence:
        pauli_string, h = trotter.pauli_terms[idx]
        u = _matrix(pauli_string, h * coeff * trotter.time / trotter.num_steps, 1.0) @ u
    return u


def _controlled(u, ctrl_values):
    """U on the block where the control qubits (first, in order) match `ctrl_values`."""
    k = len(ctrl_values)
    dim = u.shape[0]
    out = np.eye(2 ** k * dim, dtype=complex)
    block = int("".join(str(c) for c in ctrl_values), 2)
    out[block * dim:(block + 1) * dim, block * dim:(block + 1) * dim] = u
    return out


def _controlled_bloq(bloq, ctrl_values):
    for value in reversed(ctrl_values):  # each `.controlled` prepends one control
        bloq = bloq.controlled(CtrlSpec(cvs=value))
    return bloq


def _cases():
    """(id, bloq factory, expected matrix) for U, U^dagger, U**3 and U**-2."""
    pse = PauliStringEvolution("XYZ", 0.7, TIME)
    cpse = CommutingPauliStringEvolution(COMMUTING_TERMS, TIME)
    trotter = Trotterization.from_method(TERMS, "second order", TIME, 3)
    u_pse = _matrix("XYZ", 0.7)
    u_cpse = expm(-1j * TIME * sum(h * pauli_string_to_matrix(p) for p, h in COMMUTING_TERMS))
    u_tr = _trotter_matrix(trotter)
    for name, bloq, u in (("pse", pse, u_pse), ("cpse", cpse, u_cpse), ("trotter", trotter, u_tr)):
        yield f"{name}-U", bloq, u
        yield f"{name}-adjoint", bloq.adjoint(), u.conj().T
        yield f"{name}-pow3", bloq ** 3, np.linalg.matrix_power(u, 3)
        yield f"{name}-pow-2", bloq ** -2, np.linalg.matrix_power(u.conj().T, 2)


CASES = list(_cases())
IDS = [c[0] for c in CASES]
CTRL_VALUES = [(), (1,), (0,), (1, 0)]


def _routes(bloq):
    """The three Cirq unitaries of `bloq`, with qubits ordered as the signature lists them."""
    qubits = get_named_qubits(bloq.signature.lefts())
    order = [q for reg in bloq.signature for q in qubits[reg.name].ravel()]
    ops, _ = bloq.as_cirq_op(cirq.SimpleQubitManager(), **qubits)
    return {
        "BloqAsCirqGate": cirq.unitary(BloqAsCirqGate(bloq)),
        "as_cirq_op": cirq.Circuit(ops).unitary(qubit_order=order),
        "decomposition": bloq.decompose_bloq().to_cirq_circuit(cirq_quregs=qubits).unitary(
            qubit_order=order),
    }


@pytest.mark.parametrize("ctrl_values", CTRL_VALUES, ids=str)
@pytest.mark.parametrize("name, bloq, u", CASES, ids=IDS)
def test_cirq_unitary_is_exact(name, bloq, u, ctrl_values):
    if len(ctrl_values) == 2:
        bloq = _with_direct(bloq)
    expected = _controlled(u, ctrl_values) if ctrl_values else u
    for route, actual in _routes(_controlled_bloq(bloq, ctrl_values)).items():
        assert actual.shape == expected.shape, route
        np.testing.assert_allclose(actual, expected, atol=1e-10, err_msg=route)


def _with_direct(bloq):
    """Two controls with the "direct" method, the ancilla-free one (see the test below)."""
    import attrs
    return attrs.evolve(bloq, multi_control_method="direct")


@pytest.mark.parametrize("name, bloq, u", CASES[:4], ids=IDS[:4])
def test_controlled_adjoint_is_exact(name, bloq, u):
    controlled = bloq.controlled().adjoint()
    expected = _controlled(u, (1,)).conj().T
    for route, actual in _routes(controlled).items():
        np.testing.assert_allclose(actual, expected, atol=1e-10, err_msg=route)


def test_wrong_without_the_s_dagger_fix_is_detected():
    """The Y basis change must be right in Cirq: `Y = S H Z H S^dagger`."""
    for bloq in (PauliStringEvolution("Y", 0.7, TIME), PauliStringEvolution("IYX", 0.7, TIME)):
        expected = _matrix(bloq.pauli_string, 0.7)
        np.testing.assert_allclose(cirq.unitary(BloqAsCirqGate(bloq)), expected, atol=1e-10)


@pytest.mark.parametrize("is_adjoint, matrix", [(False, np.diag([1, 1j])), (True, np.diag([1, -1j]))])
def test_cirq_safe_s_gate(is_adjoint, matrix):
    """`CirqSafeSGate` is S or S^dagger in Qualtran and in Cirq, and 1 Clifford in every counter
    (plain `SGate(is_adjoint=True)` becomes `cirq.S` in Cirq)."""
    from pyLIQTR.utils.resource_analysis import estimate_resources
    from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity
    from qualtran.resource_counting import GateCounts, QECGatesCost, get_cost_value
    from qhat.common.pauli_string_evolution import CirqSafeSGate
    gate = CirqSafeSGate(is_adjoint=is_adjoint)
    np.testing.assert_allclose(gate.tensor_contract(), matrix, atol=1e-12)
    np.testing.assert_allclose(cirq.unitary(BloqAsCirqGate(gate)), matrix, atol=1e-12)
    q = cirq.LineQubit(0)
    ops, _ = gate.as_cirq_op(cirq.SimpleQubitManager(), q=np.array([q]))
    np.testing.assert_allclose(cirq.unitary(cirq.Circuit(ops)), matrix, atol=1e-12)
    assert gate.adjoint() == CirqSafeSGate(is_adjoint=not is_adjoint)
    assert t_complexity(gate) == TComplexity(clifford=1)
    assert get_cost_value(gate, QECGatesCost()) == GateCounts(clifford=1)
    resources = estimate_resources(gate, profile=True)
    assert (resources["T"], resources["Clifford"], resources["Rotations"]) == (0, 1, 0)


def test_y_basis_change_uses_one_s_dagger():
    """A Y costs S^dagger, H, H, S (4 Cliffords), with no plain `SGate(is_adjoint=True)`."""
    from collections import Counter
    from qualtran.bloqs.basic_gates import Hadamard, SGate
    from qhat.common.pauli_string_evolution import CirqSafeSGate
    bloq = PauliStringEvolution("IYX", 0.7, TIME)
    gates = Counter(bi.bloq for bi in bloq.decompose_bloq().bloq_instances)
    assert gates[CirqSafeSGate(is_adjoint=True)] == 1
    assert gates[CirqSafeSGate()] == 1
    assert gates[Hadamard()] == 4
    assert not any(type(g) is SGate for g in gates)
    assert bloq.t_complexity().clifford == 4 + 2 + 2    # Y, X, 2 CNOTs


def test_ancilla_method_has_no_cirq_unitary():
    """Known limitation: the AND ladder allocates an ancilla, which `cirq.unitary` can't contract.

    The circuit itself is still built (for resource estimation and diagrams)."""
    bloq = PauliStringEvolution("XYZ", 0.7, TIME).controlled().controlled()
    assert bloq.multi_control_method == "ancilla"
    qubits = get_named_qubits(bloq.signature.lefts())
    ops, _ = bloq.as_cirq_op(cirq.SimpleQubitManager(), **qubits)
    assert len(cirq.Circuit(ops).all_qubits()) == 5
    with pytest.raises(TypeError):
        cirq.unitary(BloqAsCirqGate(bloq))


def test_bloq_as_cirq_gate_has_the_right_qubits():
    trotter = Trotterization.from_method(TERMS, "second order", TIME, 3).controlled()
    gate = BloqAsCirqGate(trotter)
    assert cirq.num_qubits(gate) == N + 1
    assert gate.on(*cirq.LineQubit.range(N + 1)).qubits == tuple(cirq.LineQubit.range(N + 1))
