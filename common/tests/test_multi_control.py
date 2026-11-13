"""Two or more controls on the evolution bloqs, by both multi-control methods (plan D1).

Every matrix check is exact (no phase freedom), against a block matrix built here: U on the block
matching the control values (controls most significant), identity elsewhere.
"""

import attrs
import numpy as np
import pytest
from scipy.linalg import expm

from qualtran import BloqBuilder, CtrlSpec
from qualtran.bloqs.mcmt import And, MultiAnd
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity
from qualtran.resource_counting import QubitCount, get_cost_value
from qualtran.resource_counting._qubit_counts import _cbloq_max_width

from qhat.common import controls
from qhat.common.commuting_pauli_string_evolution import (CommutingPauliStringEvolution,
                                                          all_commute)
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.pauli_utils import pauli_string_to_matrix
from qhat.common.trotter_flattened import Trotterization, get_trotterization_coefficients

METHODS = ["ancilla", "direct"]
# 2 and 3 controls, mixed values
CTRL_VALUES = [(1, 0), (0, 1), (1, 1), (0, 0), (1, 0, 1), (0, 1, 1), (0, 0, 0)]


def block_matrix(unitary, ctrl_values):
    dim = unitary.shape[0]
    block = int(''.join(str(v) for v in ctrl_values), 2)
    result = np.eye(dim * 2 ** len(ctrl_values), dtype=complex)
    result[block * dim:(block + 1) * dim, block * dim:(block + 1) * dim] = unitary
    return result


def decomposition_width(bloq):
    """`QubitCount` worked out from the decomposition, bypassing the bloq's `my_static_costs`."""
    return _cbloq_max_width(bloq.decompose_bloq()._binst_graph,
                            lambda b: get_cost_value(b, QubitCount()))


def check_exact(bloq, expected):
    np.testing.assert_allclose(bloq.tensor_contract(), expected, atol=1e-11)
    np.testing.assert_allclose(bloq.decompose_bloq().tensor_contract(), expected, atol=1e-11)


def and_cost(ctrl_values):
    """The AND ladder and its uncompute (Qualtran 0.4: `And` 4 T + 9 Cliff, uncompute 0 T)."""
    ladder = And(*ctrl_values) if len(ctrl_values) == 2 else MultiAnd(cvs=ctrl_values)
    return t_complexity(ladder) + t_complexity(ladder.adjoint())


# -------------------------------------------------------------------------------------------------
# Helpers in common/controls.py
# -------------------------------------------------------------------------------------------------

class TestControlledRotationHelpers:
    @pytest.mark.parametrize("ctrl_values", [(1,), (0,)] + CTRL_VALUES + [(1, 0, 0, 1)])
    def test_controlled_rz_is_exact(self, ctrl_values):
        theta = 0.737
        k = len(ctrl_values)
        bb = BloqBuilder()
        ctrl = [bb.add_register(f'c{i}', 1) for i in range(k)]
        target = bb.add_register('t', 1)
        target, ctrl = controls.add_controlled_rz(bb, theta, target, ctrl, ctrl_values)
        cbloq = bb.finalize(**{f'c{i}': ctrl[i] for i in range(k)}, t=target)
        rz = np.diag([np.exp(-0.5j * theta), np.exp(0.5j * theta)])
        np.testing.assert_allclose(cbloq.tensor_contract(), block_matrix(rz, ctrl_values),
                                   atol=1e-12)
        assert t_complexity(cbloq) == TComplexity(clifford=2 ** k, rotations=2 ** k)
        assert t_complexity(cbloq) == controls.controlled_rz_t_complexity(k)

    @pytest.mark.parametrize("ctrl_values", [(1,), (0,)] + CTRL_VALUES + [(1, 0, 0, 1)])
    def test_controlled_phase_is_exact(self, ctrl_values):
        phi = 0.41
        k = len(ctrl_values)
        bb = BloqBuilder()
        ctrl = [bb.add_register(f'c{i}', 1) for i in range(k)]
        ctrl = controls.add_controlled_phase(bb, phi, ctrl, ctrl_values)
        cbloq = bb.finalize(**{f'c{i}': ctrl[i] for i in range(k)})
        expected = block_matrix(np.array([[np.exp(-1j * phi)]]), ctrl_values)
        np.testing.assert_allclose(cbloq.tensor_contract(), expected, atol=1e-12)
        assert t_complexity(cbloq) == TComplexity(clifford=2 ** k - 2, rotations=2 ** k - 1)
        assert t_complexity(cbloq) == controls.controlled_phase_t_complexity(k)

    def test_single_control_is_the_textbook_circuit(self):
        """One control on |1>: Rz(theta/2) . CNOT . Rz(-theta/2) . CNOT (plan D4)."""
        from qualtran.bloqs.basic_gates import CNOT, Rz
        bb = BloqBuilder()
        ctrl, target = bb.add_register('c', 1), bb.add_register('t', 1)
        target, (ctrl,) = controls.add_controlled_rz(bb, 0.6, target, [ctrl], (1,))
        cbloq = bb.finalize(c=ctrl, t=target)
        import networkx as nx
        from qualtran import BloqInstance
        order = [node.bloq for node in nx.topological_sort(cbloq._binst_graph)
                 if isinstance(node, BloqInstance)]
        assert order == [Rz(0.3), CNOT(), Rz(-0.3), CNOT()]

    def test_peak_qubits(self):
        assert controls.peak_qubits(10, (1, 0, 1), "ancilla") == 15
        assert controls.peak_qubits(10, (1, 0, 1), "direct") == 13
        assert controls.peak_qubits(10, (0,), "ancilla") == 11
        assert controls.peak_qubits(10, (), "ancilla") == 10


# -------------------------------------------------------------------------------------------------
# PauliStringEvolution
# -------------------------------------------------------------------------------------------------

PAULI_STRINGS = ["Z", "Y", "XZ", "IYI", "XIZ", "II"]


class TestPauliStringEvolution:
    @pytest.mark.parametrize("pauli_string", PAULI_STRINGS)
    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    @pytest.mark.parametrize("method", METHODS)
    @pytest.mark.parametrize("is_adjoint", [False, True])
    def test_matches_block_matrix(self, pauli_string, ctrl_values, method, is_adjoint):
        pse = PauliStringEvolution(pauli_string, -0.62, 1.37, 1.7, is_adjoint=is_adjoint,
                                   ctrl_values=ctrl_values, multi_control_method=method)
        forward = expm(-1j * -0.62 * 1.37 / 1.7 * pauli_string_to_matrix(pauli_string))
        check_exact(pse, block_matrix(forward.conj().T if is_adjoint else forward, ctrl_values))

    @pytest.mark.parametrize("pauli_string", PAULI_STRINGS)
    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    @pytest.mark.parametrize("method", METHODS)
    def test_count_matches_decomposition(self, pauli_string, ctrl_values, method):
        pse = PauliStringEvolution(pauli_string, 0.37, 0.91, ctrl_values=ctrl_values,
                                   multi_control_method=method)
        assert t_complexity(pse) == t_complexity(pse.decompose_bloq())

    @pytest.mark.parametrize("ctrl_values", [(1, 0), (0, 1, 1)])
    def test_hand_count(self, ctrl_values):
        """`XYZI`: basis changes 2 (X) + 4 (Y) and a 4-CNOT ladder are 10 Cliffords."""
        k = len(ctrl_values)
        pse = PauliStringEvolution("XYZI", 0.37, 0.91, ctrl_values=ctrl_values)
        direct = attrs.evolve(pse, multi_control_method="direct")
        # "direct": the rotation becomes 2^k rotations and 2^k CNOTs
        assert t_complexity(direct) == TComplexity(clifford=10 + 2 ** k, rotations=2 ** k)
        # "ancilla": the AND ladder, then the singly-controlled term (2 rotations, 2 CNOTs)
        assert t_complexity(pse) == and_cost(ctrl_values) + TComplexity(clifford=12, rotations=2)
        assert and_cost(ctrl_values).t == 4 * (k - 1)

    @pytest.mark.parametrize("ctrl_values", [(1, 0), (0, 1, 1)])
    def test_hand_count_all_identity(self, ctrl_values):
        k = len(ctrl_values)
        pse = PauliStringEvolution("II", 0.37, 0.91, ctrl_values=ctrl_values)
        direct = attrs.evolve(pse, multi_control_method="direct")
        assert t_complexity(direct) == TComplexity(clifford=2 ** k - 2, rotations=2 ** k - 1)
        assert t_complexity(pse) == and_cost(ctrl_values) + TComplexity(rotations=1)

    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    @pytest.mark.parametrize("method", METHODS)
    def test_qubit_count(self, ctrl_values, method):
        k = len(ctrl_values)
        pse = PauliStringEvolution("XIZ", 0.37, 0.91, ctrl_values=ctrl_values,
                                   multi_control_method=method)
        expected = 3 + k + (k - 1 if method == "ancilla" else 0)
        assert get_cost_value(pse, QubitCount()) == expected
        assert decomposition_width(pse) == expected

    def test_direct_uses_no_and(self):
        pse = PauliStringEvolution("XZ", 0.37, 0.91, ctrl_values=(1, 0, 1),
                                   multi_control_method="direct")
        names = {type(binst.bloq).__name__ for binst in pse.decompose_bloq().bloq_instances}
        assert names <= {"Rz", "CNOT", "Hadamard"}

    def test_ancilla_ands_once_around_one_control(self):
        pse = PauliStringEvolution("XZ", 0.37, 0.91, ctrl_values=(1, 0, 1))
        bloqs = [binst.bloq for binst in pse.decompose_bloq().bloq_instances]
        assert len(bloqs) == 3
        assert attrs.evolve(pse, ctrl_values=(1,)) in bloqs
        assert MultiAnd(cvs=(1, 0, 1)) in bloqs


# -------------------------------------------------------------------------------------------------
# CommutingPauliStringEvolution
# -------------------------------------------------------------------------------------------------

class TestCommutingPauliStringEvolution:
    TERMS = (("XZI", 0.4), ("ZXI", -0.7), ("IIY", 0.3), ("III", 0.9))  # commuting, with I and III

    def forward(self):
        result = np.eye(8, dtype=complex)
        for pauli_string, coefficient in self.TERMS:
            result = expm(-1j * coefficient * 0.83 * pauli_string_to_matrix(pauli_string)) @ result
        return result

    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    @pytest.mark.parametrize("method", METHODS)
    @pytest.mark.parametrize("is_adjoint", [False, True])
    def test_matches_block_matrix(self, ctrl_values, method, is_adjoint):
        cpse = CommutingPauliStringEvolution(self.TERMS, 0.83, is_adjoint=is_adjoint,
                                             ctrl_values=ctrl_values, multi_control_method=method)
        forward = self.forward()
        check_exact(cpse, block_matrix(forward.conj().T if is_adjoint else forward, ctrl_values))

    @pytest.mark.parametrize("ctrl_values", [(1, 0), (0, 1, 1)])
    @pytest.mark.parametrize("method", METHODS)
    def test_counts(self, ctrl_values, method):
        k = len(ctrl_values)
        cpse = CommutingPauliStringEvolution(self.TERMS, 0.83, ctrl_values=ctrl_values,
                                             multi_control_method=method)
        assert t_complexity(cpse) == t_complexity(cpse.decompose_bloq())
        one = CommutingPauliStringEvolution(self.TERMS, 0.83, ctrl_values=(1,))
        if method == "ancilla":
            # One AND ladder for the whole bloq, not one per term
            assert t_complexity(cpse) == and_cost(ctrl_values) + t_complexity(one)
        expected = 3 + k + (k - 1 if method == "ancilla" else 0)
        assert get_cost_value(cpse, QubitCount()) == expected
        assert decomposition_width(cpse) == expected

    def test_direct_terms_carry_all_controls(self):
        cpse = CommutingPauliStringEvolution(self.TERMS, 0.83, ctrl_values=(1, 0),
                                             multi_control_method="direct")
        for binst in cpse.decompose_bloq().bloq_instances:
            assert binst.bloq.ctrl_values == (1, 0)
            assert binst.bloq.multi_control_method == "direct"


# -------------------------------------------------------------------------------------------------
# Trotterization
# -------------------------------------------------------------------------------------------------

def reference_trotter_matrix(pauli_terms, coefficients, time, num_steps, hbar):
    """Trotter product from expm of each term; ramps alternate forward and backward term order."""
    dt = time / num_steps
    result = np.eye(2 ** len(pauli_terms[0][0]), dtype=complex)
    for _ in range(num_steps):
        for i, coeff in enumerate(coefficients):
            for pauli_string, h in (pauli_terms if i % 2 == 0 else pauli_terms[::-1]):
                result = expm(-1j * h * coeff * dt * pauli_string_to_matrix(pauli_string)
                              / hbar) @ result
    return result


class TestTrotterization:
    PAULI_TERMS = (("XZI", 0.5), ("ZIY", 0.3), ("YXX", -0.2), ("III", 0.4))
    ORDERS = {"first order": (1.0,), "third order": get_trotterization_coefficients(3),
              "custom": (0.3, 0.7)}

    def make(self, order="third order", **kwargs):
        return Trotterization(self.PAULI_TERMS, self.ORDERS[order], 0.43, 3, hbar=0.9, **kwargs)

    def test_terms_do_not_commute(self):
        assert not all_commute([p for p, _ in self.PAULI_TERMS])

    @pytest.mark.parametrize("order", ["first order", "third order", "custom"])
    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    @pytest.mark.parametrize("method", METHODS)
    @pytest.mark.parametrize("is_adjoint", [False, True])
    def test_matches_block_matrix(self, order, ctrl_values, method, is_adjoint):
        trotter = self.make(order, is_adjoint=is_adjoint, ctrl_values=ctrl_values,
                            multi_control_method=method)
        u = reference_trotter_matrix(self.PAULI_TERMS, self.ORDERS[order], 0.43, 3, 0.9)
        check_exact(trotter, block_matrix(u.conj().T if is_adjoint else u, ctrl_values))

    @pytest.mark.parametrize("contraction", ["incremental", "structured", "qualtran"])
    @pytest.mark.parametrize("method", METHODS)
    def test_every_contraction_method(self, contraction, method):
        trotter = self.make(ctrl_values=(0, 1, 0), multi_control_method=method,
                            tensor_contraction_method=contraction)
        u = reference_trotter_matrix(self.PAULI_TERMS, self.ORDERS["third order"], 0.43, 3, 0.9)
        np.testing.assert_allclose(trotter.tensor_contract(), block_matrix(u, (0, 1, 0)),
                                   atol=1e-11)

    @pytest.mark.parametrize("ctrl_values", [(1, 0), (0, 1, 1)])
    @pytest.mark.parametrize("method", METHODS)
    def test_count_matches_decomposition(self, ctrl_values, method):
        trotter = self.make(ctrl_values=ctrl_values, multi_control_method=method)
        assert t_complexity(trotter) == t_complexity(trotter.decompose_bloq())

    @pytest.mark.parametrize("ctrl_values", [(1, 0), (0, 1, 1)])
    def test_hand_count(self, ctrl_values):
        """Per term occurrence: "direct" replaces each rotation by 2^k rotations and 2^k CNOTs;
        "ancilla" is one AND ladder plus the singly-controlled Trotterization."""
        from collections import Counter
        k = len(ctrl_values)
        trotter = self.make(ctrl_values=ctrl_values)
        direct = attrs.evolve(trotter, multi_control_method="direct")
        expected = TComplexity()
        for idx, n in Counter(i for i, _ in trotter.expanded_sequence).items():
            active = [p for p in self.PAULI_TERMS[idx][0] if p != 'I']
            if not active:
                expected += n * TComplexity(clifford=2 ** k - 2, rotations=2 ** k - 1)
                continue
            basis = sum({'X': 2, 'Y': 4, 'Z': 0}[p] for p in active)
            expected += n * TComplexity(clifford=basis + 2 * (len(active) - 1) + 2 ** k,
                                        rotations=2 ** k)
        assert t_complexity(direct) == expected
        assert t_complexity(trotter) == (and_cost(ctrl_values)
                                         + t_complexity(self.make(ctrl_values=(1,))))

    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    @pytest.mark.parametrize("method", METHODS)
    def test_qubit_count(self, ctrl_values, method):
        k = len(ctrl_values)
        trotter = self.make(ctrl_values=ctrl_values, multi_control_method=method)
        expected = 3 + k + (k - 1 if method == "ancilla" else 0)
        assert get_cost_value(trotter, QubitCount()) == expected
        assert decomposition_width(trotter) == expected

    def test_nested_controlled_equals_ctrl_values(self):
        trotter = self.make()
        assert trotter.controlled().controlled() == self.make(ctrl_values=(1, 1))
        assert (trotter.controlled(CtrlSpec(cvs=0)).controlled()
                == self.make(ctrl_values=(1, 0)))
        assert type(trotter.controlled().controlled().controlled()) is Trotterization

    def test_control_order_matters(self):
        assert self.make(ctrl_values=(1, 0)) != self.make(ctrl_values=(0, 1))
        assert hash(self.make(ctrl_values=(1, 0))) != hash(self.make(ctrl_values=(0, 1)))

    def test_method_is_part_of_equality(self):
        assert self.make(ctrl_values=(1, 0)) != self.make(ctrl_values=(1, 0),
                                                          multi_control_method="direct")

    def test_nested_controlled_matches_block_matrix(self):
        """The adder wires a new control into a bloq that already has a `ctrl` register."""
        inner = self.make().controlled(CtrlSpec(cvs=0))
        bb = BloqBuilder()
        c_new = bb.add_register('c_new', 1)
        c_old = bb.add_register('c_old', 1)
        from qualtran import QBit, Register
        q = bb.add_register(Register('q', QBit(), shape=(3,)))
        _, adder = inner.get_ctrl_system(CtrlSpec())
        (c_new,), (ctrl, q) = adder(bb, [c_new], {'ctrl': np.array([c_old]), 'q': q})
        cbloq = bb.finalize(c_new=c_new, c_old=ctrl[0], q=q)
        u = reference_trotter_matrix(self.PAULI_TERMS, self.ORDERS["third order"], 0.43, 3, 0.9)
        np.testing.assert_allclose(cbloq.tensor_contract(), block_matrix(u, (1, 0)), atol=1e-11)

    @pytest.mark.parametrize("method", METHODS)
    def test_adjoint_and_power_keep_controls(self, method):
        trotter = self.make(ctrl_values=(1, 0), multi_control_method=method)
        assert trotter.adjoint().ctrl_values == (1, 0)
        assert trotter.adjoint().multi_control_method == method
        assert (trotter ** -2) == (trotter ** 2).adjoint()
        assert (trotter ** 2).multi_control_method == method
        assert self.make().controlled().adjoint().controlled() == self.make(
            is_adjoint=True, ctrl_values=(1, 1))

    @pytest.mark.parametrize("method", METHODS)
    def test_terms_carry_controls(self, method):
        trotter = self.make(ctrl_values=(1, 0), multi_control_method=method)
        bloqs = [binst.bloq for binst in trotter.decompose_bloq().bloq_instances]
        terms = [b for b in bloqs if isinstance(b, CommutingPauliStringEvolution)]
        if method == "direct":
            assert len(terms) == trotter.num_operations
            assert all(t.ctrl_values == (1, 0) for t in terms)
        else:
            # The AND, one singly-controlled Trotterization, the uncompute
            assert len(bloqs) == 3
            assert attrs.evolve(trotter, ctrl_values=(1,)) in bloqs
            assert And(1, 0) in bloqs and And(1, 0).adjoint() in bloqs

    def test_large_step_count_is_fast(self):
        import time
        trotter = Trotterization.from_method(self.PAULI_TERMS, "fourth order", 0.43, 5000)
        start = time.perf_counter()
        for method in METHODS:
            controlled = attrs.evolve(trotter, ctrl_values=(1, 0, 1), multi_control_method=method)
            t_complexity(controlled)
            get_cost_value(controlled, QubitCount())
            controlled.tensor_contract()
        assert time.perf_counter() - start < 5

    def test_str_and_repr(self):
        trotter = self.make(ctrl_values=(1, 0), multi_control_method="direct")
        assert str(trotter).startswith("C[1,0](")
        assert "ctrl_values=(1, 0)" in repr(trotter)
        assert "multi_control_method='direct'" in repr(trotter)
        assert "multi_control_method" not in repr(self.make(ctrl_values=(1, 0)))

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError, match=r"\('ancilla', 'direct'\)"):
            self.make(multi_control_method="toffoli")

    @pytest.mark.parametrize("method", METHODS)
    def test_method_ignored_with_one_control(self, method):
        """With 0 or 1 controls, both methods give the same circuit and counts."""
        trotter = self.make(ctrl_values=(0,), multi_control_method=method)
        reference = self.make(ctrl_values=(0,))
        assert t_complexity(trotter) == t_complexity(reference)
        assert get_cost_value(trotter, QubitCount()) == 4
        np.testing.assert_allclose(trotter.decompose_bloq().tensor_contract(),
                                   reference.decompose_bloq().tensor_contract(), atol=1e-12)


# -------------------------------------------------------------------------------------------------
# pyLIQTR (plan D10)
# -------------------------------------------------------------------------------------------------

def test_pyliqtr_counts_follow_t_complexity():
    """pyLIQTR uses our `_t_complexity_` (notes F17): rotations become T, the AND ladder adds 4 T
    per extra control."""
    from qhat.analysis.resource_estimation import resource_estimation_pyliqtr
    trotter = Trotterization((("XZ", 0.5), ("ZY", 0.3)), (0.5, 0.5), 0.43, 3)
    one = resource_estimation_pyliqtr(trotter.controlled())
    for ctrl_values in [(1, 0), (0, 1, 1)]:
        k = len(ctrl_values)
        ancilla = resource_estimation_pyliqtr(attrs.evolve(trotter, ctrl_values=ctrl_values))
        assert ancilla["T_count"] == one["T_count"] + 4 * (k - 1)
        direct = resource_estimation_pyliqtr(
            attrs.evolve(trotter, ctrl_values=ctrl_values, multi_control_method="direct"))
        assert direct["T_count"] > ancilla["T_count"]
