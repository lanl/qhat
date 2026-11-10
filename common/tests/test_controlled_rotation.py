"""The "ancilla" controlled-rotation method, option (c) of plan T10.

Each controlled rotation ANDs the control and the parity qubit into an ancilla, rotates it and
un-ANDs, leaving a phase on the control. A parent bloq applies the sum of its terms' phases once.
Every matrix check is exact (no phase freedom), against a block matrix built here, so a missing
or doubled phase fails.
"""

import time

import attrs
import numpy as np
import pytest
from scipy.linalg import expm

from qualtran import BloqBuilder, CtrlSpec
from qualtran.bloqs.basic_gates import GlobalPhase, Rz
from qualtran.bloqs.mcmt import And
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity
from qualtran.resource_counting import QubitCount, get_cost_value
from qualtran.resource_counting._qubit_counts import _cbloq_max_width

from qhat.common import controls
from qhat.common.commuting_pauli_string_evolution import CommutingPauliStringEvolution
from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.pauli_utils import pauli_string_to_matrix
from qhat.common.tests.test_multi_control import (block_matrix, check_exact, decomposition_width,
                                                  reference_trotter_matrix)
from qhat.common.trotter_flattened import Trotterization

ANCILLA = "ancilla"
CTRL_VALUES = [(1,), (0,), (1, 0), (0, 1, 1)]
PAULI_STRINGS = ["Z", "Y", "XZ", "IYI", "XYZI", "II"]
# Terms that don't commute, with an all-identity term; odd time and hbar
TROTTER_TERMS = (("XZI", 0.4), ("ZIY", -0.3), ("YXX", 0.25), ("III", 0.6))
TIME, HBAR = 0.43, 0.9


def and_pair(ctrl_value):
    """`And` and its uncompute: 4 T, and 13 Cliffords (control on |1>) or 17 (on |0>)."""
    return t_complexity(And(ctrl_value, 1)) + t_complexity(And(ctrl_value, 1).adjoint())


def trotter(order=(0.3, 0.7), num_steps=3, **kwargs):
    return Trotterization(TROTTER_TERMS, order, TIME, num_steps, hbar=HBAR,
                          controlled_rotation_method=ANCILLA, **kwargs)


def count_bloqs(bloq, kind):
    return sum(isinstance(binst.bloq, kind) for binst in bloq.decompose_bloq().bloq_instances)


# -------------------------------------------------------------------------------------------------
# The helper in common/controls.py
# -------------------------------------------------------------------------------------------------

class TestAndAncillaRz:
    @pytest.mark.parametrize("ctrl_value", [1, 0])
    def test_is_controlled_rz_times_the_left_out_phase(self, ctrl_value):
        """Without the caller's correction it is C-Rz(theta) times exp(+i theta/2) on the control."""
        theta = 0.737
        bb = BloqBuilder()
        ctrl, target = bb.add_register('c', 1), bb.add_register('t', 1)
        target, ctrl = controls.add_and_ancilla_rz(bb, theta, target, ctrl, ctrl_value)
        cbloq = bb.finalize(c=ctrl, t=target)
        rz = np.diag([np.exp(-0.5j * theta), np.exp(0.5j * theta)])
        expected = block_matrix(rz * np.exp(0.5j * theta), (ctrl_value,))
        np.testing.assert_allclose(cbloq.tensor_contract(), expected, atol=1e-12)
        assert t_complexity(cbloq) == controls.and_ancilla_rz_t_complexity(ctrl_value)
        assert t_complexity(cbloq) == and_pair(ctrl_value) + TComplexity(rotations=1)
        assert t_complexity(cbloq).t == 4
        # 2 qubits + the ancilla at the peak
        assert _cbloq_max_width(cbloq._binst_graph,
                                lambda b: get_cost_value(b, QubitCount())) == 3


# -------------------------------------------------------------------------------------------------
# PauliStringEvolution
# -------------------------------------------------------------------------------------------------

class TestPauliStringEvolution:
    @pytest.mark.parametrize("pauli_string", PAULI_STRINGS)
    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    @pytest.mark.parametrize("is_adjoint", [False, True])
    def test_matches_block_matrix(self, pauli_string, ctrl_values, is_adjoint):
        pse = PauliStringEvolution(pauli_string, -0.62, 1.37, 1.7, is_adjoint=is_adjoint,
                                   ctrl_values=ctrl_values, controlled_rotation_method=ANCILLA)
        forward = expm(-1j * -0.62 * 1.37 / 1.7 * pauli_string_to_matrix(pauli_string))
        check_exact(pse, block_matrix(forward.conj().T if is_adjoint else forward, ctrl_values))

    @pytest.mark.parametrize("pauli_string", PAULI_STRINGS)
    @pytest.mark.parametrize("ctrl_value", [1, 0])
    @pytest.mark.parametrize("is_adjoint", [False, True])
    def test_omitted_phase(self, pauli_string, ctrl_value, is_adjoint):
        """With `omit_control_phase`, the active block is U times exp(+i * control_phase)."""
        pse = PauliStringEvolution(pauli_string, -0.62, 1.37, 1.7, is_adjoint=is_adjoint,
                                   ctrl_values=(ctrl_value,), controlled_rotation_method=ANCILLA,
                                   omit_control_phase=True)
        forward = expm(-1j * -0.62 * 1.37 / 1.7 * pauli_string_to_matrix(pauli_string))
        u = forward.conj().T if is_adjoint else forward
        phi = (-1 if is_adjoint else 1) * -0.62 * 1.37 / 1.7
        assert pse.control_phase == pytest.approx(phi)
        check_exact(pse, block_matrix(u * np.exp(1j * phi), (ctrl_value,)))

    @pytest.mark.parametrize("pauli_string", PAULI_STRINGS)
    @pytest.mark.parametrize("ctrl_values, omit", [(cvs, False) for cvs in CTRL_VALUES]
                             + [((1,), True), ((0,), True)])
    def test_count_and_call_graph_match_decomposition(self, pauli_string, ctrl_values, omit):
        pse = PauliStringEvolution(pauli_string, 0.37, 0.91, ctrl_values=ctrl_values,
                                   controlled_rotation_method=ANCILLA, omit_control_phase=omit)
        decomposition = pse.decompose_bloq()
        # Qualtran can't count an empty CompositeBloq (an all-identity term with its phase left out)
        expected = (t_complexity(decomposition) if decomposition.bloq_instances else TComplexity())
        assert t_complexity(pse) == expected
        callees = {}
        for binst in decomposition.bloq_instances:
            callees[binst.bloq] = callees.get(binst.bloq, 0) + 1
        assert dict(pse.build_call_graph(None)) == callees

    @pytest.mark.parametrize("ctrl_value", [1, 0])
    def test_hand_count(self, ctrl_value):
        """`XYZI`: 10 Cliffords of basis changes and ladder, the `And` pair, 1 rotation, and
        1 more for the phase on the control (or none, if a parent applies it)."""
        pse = PauliStringEvolution("XYZI", 0.37, 0.91, ctrl_values=(ctrl_value,),
                                   controlled_rotation_method=ANCILLA)
        base = TComplexity(clifford=10) + and_pair(ctrl_value)
        assert t_complexity(pse) == base + TComplexity(rotations=2)
        assert t_complexity(attrs.evolve(pse, omit_control_phase=True)) == base + TComplexity(
            rotations=1)

    def test_all_identity(self):
        """No rotation and no ancilla: just the phase on the control, or nothing if omitted."""
        pse = PauliStringEvolution("II", 0.37, 0.91, ctrl_values=(1,),
                                   controlled_rotation_method=ANCILLA)
        assert t_complexity(pse) == TComplexity(rotations=1)
        assert get_cost_value(pse, QubitCount()) == 3
        omitted = attrs.evolve(pse, omit_control_phase=True)
        assert t_complexity(omitted) == TComplexity()
        assert omitted.decompose_bloq().bloq_instances == set()

    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    def test_qubit_count(self, ctrl_values):
        """One rotation ancilla on top of the controls and the multi-control AND ancillas."""
        k = len(ctrl_values)
        pse = PauliStringEvolution("XIZ", 0.37, 0.91, ctrl_values=ctrl_values,
                                   controlled_rotation_method=ANCILLA)
        expected = 3 + k + (k - 1) + 1
        assert get_cost_value(pse, QubitCount()) == expected
        assert decomposition_width(pse) == expected
        if k > 1:
            assert decomposition_width(attrs.evolve(pse, ctrl_values=(1,))) == 3 + 1 + 1

    def test_uncontrolled_is_unchanged(self):
        pse = PauliStringEvolution("XYZI", 0.37, 0.91)
        with_method = attrs.evolve(pse, controlled_rotation_method=ANCILLA)
        assert t_complexity(with_method) == t_complexity(pse)
        assert get_cost_value(with_method, QubitCount()) == 4
        assert (sorted(map(str, (b.bloq for b in with_method.decompose_bloq().bloq_instances)))
                == sorted(map(str, (b.bloq for b in pse.decompose_bloq().bloq_instances))))

    def test_controlled_adjoint_power_keep_the_method(self):
        pse = PauliStringEvolution("XZ", 0.37, 0.91, controlled_rotation_method=ANCILLA)
        controlled = pse.controlled(CtrlSpec(cvs=0))
        assert controlled.controlled_rotation_method == ANCILLA
        assert controlled.adjoint() == pse.adjoint().controlled(CtrlSpec(cvs=0))
        assert (controlled ** -2).controlled_rotation_method == ANCILLA
        assert controlled != attrs.evolve(controlled, controlled_rotation_method="textbook")

    def test_uses_and_on_the_parity_qubit(self):
        pse = PauliStringEvolution("XZ", 0.37, 0.91, ctrl_values=(0,),
                                   controlled_rotation_method=ANCILLA)
        bloqs = [binst.bloq for binst in pse.decompose_bloq().bloq_instances]
        assert And(0, 1) in bloqs and And(0, 1).adjoint() in bloqs
        assert sum(isinstance(b, Rz) for b in bloqs) == 2    # the rotation and the phase


# -------------------------------------------------------------------------------------------------
# CommutingPauliStringEvolution
# -------------------------------------------------------------------------------------------------

COMMUTING_TERMS = (("XZI", 0.4), ("XIZ", -0.3), ("ZYY", 0.2), ("III", 0.7))


class TestCommutingPauliStringEvolution:
    def cpse(self, **kwargs):
        return CommutingPauliStringEvolution(COMMUTING_TERMS, 0.83, hbar=1.1,
                                             controlled_rotation_method=ANCILLA, **kwargs)

    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    @pytest.mark.parametrize("is_adjoint", [False, True])
    def test_matches_block_matrix(self, ctrl_values, is_adjoint):
        hamiltonian = sum(h * pauli_string_to_matrix(s) for s, h in COMMUTING_TERMS)
        forward = expm(-1j * hamiltonian * 0.83 / 1.1)
        check_exact(self.cpse(ctrl_values=ctrl_values, is_adjoint=is_adjoint),
                    block_matrix(forward.conj().T if is_adjoint else forward, ctrl_values))

    @pytest.mark.parametrize("ctrl_value", [1, 0])
    def test_one_merged_phase(self, ctrl_value):
        """The terms leave out their phases, and the bloq applies their sum as one phase."""
        cpse = self.cpse(ctrl_values=(ctrl_value,))
        terms = [binst.bloq for binst in cpse.decompose_bloq().bloq_instances
                 if isinstance(binst.bloq, PauliStringEvolution)]
        assert len(terms) == 4 and all(term.omit_control_phase for term in terms)
        assert cpse.control_phase == pytest.approx(sum(t.control_phase for t in terms))
        assert count_bloqs(cpse, Rz) == 1 and count_bloqs(cpse, GlobalPhase) == 1

    def test_counts(self):
        """3 rotating terms: their Cliffords, 3 `And` pairs and rotations, and 1 merged phase."""
        cpse = self.cpse(ctrl_values=(1,))
        # XZI: 2 + 2; XIZ: 2 + 2; ZYY: 8 + 4 Cliffords of basis changes and ladder
        expected = TComplexity(clifford=20) + 3 * and_pair(1) + TComplexity(rotations=3 + 1)
        assert t_complexity(cpse) == expected == t_complexity(cpse.decompose_bloq())
        assert get_cost_value(cpse, QubitCount()) == 3 + 1 + 1
        assert decomposition_width(cpse) == 5

    def test_only_identities_need_no_ancilla(self):
        cpse = CommutingPauliStringEvolution((("II", 0.3),), 0.83, ctrl_values=(1,),
                                             controlled_rotation_method=ANCILLA)
        assert get_cost_value(cpse, QubitCount()) == 3 == decomposition_width(cpse)


# -------------------------------------------------------------------------------------------------
# Trotterization
# -------------------------------------------------------------------------------------------------

class TestTrotterization:
    @pytest.mark.parametrize("order", [(1.0,), (0.5, 0.5), (0.3, 0.7)])
    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    @pytest.mark.parametrize("is_adjoint", [False, True])
    def test_matches_block_matrix(self, order, ctrl_values, is_adjoint):
        bloq = trotter(order, ctrl_values=ctrl_values, is_adjoint=is_adjoint)
        forward = reference_trotter_matrix(TROTTER_TERMS, order, TIME, 3, HBAR)
        check_exact(bloq, block_matrix(forward.conj().T if is_adjoint else forward, ctrl_values))

    @pytest.mark.parametrize("power", [2, -3])
    def test_power_matches_block_matrix(self, power):
        bloq = trotter(num_steps=2, ctrl_values=(0,)) ** power
        forward = reference_trotter_matrix(TROTTER_TERMS, (0.3, 0.7), TIME * abs(power),
                                           2 * abs(power), HBAR)
        check_exact(bloq, block_matrix(forward.conj().T if power < 0 else forward, (0,)))

    @pytest.mark.parametrize("ctrl_value", [1, 0])
    def test_exactly_one_correction(self, ctrl_value):
        """Every term leaves out its phase, through both levels; the Trotterization adds one."""
        bloq = trotter(ctrl_values=(ctrl_value,))
        terms = [binst.bloq for binst in bloq.decompose_bloq().bloq_instances
                 if isinstance(binst.bloq, CommutingPauliStringEvolution)]
        assert len(terms) == bloq.num_operations
        assert all(term.omit_control_phase for term in terms)
        assert all(pse.omit_control_phase for term in terms for pse in term._term_evolutions())
        assert count_bloqs(bloq, Rz) == 1 and count_bloqs(bloq, GlobalPhase) == 1
        assert bloq._control_phase == pytest.approx(sum(term.control_phase for term in terms))

    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    def test_count_matches_decomposition(self, ctrl_values):
        bloq = trotter(ctrl_values=ctrl_values)
        assert t_complexity(bloq) == t_complexity(bloq.decompose_bloq())

    @pytest.mark.parametrize("ctrl_value", [1, 0])
    def test_hand_count(self, ctrl_value):
        """1 rotation + an `And` pair per rotating term occurrence, plus 1 merged phase; the
        textbook method has 2 rotations + 2 CNOTs per occurrence and 1 rotation per identity
        occurrence."""
        bloq = trotter(ctrl_values=(ctrl_value,))
        textbook = attrs.evolve(bloq, controlled_rotation_method="textbook")
        uncontrolled = t_complexity(attrs.evolve(bloq, ctrl_values=()))
        occurrences = bloq._occurrences()
        rotating = sum(n for (i, _), n in occurrences.items() if TROTTER_TERMS[i][0] != "III")
        identity = sum(occurrences.values()) - rotating
        assert uncontrolled.rotations == rotating
        assert t_complexity(bloq) == (TComplexity(clifford=uncontrolled.clifford)
                                      + rotating * and_pair(ctrl_value)
                                      + TComplexity(rotations=rotating + 1))
        assert t_complexity(textbook) == TComplexity(
            clifford=uncontrolled.clifford + 2 * rotating, rotations=2 * rotating + identity)

    @pytest.mark.parametrize("ctrl_values", CTRL_VALUES)
    def test_qubit_count(self, ctrl_values):
        k = len(ctrl_values)
        bloq = trotter(ctrl_values=ctrl_values)
        assert get_cost_value(bloq, QubitCount()) == 3 + k + (k - 1) + 1
        assert decomposition_width(bloq) == 3 + k + (k - 1) + 1

    def test_call_graph_matches_decomposition(self):
        bloq = trotter(ctrl_values=(0,))
        callees = {}
        for binst in bloq.decompose_bloq().bloq_instances:
            callees[binst.bloq] = callees.get(binst.bloq, 0) + 1
        assert dict(bloq.build_call_graph(None)) == callees

    def test_uncontrolled_is_unchanged(self):
        bloq = trotter()
        textbook = attrs.evolve(bloq, controlled_rotation_method="textbook")
        assert t_complexity(bloq) == t_complexity(textbook)
        assert get_cost_value(bloq, QubitCount()) == 3
        assert np.array_equal(bloq.tensor_contract(), textbook.tensor_contract())

    def test_method_is_part_of_equality(self):
        bloq = trotter(ctrl_values=(1,))
        assert bloq != attrs.evolve(bloq, controlled_rotation_method="textbook")
        assert hash(bloq) != hash(attrs.evolve(bloq, controlled_rotation_method="textbook"))

    def test_controlled_adjoint_power_keep_the_method(self):
        bloq = trotter()
        assert bloq.controlled().controlled_rotation_method == ANCILLA
        assert bloq.controlled().adjoint() == bloq.adjoint().controlled()
        assert (bloq.controlled() ** 2).controlled_rotation_method == ANCILLA
        assert bloq.controlled().controlled().controlled_rotation_method == ANCILLA

    def test_is_cheaper_than_textbook_once_rotations_are_converted(self):
        """At the default precision a rotation is about 100 T, so (c) costs about half."""
        from qhat.analysis.resource_estimation import estimate_resources
        bloq = trotter(ctrl_values=(1,))
        textbook = attrs.evolve(bloq, controlled_rotation_method="textbook")
        cheap = estimate_resources("qualtran", bloq)["qualtran"]["T_count"]
        dear = estimate_resources("qualtran", textbook)["qualtran"]["T_count"]
        assert cheap < 0.6 * dear

    def test_large_step_count_is_fast(self):
        bloq = trotter(order=(0.5, 0.5), num_steps=5000, ctrl_values=(1, 0))
        start = time.time()
        t_complexity(bloq)
        get_cost_value(bloq, QubitCount())
        bloq.tensor_contract()
        bloq.build_call_graph(None)
        assert time.time() - start < 5

    def test_str_and_repr(self):
        assert "controlled_rotation_method='ancilla'" in repr(trotter(ctrl_values=(1,)))
        assert "controlled_rotation_method" not in repr(attrs.evolve(
            trotter(ctrl_values=(1,)), controlled_rotation_method="textbook"))


# -------------------------------------------------------------------------------------------------
# Validation
# -------------------------------------------------------------------------------------------------

class TestValidation:
    @pytest.mark.parametrize("bloq_type", ["pse", "cpse", "trotter"])
    def test_unknown_method_raises(self, bloq_type):
        make = {"pse": lambda **kw: PauliStringEvolution("XZ", 0.3, 0.7, **kw),
                "cpse": lambda **kw: CommutingPauliStringEvolution((("XZ", 0.3),), 0.7, **kw),
                "trotter": lambda **kw: Trotterization((("XZ", 0.3), ("ZI", 0.2)), (1.0,), 0.7,
                                                       2, **kw)}[bloq_type]
        with pytest.raises(ValueError, match=r"\('textbook', 'ancilla'\)"):
            make(controlled_rotation_method="toffoli")

    def test_direct_with_two_controls_raises(self):
        """(c) needs one control per rotation; "direct" keeps k controls on every rotation."""
        with pytest.raises(ValueError, match="multi_control_method 'direct'"):
            trotter(ctrl_values=(1, 0), multi_control_method="direct")
        one = trotter(ctrl_values=(1,), multi_control_method="direct")    # fine with 1 control
        with pytest.raises(ValueError, match="multi_control_method 'direct'"):
            one.controlled()

    @pytest.mark.parametrize("kwargs", [{}, {"ctrl_values": (1, 0)},
                                        {"ctrl_values": (1,),
                                         "controlled_rotation_method": "textbook"}])
    def test_omit_control_phase_needs_one_control_and_ancilla(self, kwargs):
        kwargs = {"controlled_rotation_method": ANCILLA, **kwargs}
        with pytest.raises(ValueError, match="omit_control_phase"):
            PauliStringEvolution("XZ", 0.3, 0.7, omit_control_phase=True, **kwargs)


# -------------------------------------------------------------------------------------------------
# QPE and pyLIQTR
# -------------------------------------------------------------------------------------------------

class TestQPE:
    def test_flexible_qpe(self):
        """Each controlled power keeps (c); one ancilla, reused, on top of phase + system."""
        qpe = FlexibleQPE.from_num_phase_qubits(trotter(num_steps=2), 3)
        powers = [binst.bloq for binst in qpe.decompose_bloq().bloq_instances
                  if isinstance(binst.bloq, Trotterization)]
        assert len(powers) == 3
        assert all(p.controlled_rotation_method == ANCILLA and p.ctrl_values == (1,)
                   for p in powers)
        assert get_cost_value(qpe, QubitCount()) == 3 + 3 + 1
        assert t_complexity(qpe) == t_complexity(qpe.decompose_bloq())
        np.testing.assert_allclose(qpe.decompose_bloq().tensor_contract(), qpe.tensor_contract(),
                                   atol=1e-10)


def test_pyliqtr_handles_the_measurement_based_uncompute():
    """pyLIQTR agrees with the Qualtran estimator for (c), counts and qubits."""
    from qhat.analysis.resource_estimation import estimate_resources
    for bloq in [trotter(ctrl_values=(1,)), trotter(ctrl_values=(0, 1, 1)),
                 FlexibleQPE.from_num_phase_qubits(trotter(), 2)]:
        results = estimate_resources(["qualtran", "pyliqtr"], bloq)
        for key in ("T_count", "Clifford_count", "rotation_count", "qubit_count"):
            assert results["qualtran"][key] == results["pyliqtr"][key], key
