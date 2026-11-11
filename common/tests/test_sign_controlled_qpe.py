"""Sign-controlled QPE, option (d) of plan T11.

Phase bit j applies V^dagger on |0> and V on |1>, with V^2 = U^(2^j), instead of controlled
U^(2^j). V is a Trotterization with a Z on the phase qubit added to every term. Every matrix
check is against matrices built here from `expm`, with terms that don't commute.
"""

import time

import attrs
import numpy as np
import pytest
from scipy.linalg import expm

from qualtran.bloqs.bookkeeping import Join, Split
from qualtran.cirq_interop.t_complexity_protocol import t_complexity
from qualtran.resource_counting import QubitCount, get_cost_value
from qualtran.resource_counting._qubit_counts import _cbloq_max_width

from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.pauli_utils import pauli_string_to_matrix
from qhat.common.trotter_flattened import Trotterization
from qhat.common.tests.test_multi_control import block_matrix

SIGN = "sign"
# Terms that don't commute, with an all-identity term; odd time and hbar
TERMS = (("XZI", 0.4), ("ZIY", -0.3), ("YXX", 0.25), ("III", 0.6))
TIME, HBAR = 0.43, 0.9
SYMMETRIC = ["second order", "fourth order"]


def trotter(method="fourth order", num_steps=4, terms=TERMS, **kwargs):
    return Trotterization.from_method(terms, method, TIME, num_steps, hbar=HBAR, **kwargs)


def reference(tr, sign=1):
    """exp(-i sign h c dt P / hbar) for each term of the expanded sequence, built here."""
    dt = tr.time / tr.num_steps
    out = np.eye(2**tr.num_qubits, dtype=complex)
    for i, c in tr.expanded_sequence:
        p, h = tr.pauli_terms[i]
        out = expm(-1j * sign * h * c * dt / tr.hbar * pauli_string_to_matrix(p)) @ out
    return out.conj().T if tr.is_adjoint else out


def phase_distribution(matrix, m, psi):
    """Probabilities of each phase-register outcome, for phase register |0> and target psi."""
    d = len(psi)
    state = np.zeros(2**m * d, dtype=complex)
    state[:d] = psi
    return (np.abs((matrix @ state).reshape(2**m, d)) ** 2).sum(axis=1)


# -------------------------------------------------------------------------------------------------
# The Trotterization pieces
# -------------------------------------------------------------------------------------------------

class TestIsSymmetric:
    @pytest.mark.parametrize("method, expected", [
        ("first order", False), ("second order", True), ("third order", False),
        ("fourth order", True), ("bm4", True), ("opt4", False), ((0.3, 0.7), False)])
    @pytest.mark.parametrize("num_steps", [1, 3])
    def test_methods(self, method, expected, num_steps):
        if isinstance(method, tuple):
            tr = Trotterization(TERMS, method, TIME, num_steps)
        else:
            tr = trotter(method, num_steps)
        assert tr.is_symmetric is expected
        # Symmetric means U(-t) is U^dagger
        negated = np.allclose(reference(attrs.evolve(tr, time=-TIME)), reference(tr).conj().T,
                              atol=1e-12)
        assert negated is expected


class TestSquareRoot:
    @pytest.mark.parametrize("num_steps", [2, 6])
    def test_squares_to_the_original(self, num_steps):
        u = trotter(num_steps=num_steps).adjoint()
        root = u.square_root()
        assert root.num_steps == num_steps // 2 and root.time == TIME / 2 and root.is_adjoint
        assert root**2 == u
        np.testing.assert_allclose(root.tensor_contract() @ root.tensor_contract(),
                                   reference(u), atol=1e-12)

    def test_odd_step_count_has_none(self):
        assert trotter(num_steps=3).square_root() is None


class TestSignControlled:
    @pytest.mark.parametrize("method", SYMMETRIC)
    @pytest.mark.parametrize("num_steps", [1, 3])
    @pytest.mark.parametrize("adjoint", [False, True])
    def test_blocks_are_v_dagger_and_v(self, method, num_steps, adjoint):
        """Exact: V^dagger on control |0>, V on |1>, for the tensor and the decomposition."""
        v = trotter(method, num_steps)
        v = v.adjoint() if adjoint else v
        expected = np.zeros((2 * 2**v.num_qubits,) * 2, dtype=complex)
        d = 2**v.num_qubits
        expected[:d, :d] = reference(v).conj().T
        expected[d:, d:] = reference(v)
        s = v.sign_controlled()
        np.testing.assert_allclose(s.tensor_contract(), expected, atol=1e-12)
        np.testing.assert_allclose(s.decompose_bloq().tensor_contract(), expected, atol=1e-12)

    def test_is_not_controlled_v(self):
        """`.controlled()` is identity on |0>; the sign-controlled bloq isn't."""
        v = trotter(num_steps=2)
        assert v.controlled() != v.sign_controlled()
        np.testing.assert_allclose(v.controlled().tensor_contract(),
                                   block_matrix(reference(v), (1,)), atol=1e-12)
        d = 2**v.num_qubits
        assert not np.allclose(v.sign_controlled().tensor_contract()[:d, :d], np.eye(d))

    def test_terms_and_counts(self):
        """Each term gets a Z on the new first qubit and the opposite sign: 1 rotation and 2
        more CNOTs than the term; an all-identity term becomes 1 rotation on the control."""
        v = trotter(num_steps=3)
        s = v.sign_controlled()
        assert s.pauli_terms == tuple(('Z' + p, -h) for p, h in TERMS)
        assert s.num_qubits == v.num_qubits + 1
        occurrences = len(v.expanded_sequence)
        rotating = sum(1 for i, _ in v.expanded_sequence if set(TERMS[i][0]) != {'I'})
        tc, tc_v = t_complexity(s), t_complexity(v)
        assert tc.t == 0
        assert tc.rotations == occurrences
        assert tc.clifford == tc_v.clifford + 2 * rotating
        assert tc == t_complexity(s.decompose_bloq())
        assert get_cost_value(s, QubitCount()) == v.num_qubits + 1

    @pytest.mark.parametrize("method", ["first order", "third order"])
    def test_non_symmetric_raises(self, method):
        with pytest.raises(ValueError, match="symmetric"):
            trotter(method, 2).sign_controlled()

    def test_controlled_raises(self):
        with pytest.raises(ValueError, match="uncontrolled"):
            trotter().controlled().sign_controlled()


# -------------------------------------------------------------------------------------------------
# FlexibleQPE with controlled_power_method = "sign"
# -------------------------------------------------------------------------------------------------

def sign_qpe(u, m):
    return FlexibleQPE(u, m, controlled_power_method=SIGN)


class TestSignQPE:
    @pytest.mark.parametrize("method", SYMMETRIC)
    @pytest.mark.parametrize("num_steps", [2, 4])
    @pytest.mark.parametrize("adjoint", [False, True])
    @pytest.mark.parametrize("m", [1, 3])
    def test_same_phase_distribution_as_controlled(self, method, num_steps, adjoint, m):
        """The phase register reads out the same distribution as textbook QPE, for random states."""
        u = trotter(method, num_steps)
        u = u.adjoint() if adjoint else u
        standard = FlexibleQPE(u, m).tensor_contract()
        sign = sign_qpe(u, m).tensor_contract()
        rng = np.random.default_rng(11)
        for _ in range(4):
            psi = rng.normal(size=2**u.num_qubits) + 1j * rng.normal(size=2**u.num_qubits)
            psi /= np.linalg.norm(psi)
            np.testing.assert_allclose(phase_distribution(sign, m, psi),
                                       phase_distribution(standard, m, psi), atol=1e-12)

    @pytest.mark.parametrize("k", [1, 3, 6])
    def test_reads_out_an_exact_phase(self, k):
        """|00> of exp(-i (a ZI + b IZ) t) has phase k/8 turns, which is measured with certainty.
        (These terms commute, so the phase is exact; ordering is tested above.)"""
        m = 3
        a = 0.7
        b = (-2 * np.pi * k / 2**m) * HBAR / TIME - a     # (a + b) t / hbar = -2 pi k / 8
        u = trotter(num_steps=2, terms=(("ZI", a), ("IZ", b)))
        psi = np.zeros(4)
        psi[0] = 1
        probs = phase_distribution(sign_qpe(u, m).tensor_contract(), m, psi)
        assert probs[k] == pytest.approx(1.0)

    def test_phase_blocks(self):
        """W[k] = U^k (R^dagger)^(2^m - 1), with R the square root of U, built here."""
        m = 3
        u = trotter(num_steps=2)
        r = reference(u.square_root())
        w = sign_qpe(u, m)._phase_blocks()
        for k in range(2**m):
            expected = (np.linalg.matrix_power(r, 2 * k)
                        @ np.linalg.matrix_power(r.conj().T, 2**m - 1))
            np.testing.assert_allclose(w[k], expected, atol=1e-11)

    @pytest.mark.parametrize("m", [1, 2, 3])
    def test_decomposition_matches_tensor(self, m):
        qpe = sign_qpe(trotter(num_steps=2).adjoint(), m)
        np.testing.assert_allclose(qpe.decompose_bloq().tensor_contract(), qpe.tensor_contract(),
                                   atol=1e-11)

    def test_decomposition_and_counts(self):
        """Each phase bit applies the sign-controlled square root of U^(2^j): no T, no ancilla,
        1 rotation per term occurrence."""
        m = 3
        u = trotter(num_steps=2)
        qpe = sign_qpe(u, m)
        powers = sorted((binst.bloq for binst in qpe.decompose_bloq().bloq_instances
                         if isinstance(binst.bloq, Trotterization)), key=lambda b: b.num_steps)
        assert powers == [(u**2**j).square_root().sign_controlled() for j in range(m)]
        assert [p.num_steps for p in powers] == [1, 2, 4]
        expected = t_complexity(qpe.ancilla_prep) + t_complexity(qpe.qft_inv)
        for p in powers:
            expected += t_complexity(p)
        assert t_complexity(qpe) == expected == t_complexity(qpe.decompose_bloq())
        assert t_complexity(qpe).t == t_complexity(qpe.qft_inv).t
        counts = dict(qpe.build_call_graph(None))
        assert all(counts[p] == 1 for p in powers)
        # FlexibleQPE's call graph leaves out the decomposition's Split and Join
        assert counts == {b: n for b, n in qpe.decompose_bloq().build_call_graph(None)
                          if not isinstance(b, (Split, Join))}

    def test_quarter_of_the_rotations(self):
        """Per term occurrence of U^(2^j): (b) 2 rotations, (d) 1 rotation on half as many."""
        m = 3
        u = trotter(num_steps=4, terms=TERMS[:3])    # no all-identity term
        qft = t_complexity(FlexibleQPE(u, m).qft_inv).rotations
        textbook = t_complexity(FlexibleQPE(u, m)).rotations - qft
        sign = t_complexity(sign_qpe(u, m)).rotations - qft
        occurrences = sum(len((u**2**j).expanded_sequence) for j in range(m))
        assert textbook == 2 * occurrences
        assert sign == sum(len((u**2**j).square_root().expanded_sequence) for j in range(m))
        assert sign < 0.27 * textbook

    def test_qubit_count(self):
        m = 3
        qpe = sign_qpe(trotter(num_steps=2), m)
        width = _cbloq_max_width(qpe.decompose_bloq()._binst_graph,
                                 lambda b: get_cost_value(b, QubitCount()))
        assert get_cost_value(qpe, QubitCount()) == width == m + 3

    def test_large_step_count_is_fast(self):
        qpe = sign_qpe(trotter(num_steps=5000), 6)
        start = time.time()
        tc = t_complexity(qpe)
        qubits = get_cost_value(qpe, QubitCount())
        qpe.tensor_contract()
        assert time.time() - start < 10
        assert tc.rotations > 0 and qubits == 6 + 3

    def test_default_is_controlled(self):
        u = trotter()
        assert FlexibleQPE(u, 3) == FlexibleQPE(u, 3, controlled_power_method="controlled")
        assert FlexibleQPE(u, 3) != sign_qpe(u, 3)
        assert str(sign_qpe(u, 3)) == "FlexibleQPE[3, sign]"

    def test_pyliqtr_agrees(self):
        from qhat.analysis.resource_estimation import estimate_resources
        results = estimate_resources(["qualtran", "pyliqtr"], sign_qpe(trotter(num_steps=2), 2))
        for key in ("T_count", "Clifford_count", "rotation_count", "qubit_count"):
            assert results["qualtran"][key] == results["pyliqtr"][key], key


class TestSignQPEValidation:
    def test_unknown_method_raises(self):
        with pytest.raises(ValueError, match="controlled_power_method"):
            FlexibleQPE(trotter(), 3, controlled_power_method="signed")

    def test_odd_step_count_raises(self):
        with pytest.raises(ValueError, match="even number of Trotter steps"):
            sign_qpe(trotter(num_steps=3), 3)

    def test_non_symmetric_raises(self):
        with pytest.raises(ValueError, match="symmetric"):
            sign_qpe(trotter("first order", num_steps=2), 3)

    def test_other_unitary_raises(self):
        from qhat.common.tests.test_flexible_qpe import TwoQubitPlainBloq
        with pytest.raises(ValueError, match="flattened Trotterization"):
            sign_qpe(TwoQubitPlainBloq(), 3)
