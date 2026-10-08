"""
Tests for FlexibleQPE and the QPE window states.
"""

from functools import cached_property
from unittest import mock

import attrs
import numpy as np
import pytest

from qualtran import Bloq, QBit, Register, Signature
from qualtran.bloqs.basic_gates import CNOT, Hadamard, OnEach, ZPowGate
from qualtran.bloqs.phase_estimation import TextbookQPE
from qualtran.bloqs.qft.qft_text_book import QFTTextBook
from qualtran.cirq_interop.t_complexity_protocol import t_complexity
from qualtran.resource_counting import get_cost_value, QubitCount
from qualtran.resource_counting._qubit_counts import _cbloq_max_width

from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.qpe_window_state import RectangularWindowState
from qhat.common.trotter_flattened import Trotterization


@attrs.frozen
class TwoQubitPlainBloq(Bloq):
    """A plain Bloq (no __pow__, not a cirq.Gate) acting on two qubits."""

    @cached_property
    def signature(self):
        return Signature([Register('a', QBit()), Register('b', QBit())])

    def build_composite_bloq(self, bb, a, b):
        a = bb.add(ZPowGate(exponent=0.37), q=a)
        a, b = bb.add(CNOT(), ctrl=a, target=b)
        b = bb.add(ZPowGate(exponent=0.11), q=b)
        return {'a': a, 'b': b}


def small_trotterization():
    return Trotterization.from_method(
        pauli_terms=[("XZ", 0.5), ("ZY", 0.3)], method="second order", time=0.3, num_steps=2
    )


def generic_tensor_contract(bloq):
    """Contract the decomposition, bypassing the optimized tensor_contract."""
    return bloq.decompose_bloq().tensor_contract()


def phase_blocks(qpe):
    """Return the target-register blocks W[k] of the controlled-power part of a FlexibleQPE."""
    return qpe._phase_blocks()


# -------------------------------------------------------------------------------------------------
# Window states
# -------------------------------------------------------------------------------------------------

class TestRectangularWindowState:

    @pytest.mark.parametrize("m", [1, 2, 4])
    def test_tensor_is_hadamards(self, m):
        expected = OnEach(m, Hadamard()).tensor_contract()
        np.testing.assert_allclose(RectangularWindowState(m).tensor_contract(), expected)

    def test_prepares_uniform_state(self):
        m = 4
        state = RectangularWindowState(m).tensor_contract()[:, 0]
        np.testing.assert_allclose(state, np.full(2**m, 2**(-m / 2)))

    def test_signature(self):
        reg = RectangularWindowState(3).signature[0]
        assert reg.name == 'qpe_reg'
        assert reg.total_bits() == 3

    def test_from_precision_and_delta(self):
        # m = n + ceil(log2(2 + 1/(2*delta)))
        assert RectangularWindowState.from_precision_and_delta(5, 0.01).m_bits == 5 + 6

    def test_from_standard_deviation_eps(self):
        # m = ceil(2*log2(pi/eps))
        assert RectangularWindowState.from_standard_deviation_eps(0.1).m_bits == 10

    @pytest.mark.parametrize("phase_error, p_fail, expected", [
        (1 / 8, 0.1, 3 + 3),      # 2^-3 = 1/8 exactly; ceil(log2(2 + 5)) = 3
        (0.1, 0.1, 4 + 3),        # 2^-4 < 0.1 < 2^-3
        (1 / 1024, 0.01, 10 + 6),  # ceil(log2(2 + 50)) = 6
    ])
    def test_from_requirements(self, phase_error, p_fail, expected):
        assert RectangularWindowState.from_requirements(phase_error, p_fail).m_bits == expected

    def test_from_requirements_tolerates_round_off(self):
        """A phase error of 2^-n (up to round-off) needs exactly n precision bits."""
        phase_error = (1 / 3) * (3 / 2**7)
        assert RectangularWindowState.from_requirements(phase_error, 0.5).m_bits == 7 + 2

    @pytest.mark.parametrize("phase_error, p_fail", [(0.0, 0.1), (1.0, 0.1), (0.1, 0.0), (0.1, 1.0)])
    def test_from_requirements_rejects_out_of_range(self, phase_error, p_fail):
        with pytest.raises(ValueError):
            RectangularWindowState.from_requirements(phase_error, p_fail)

    def test_from_num_phase_qubits(self):
        assert RectangularWindowState.from_num_phase_qubits(5) == RectangularWindowState(5)

    @pytest.mark.parametrize("phase_error, p_fail", [(1 / 8, 0.2), (1 / 16, 0.05)])
    def test_from_requirements_meets_requirements(self, phase_error, p_fail):
        """Worst-case textbook QPE (phase halfway between grid points) meets the requirement."""
        prep = RectangularWindowState.from_requirements(phase_error, p_fail)
        m = prep.m_bits
        phi = 0.5 / 2**m
        z = ZPowGate(exponent=2 * phi)
        mat = FlexibleQPE(z, m, prep).tensor_contract()
        psi_in = np.zeros(2**(m + 1))
        psi_in[1] = 1.0
        probs = np.abs(mat @ psi_in).reshape(2**m, 2)[:, 1] ** 2
        estimates = np.arange(2**m) / 2**m
        dist = np.abs(estimates - phi)
        dist = np.minimum(dist, 1 - dist)
        assert probs[dist > phase_error].sum() <= p_fail


# -------------------------------------------------------------------------------------------------
# FlexibleQPE construction
# -------------------------------------------------------------------------------------------------

class TestConstruction:

    def test_defaults_and_signature(self):
        u = small_trotterization()
        qpe = FlexibleQPE(u, 3)
        assert qpe.num_ancilla_qubits == 3
        assert qpe.ancilla_prep == RectangularWindowState(3)
        assert qpe.qft_inv == QFTTextBook(3, with_reverse=True).adjoint()
        assert [reg.name for reg in qpe.signature] == ['qpe_reg', 'q']
        assert qpe.signature.n_qubits() == 3 + 2
        assert str(qpe) == 'FlexibleQPE[3]'

    def test_qubit_counts(self):
        qpe = FlexibleQPE(TwoQubitPlainBloq(), 3)
        assert qpe.num_ancilla_qubits == 3
        assert qpe.num_state_qubits == 2
        assert qpe.num_total_qubits == qpe.signature.n_qubits() == 5

    def test_num_ancilla_qubits_keyword(self):
        u = small_trotterization()
        assert FlexibleQPE(u, num_ancilla_qubits=3) == FlexibleQPE(u, 3)

    def test_explicit_components_equal_defaults(self):
        u = small_trotterization()
        assert FlexibleQPE(u, 3) == FlexibleQPE(
            u, 3, RectangularWindowState(3), QFTTextBook(3, with_reverse=True).adjoint())

    def test_mismatched_ancilla_prep_raises(self):
        with pytest.raises(ValueError, match="ancilla_prep"):
            FlexibleQPE(ZPowGate(exponent=0.1), 3, ancilla_prep=RectangularWindowState(2))

    def test_mismatched_qft_inv_raises(self):
        with pytest.raises(ValueError, match="qft_inv"):
            FlexibleQPE(ZPowGate(exponent=0.1), 3, qft_inv=QFTTextBook(2).adjoint())

    def test_unitary_power_fast_forwards_when_available(self):
        qpe = FlexibleQPE(small_trotterization(), 3)
        assert qpe.unitary_power(1) is qpe.unitary
        assert isinstance(qpe.unitary_power(4), Trotterization)
        assert qpe.unitary_power(4).num_steps == 4 * qpe.unitary.num_steps

    def test_unitary_power_falls_back_to_repetition(self):
        from qualtran.bloqs.basic_gates import Power
        qpe = FlexibleQPE(TwoQubitPlainBloq(), 3)
        assert qpe.unitary_power(4) == Power(TwoQubitPlainBloq(), 4)

    def test_from_requirements_dispatches_to_window_state(self):
        u = small_trotterization()
        qpe = FlexibleQPE.from_requirements(u, 0.1, 0.1)
        prep = RectangularWindowState.from_requirements(0.1, 0.1)
        assert qpe == FlexibleQPE(u, prep.m_bits, prep)

    def test_from_requirements_uses_given_window_state_class(self):
        calls = []

        @attrs.frozen
        class RecordingWindowState(RectangularWindowState):
            @classmethod
            def from_requirements(cls, phase_error, probability_of_failure):
                calls.append((phase_error, probability_of_failure))
                return cls(2)

        qpe = FlexibleQPE.from_requirements(
            small_trotterization(), 0.1, 0.2, ancilla_prep=RecordingWindowState)
        assert calls == [(0.1, 0.2)]
        assert qpe.ancilla_prep == RecordingWindowState(2)
        assert qpe.num_ancilla_qubits == 2

    def test_from_requirements_qft_inv_factory_receives_register_size(self):
        sizes = []

        def factory(m):
            sizes.append(m)
            return QFTTextBook(m, with_reverse=True).adjoint()

        qpe = FlexibleQPE.from_requirements(small_trotterization(), 1 / 8, 0.1, qft_inv=factory)
        assert sizes == [qpe.num_ancilla_qubits] == [6]

    def test_from_num_phase_qubits(self):
        u = small_trotterization()
        assert FlexibleQPE.from_num_phase_qubits(u, 4) == FlexibleQPE(u, 4)


# -------------------------------------------------------------------------------------------------
# Tensor contraction
# -------------------------------------------------------------------------------------------------

class TestTensorContraction:

    @pytest.mark.parametrize("m", [1, 2, 3])
    def test_matches_qualtran_textbook_qpe(self, m):
        """For a cirq-compatible unitary, 0.4.0's TextbookQPE works and must agree."""
        z = ZPowGate(exponent=2 * 0.234)
        np.testing.assert_allclose(
            FlexibleQPE(z, m).tensor_contract(),
            TextbookQPE(z, m).tensor_contract(),
            atol=1e-12,
        )

    @pytest.mark.parametrize("unitary", [
        ZPowGate(exponent=2 * 0.234),
        small_trotterization(),
        TwoQubitPlainBloq(),
    ], ids=["zpow", "trotterization", "plain-bloq"])
    @pytest.mark.parametrize("m", [1, 3])
    def test_optimized_matches_decomposition(self, unitary, m):
        qpe = FlexibleQPE(unitary, m)
        np.testing.assert_allclose(
            qpe.tensor_contract(), generic_tensor_contract(qpe), atol=1e-12)

    def test_trotterization_tensor_contract_works(self):
        """Qualtran 0.4.0's TextbookQPE fails here (cirq.pow on a non-cirq Bloq)."""
        qpe = FlexibleQPE(small_trotterization(), 3)
        mat = qpe.tensor_contract()
        assert mat.shape == (2**5, 2**5)
        np.testing.assert_allclose(mat.conj().T @ mat, np.eye(2**5), atol=1e-10)

    def test_plain_bloq_blocks_are_unitary_powers(self):
        qpe = FlexibleQPE(TwoQubitPlainBloq(), 3)
        u = TwoQubitPlainBloq().tensor_contract()
        for k, block in enumerate(phase_blocks(qpe)):
            np.testing.assert_allclose(block, np.linalg.matrix_power(u, k), atol=1e-12)

    @pytest.mark.parametrize("combine_terms", [True, False])
    def test_trotterization_unitary_powers_are_matrix_powers(self, combine_terms):
        """tensor_contract squares U's matrix, so Trotterization.__pow__ must give exactly U^k."""
        u = attrs.evolve(small_trotterization(), combine_terms=combine_terms)
        qpe = FlexibleQPE(u, 3)
        umat = u.tensor_contract()
        for j in range(3):
            np.testing.assert_allclose(
                qpe.unitary_power(2**j).tensor_contract(),
                np.linalg.matrix_power(umat, 2**j), atol=1e-10)

    def test_contracts_unitary_once(self):
        """The U^(2^j) bloqs are not contracted; their matrices come from squaring U."""
        qpe = FlexibleQPE(small_trotterization(), 4)
        with mock.patch.object(
                Trotterization, 'tensor_contract', autospec=True,
                side_effect=Trotterization.tensor_contract) as contract:
            qpe.tensor_contract()
        assert contract.call_count == 1

    def test_nested_uses_custom_tensor(self):
        """Inside a larger bloq, FlexibleQPE adds its dense matrix rather than decomposing."""
        qpe = FlexibleQPE(TwoQubitPlainBloq(), 3)
        with mock.patch.object(FlexibleQPE, 'decompose_bloq',
                               side_effect=AssertionError("decomposed")):
            nested = qpe.as_composite_bloq().tensor_contract()
        np.testing.assert_allclose(nested, generic_tensor_contract(qpe), atol=1e-12)

    def test_trotterization_matches_unoptimized_reference(self):
        """Full QPE matrix vs. (iQFT ⊗ I) · diag_k(U^k) · (H^⊗m ⊗ I) built from U's matrix."""
        m = 3
        u = small_trotterization()
        umat = u.tensor_contract()
        d = umat.shape[0]
        blocks = np.zeros((2**m * d, 2**m * d), dtype=complex)
        for k in range(2**m):
            blocks[k*d:(k+1)*d, k*d:(k+1)*d] = np.linalg.matrix_power(umat, k)
        prep = np.kron(OnEach(m, Hadamard()).tensor_contract(), np.eye(d))
        qft_inv = np.kron(QFTTextBook(m, with_reverse=True).adjoint().tensor_contract(), np.eye(d))
        np.testing.assert_allclose(
            FlexibleQPE(u, m).tensor_contract(),
            qft_inv @ blocks @ prep, atol=1e-10)

    @pytest.mark.parametrize("k", [0, 1, 5, 7])
    def test_estimates_exact_phase(self, k):
        """An eigenphase of exactly k/2^m is measured as k with certainty."""
        m = 3
        z = ZPowGate(exponent=2 * k / 2**m)  # |1> has eigenvalue exp(2 pi i k / 2^m)
        mat = FlexibleQPE(z, m).tensor_contract()
        psi_in = np.zeros(2**(m + 1))
        psi_in[1] = 1.0  # qpe_reg = |0>, target = |1>
        probs = np.abs(mat @ psi_in).reshape(2**m, 2)[:, 1] ** 2
        assert probs[k] == pytest.approx(1.0)


# -------------------------------------------------------------------------------------------------
# Resource estimation
# -------------------------------------------------------------------------------------------------

class TestResourceEstimation:

    def test_call_graph_has_one_entry_per_phase_bit(self):
        u = small_trotterization()
        qpe = FlexibleQPE(u, 3)
        counts = dict(qpe.build_call_graph(None))
        assert counts[qpe.ancilla_prep] == 1
        assert counts[qpe.qft_inv] == 1
        for j in range(3):
            assert counts[(u ** 2**j).controlled()] == 1

    def test_t_complexity_without_fast_forwarding(self):
        u = TwoQubitPlainBloq()
        m = 3
        qpe = FlexibleQPE(u, m)
        expected = (t_complexity(RectangularWindowState(m))
                    + (2**m - 1) * t_complexity(u.controlled())
                    + t_complexity(qpe.qft_inv))
        assert t_complexity(qpe) == expected
        assert t_complexity(qpe) == t_complexity(TextbookQPE(u, m))

    def test_t_complexity_with_fast_forwarding(self):
        u = small_trotterization()
        m = 3
        qpe = FlexibleQPE(u, m)
        expected = t_complexity(RectangularWindowState(m)) + t_complexity(qpe.qft_inv)
        for j in range(m):
            expected = expected + t_complexity((u ** 2**j).controlled())
        assert t_complexity(qpe) == expected

    @pytest.mark.parametrize("unitary", [small_trotterization(), TwoQubitPlainBloq()],
                             ids=["trotterization", "plain-bloq"])
    def test_static_qubit_count_matches_decomposition(self, unitary):
        qpe = FlexibleQPE(unitary, 3)
        from_decomposition = _cbloq_max_width(
            qpe.decompose_bloq()._binst_graph, lambda b: get_cost_value(b, QubitCount()))
        assert get_cost_value(qpe, QubitCount()) == from_decomposition == 3 + 2
