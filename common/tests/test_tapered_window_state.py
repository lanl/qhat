"""
Tests for the tapered (DPSS-family) QPE window states.

Two definitions of failure appear here:
- QHAT's (used for sizing): the estimate k/N is more than delta = 2^-(l+1) from the phase.
- The paper's (eq. 15, used by its theorems): the estimate is not one of the 2K + 1 closest to the
  phase. Those 2K + 1 are all within delta, so QHAT's failure is at most the paper's.
"""

from collections import Counter
from unittest import mock

import numpy as np
import pytest

from qualtran.bloqs.basic_gates import CNOT, Ry, ZPowGate
from qualtran.bloqs.qft.qft_text_book import QFTTextBook
from qualtran.cirq_interop.t_complexity_protocol import t_complexity, TComplexity
from qualtran.resource_counting import get_cost_value, QubitCount, SympySymbolAllocator

from qhat.common.exact_state_preparation import ExactStatePreparation
from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.qpe_window_state import precision_bits, RectangularWindowState
from qhat.common.tapered_window_state import (
    dpss_average_failure_bound,
    DPSSWindowState,
    half_bandwidth,
)

TAPERED = [DPSSWindowState]
TAPERED_IDS = ["dpss"]


def window(prep):
    """The state the window circuit prepares from |0...0>, contracted gate by gate."""
    return prep.decompose_bloq().tensor_contract()[:, 0]


def outcome_probs(w, theta):
    """QPE outcome probabilities for window w and phase theta: |QFT^-1 (w e^{2 pi i n theta})|^2."""
    n = np.arange(len(w))
    return np.abs(np.fft.fft(w * np.exp(2j * np.pi * n * theta)))**2 / len(w)


def distances(num_outcomes, theta):
    """Circular distance (in turns) from each estimate k/N to theta."""
    d = np.abs(np.arange(num_outcomes) / num_outcomes - theta) % 1
    return np.minimum(d, 1 - d)


def qhat_failure(w, theta, phase_error):
    probs = outcome_probs(w, theta)
    return probs[distances(len(w), theta) > phase_error].sum()


def paper_failure(w, theta, k_max):
    """1 - probability of the 2K + 1 estimates closest to theta (paper eq. 15)."""
    probs = outcome_probs(w, theta)
    closest = np.argsort(distances(len(w), theta), kind="stable")[:2 * k_max + 1]
    return 1 - probs[closest].sum()


def thetas(num_outcomes, count=33):
    """Phases spread over one grid spacing (failure has period 1/N in theta), off the grid."""
    return (np.arange(count) + 0.37) / (count * num_outcomes)


def sinc_kernel(n, k_max):
    """C[a, b] = sin(pi (a - b) (2K + 1)/N) / (pi (a - b)): w^T C w is the energy of w at
    frequencies |f| <= (2K + 1)/(2N), i.e. the paper's average success probability."""
    d = np.subtract.outer(np.arange(n), np.arange(n))
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.sin(np.pi * d * (2 * k_max + 1) / n) / (np.pi * d)
    c[d == 0] = (2 * k_max + 1) / n
    return c


def average_paper_failure(prep):
    w = prep.taper()
    return 1 - w @ sinc_kernel(len(w), prep.half_bandwidth) @ w


# -------------------------------------------------------------------------------------------------
# Tapers
# -------------------------------------------------------------------------------------------------

class TestTapers:

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("l, m", [(0, 1), (2, 2), (1, 3), (3, 3), (2, 4)])
    def test_real_unit_norm_symmetric_positive(self, cls, l, m):
        w = cls(l, m).taper()
        assert w.shape == (2**(l + m),) and np.isrealobj(w)
        assert np.linalg.norm(w) == pytest.approx(1.0)
        np.testing.assert_allclose(w, w[::-1], atol=1e-14)
        assert w.sum() > 0

    @pytest.mark.parametrize("l, m", [(0, 2), (2, 2), (1, 3), (3, 3)])
    def test_dpss_is_top_eigenvector_of_sinc_kernel(self, l, m):
        """The DPSS maximizes the energy in |f| <= (2K + 1)/(2N) (paper Sec. IV C)."""
        _, vectors = np.linalg.eigh(sinc_kernel(2**(l + m), half_bandwidth(m)))
        top = vectors[:, -1] * np.sign(vectors[:, -1].sum())
        np.testing.assert_allclose(DPSSWindowState(l, m).taper(), top, atol=1e-8)


# -------------------------------------------------------------------------------------------------
# Circuits
# -------------------------------------------------------------------------------------------------

class TestCircuits:

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("l, m", [(0, 1), (2, 1), (0, 2), (1, 2), (2, 3), (3, 3), (1, 4)])
    def test_circuit_prepares_taper(self, cls, l, m):
        """Including the phase (the taper is real and positive-sum, so no global phase)."""
        prep = cls(l, m)
        np.testing.assert_allclose(window(prep), prep.taper(), atol=1e-13)

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    def test_signature(self, cls):
        prep = cls(2, 3)
        (reg,) = prep.signature
        assert reg.name == 'qpe_reg' and reg.total_bits() == 5
        assert prep.m_bits == 5 and prep.half_bandwidth == 3

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("theta", [0.3141, 0.8])
    def test_flexible_qpe_outcomes_match_window(self, cls, theta):
        """In FlexibleQPE, the outcome distribution is that of the window (eq. 13)."""
        prep = cls(2, 3)
        mat = FlexibleQPE(ZPowGate(exponent=2 * theta), 5, prep).tensor_contract()
        probs = np.abs(mat[:, 1].reshape(32, 2)[:, 1])**2  # qpe_reg = |0>, target = |1>
        np.testing.assert_allclose(probs, outcome_probs(prep.taper(), theta), atol=1e-13)


# -------------------------------------------------------------------------------------------------
# Failure probabilities
# -------------------------------------------------------------------------------------------------

class TestFailureProbability:

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("phase_error, p_fail", [(0.3, 0.1), (0.1, 1e-3), (0.07, 1e-8)])
    def test_from_requirements_meets_requirements_in_worst_case(self, cls, phase_error, p_fail):
        """Over phases across a grid spacing, Pr[|estimate - phase| > delta] <= p_fail, where
        delta = 2^-(l+1) <= phase_error (so the requirement itself holds too)."""
        prep = cls.from_requirements(phase_error, p_fail)
        delta = 2**-(prep.num_precision_qubits + 1)
        assert delta <= phase_error
        w = window(prep)
        worst = max(qhat_failure(w, t, delta) for t in thetas(len(w)))
        assert worst <= p_fail

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("l, m", [(0, 3), (1, 3), (3, 3), (1, 4), (3, 4)])
    def test_worst_case_meets_sizing_bound(self, cls, l, m):
        """Sizing uses factor * Thm 1, a bound on the average; check it holds for the worst case."""
        prep = cls(l, m)
        w = window(prep)
        worst = max(qhat_failure(w, t, 2**-(l + 1)) for t in thetas(len(w)))
        assert worst <= cls.failure_bound_factor * dpss_average_failure_bound(l, m)

    @pytest.mark.parametrize("l, m", [(0, 2), (2, 2), (1, 3), (3, 3)])
    def test_qhat_failure_at_most_paper_failure(self, l, m):
        w = DPSSWindowState(l, m).taper()
        for t in thetas(len(w), 9):
            assert qhat_failure(w, t, 2**-(l + 1)) <= paper_failure(w, t, half_bandwidth(m)) + 1e-15

    @pytest.mark.parametrize("l, m", [(0, 3), (1, 3), (3, 3), (1, 4), (6, 4)])
    def test_theorem_1_bounds_dpss_average(self, l, m):
        assert average_paper_failure(DPSSWindowState(l, m)) <= dpss_average_failure_bound(l, m)


# -------------------------------------------------------------------------------------------------
# Sizing
# -------------------------------------------------------------------------------------------------

class TestSizing:

    @pytest.mark.parametrize("m, k", [(1, 0), (2, 1), (3, 3), (4, 7), (6, 31)])
    def test_half_bandwidth(self, m, k):
        assert half_bandwidth(m) == k

    def test_theorem_1_example(self):
        # min(8 exp(-29 / ((2/pi^2) ln 1024)), 10 exp(-24 / ((2/pi^2) ln 1575))) = 8.64e-9
        assert dpss_average_failure_bound(3, 5) == pytest.approx(8.640019652e-9)

    def test_theorem_1_decreases_with_confidence_qubits(self):
        bounds = [dpss_average_failure_bound(3, m) for m in range(1, 9)]
        assert all(a > b for a, b in zip(bounds, bounds[1:]))

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("l", [0, 3, 10])
    @pytest.mark.parametrize("p_fail", [0.3, 1e-3, 1e-12])
    def test_confidence_qubits_is_smallest_meeting_bound(self, cls, l, p_fail):
        m = cls.confidence_qubits(l, p_fail)
        factor = cls.failure_bound_factor
        assert factor * dpss_average_failure_bound(l, m) <= p_fail
        assert m == 1 or factor * dpss_average_failure_bound(l, m - 1) > p_fail

    def test_confidence_qubits_example(self):
        assert [cls.confidence_qubits(3, 1e-8) for cls in TAPERED] == [5]
        assert [cls.failure_bound_factor for cls in TAPERED] == [1]

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("phase_error, l", [
        (0.75, 0), (0.3, 1), (0.25, 1), (0.1, 3), (1 / 16, 3), (0.12499999999999999, 2),
    ])
    def test_from_requirements_precision(self, cls, phase_error, l):
        """l = precision_bits(phase_error) - 1, so delta = 2^-(l+1) <= phase_error."""
        prep = cls.from_requirements(phase_error, 1e-4)
        assert prep.num_precision_qubits == l == precision_bits(phase_error) - 1
        assert prep.num_confidence_qubits == cls.confidence_qubits(l, 1e-4)

    def test_from_requirements_beats_rectangular(self):
        """Module docstring example: delta = 1/16 with p_fail 1e-8."""
        assert [cls.from_requirements(1 / 16, 1e-8).m_bits for cls in TAPERED] == [8]
        assert RectangularWindowState.from_requirements(1 / 16, 1e-8).m_bits == 4 + 26

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("phase_error, p_fail", [(0.0, 0.1), (1.0, 0.1), (0.1, 0.0),
                                                     (0.1, 1.0), (0.1, -0.5)])
    def test_from_requirements_rejects_out_of_range(self, cls, phase_error, p_fail):
        with pytest.raises(ValueError, match="phase_error|probability of failure"):
            cls.from_requirements(phase_error, p_fail)

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    def test_invalid_sizes_raise(self, cls):
        with pytest.raises(ValueError, match="num_confidence_qubits must be >= 1"):
            cls(3, 0)
        with pytest.raises(ValueError, match="num_precision_qubits must be >= 0"):
            cls(-1, 3)


class TestFromQubitCounts:

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("counts", [
        dict(num_precision_qubits=2, num_confidence_qubits=3),
        dict(num_phase_qubits=5, num_confidence_qubits=3),
        dict(num_phase_qubits=5, num_precision_qubits=2),
        dict(num_phase_qubits=5, num_precision_qubits=2, num_confidence_qubits=3),
        dict(num_phase_qubits=np.int64(5), num_precision_qubits=np.int64(2)),
    ], ids=["l, m", "p, m", "p, l", "p, l, m", "numpy ints"])
    def test_tapered_combinations(self, cls, counts):
        assert cls.from_qubit_counts(**counts) == cls(2, 3)

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    def test_tapered_rejects_num_phase_qubits_alone(self, cls):
        with pytest.raises(ValueError, match="cannot be sized from num_phase_qubits alone"):
            cls.from_qubit_counts(num_phase_qubits=5)

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("counts, match", [
        (dict(), "no qubit counts"),
        (dict(num_precision_qubits=2), r"Cannot size .* \['num_precision_qubits'\]"),
        (dict(num_confidence_qubits=3), r"Cannot size .* \['num_confidence_qubits'\]"),
        (dict(num_phase_qubits=5, num_precision_qubits=3, num_confidence_qubits=1),
         r"num_phase_qubits \(5\) must equal num_precision_qubits \(3\) \+ "
         r"num_confidence_qubits \(1\)"),
        (dict(num_phase_qubits=3, num_confidence_qubits=4),
         r"num_confidence_qubits \(4\) exceeds num_phase_qubits \(3\)"),
        (dict(num_phase_qubits=3, num_precision_qubits=4),
         r"num_precision_qubits \(4\) exceeds num_phase_qubits \(3\)"),
        (dict(num_phase_qubits=-1), "num_phase_qubits must be a non-negative integer"),
        (dict(num_precision_qubits=2.0, num_confidence_qubits=3),
         "num_precision_qubits must be a non-negative integer, got 2.0"),
        (dict(num_precision_qubits=2, num_confidence_qubits=True),
         "num_confidence_qubits must be a non-negative integer, got True"),
    ], ids=["none", "l", "m", "inconsistent", "m > p", "l > p", "negative", "float", "bool"])
    def test_invalid_counts_raise(self, cls, counts, match):
        with pytest.raises(ValueError, match=match):
            cls.from_qubit_counts(**counts)

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    def test_no_confidence_qubits_left_raises(self, cls):
        with pytest.raises(ValueError, match="num_confidence_qubits must be >= 1, got 0"):
            cls.from_qubit_counts(num_phase_qubits=4, num_precision_qubits=4)

    def test_flexible_qpe_from_qubit_counts(self):
        u = ZPowGate(exponent=0.47)
        qpe = FlexibleQPE.from_qubit_counts(u, num_precision_qubits=3, num_confidence_qubits=2,
                                            ancilla_prep=DPSSWindowState)
        assert qpe == FlexibleQPE(u, 5, DPSSWindowState(3, 2))

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    def test_flexible_qpe_from_requirements(self, cls):
        qpe = FlexibleQPE.from_requirements(ZPowGate(exponent=0.47), 1 / 16, 1e-8,
                                            ancilla_prep=cls)
        assert qpe.ancilla_prep == cls.from_requirements(1 / 16, 1e-8)
        assert qpe.num_ancilla_qubits == qpe.ancilla_prep.m_bits


# -------------------------------------------------------------------------------------------------
# Resources
# -------------------------------------------------------------------------------------------------

def gates_by_type(prep):
    """Gate counts of the decomposition, with ExactStatePreparation expanded (no Split/Join)."""
    cbloq = prep.decompose_bloq().flatten(lambda b: isinstance(b.bloq, ExactStatePreparation))
    return Counter(type(b.bloq).__name__ for b in cbloq.bloq_instances
                   if isinstance(b.bloq, (Ry, CNOT, ZPowGate, QFTTextBook)))


def call_graph_by_type(prep):
    counts = Counter()
    for callee, count in prep.build_call_graph(SympySymbolAllocator()):
        counts[type(callee).__name__] += count
    return counts


class TestResources:

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("l, m", [(0, 1), (2, 1), (1, 2), (2, 3), (1, 4)])
    def test_call_graph_matches_decomposition(self, cls, l, m):
        assert call_graph_by_type(cls(l, m)) == gates_by_type(cls(l, m))

    @pytest.mark.parametrize("l, m", [(2, 3), (1, 4)])
    def test_dpss_cost_is_generic_state_preparation(self, l, m):
        """Symmetry of the taper zeroes some angles, but every rotation is still counted."""
        p = l + m
        assert t_complexity(DPSSWindowState(l, m)) == TComplexity(rotations=2**p - 1,
                                                                  clifford=2**p - 2)

    def test_large_windows_are_costed_without_computing_taper(self):
        prep = DPSSWindowState(40, 4)
        with mock.patch.object(DPSSWindowState, "taper",
                               side_effect=AssertionError("taper computed")):
            assert get_cost_value(prep, QubitCount()) == 44
            assert t_complexity(prep).rotations == 2**44 - 1

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    def test_flexible_qpe_qubit_count(self, cls):
        qpe = FlexibleQPE(ZPowGate(exponent=0.47), 5, cls(2, 3))
        assert get_cost_value(qpe, QubitCount()) == 5 + 1
