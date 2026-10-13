"""
Tests for the tapered (DPSS-family) QPE window states.

Two definitions of failure appear here:
- QHAT's (used for sizing): the estimate k/N is more than delta = 2^-(num_precision_qubits + 1)
  from the phase.
- The paper's (eq. 15, used by its theorems): the estimate is not one of the 2K + 1 closest to the
  phase. Those 2K + 1 are all within delta, so QHAT's failure is at most the paper's.
"""

from collections import Counter
from unittest import mock

import numpy as np
import pytest
from scipy.signal.windows import dpss

from qualtran.bloqs.basic_gates import CNOT, Ry, ZPowGate
from qualtran.bloqs.qft.qft_text_book import QFTTextBook
from qualtran.cirq_interop.t_complexity_protocol import t_complexity, TComplexity
from qualtran.resource_counting import get_cost_value, QubitCount, SympySymbolAllocator

from qhat.common.exact_state_preparation import ExactStatePreparation
from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.qpe_window_state import precision_bits, RectangularWindowState
from qhat.common.tapered_window_state import (
    BandLimitedWindowStateBase,
    BandwidthOptimizedWindowState,
    dpss_average_failure_bound,
    DPSSWindowState,
    half_bandwidth,
    TruncatedDPSSWindowState,
)

TAPERED = [DPSSWindowState, TruncatedDPSSWindowState, BandwidthOptimizedWindowState]
BAND_LIMITED = [TruncatedDPSSWindowState, BandwidthOptimizedWindowState]
TAPERED_IDS = ["dpss", "truncated", "bandwidth-optimized"]
BAND_LIMITED_IDS = ["truncated", "bandwidth-optimized"]


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
    @pytest.mark.parametrize("num_precision, num_confidence",
                             [(0, 1), (2, 2), (1, 3), (3, 3), (2, 4)])
    def test_real_unit_norm_symmetric_positive(self, cls, num_precision, num_confidence):
        w = cls(num_precision, num_confidence).taper()
        assert w.shape == (2**(num_precision + num_confidence),) and np.isrealobj(w)
        assert np.linalg.norm(w) == pytest.approx(1.0)
        np.testing.assert_allclose(w, w[::-1], atol=1e-14)
        assert w.sum() > 0

    @pytest.mark.parametrize("num_precision, num_confidence", [(0, 2), (2, 2), (1, 3), (3, 3)])
    def test_dpss_is_top_eigenvector_of_sinc_kernel(self, num_precision, num_confidence):
        """The DPSS maximizes the energy in |f| <= (2K + 1)/(2N) (paper Sec. IV C)."""
        n, k_max = 2**(num_precision + num_confidence), half_bandwidth(num_confidence)
        _, vectors = np.linalg.eigh(sinc_kernel(n, k_max))
        top = vectors[:, -1] * np.sign(vectors[:, -1].sum())
        np.testing.assert_allclose(DPSSWindowState(num_precision, num_confidence).taper(), top,
                                   atol=1e-8)

    @pytest.mark.parametrize("num_precision, num_confidence",
                             [(0, 2), (2, 2), (1, 3), (3, 3), (2, 4)])
    def test_truncated_is_dpss_with_high_frequencies_removed(self, num_precision, num_confidence):
        """Paper eq. 92: keep DFT coefficients -K..K of the DPSS, renormalize."""
        n, k_max = 2**(num_precision + num_confidence), half_bandwidth(num_confidence)
        spectrum = np.fft.fft(dpss(n, k_max + 0.5, norm=2))
        frequency = np.fft.fftfreq(n, 1 / n)
        expected = np.fft.ifft(np.where(np.abs(frequency) <= k_max, spectrum, 0))
        expected /= np.linalg.norm(expected) * np.sign(expected.real.sum())
        np.testing.assert_allclose(TruncatedDPSSWindowState(num_precision, num_confidence).taper(),
                                   expected, atol=1e-14)

    @pytest.mark.parametrize("cls", BAND_LIMITED, ids=BAND_LIMITED_IDS)
    @pytest.mark.parametrize("num_precision, num_confidence", [(0, 1), (2, 2), (1, 3), (3, 4)])
    def test_band_limited_tapers_have_only_low_frequencies(self, cls, num_precision,
                                                           num_confidence):
        w = cls(num_precision, num_confidence).taper()
        frequency = np.fft.fftfreq(len(w), 1 / len(w))
        high = np.abs(frequency) > half_bandwidth(num_confidence)
        np.testing.assert_allclose(np.fft.fft(w)[high], 0, atol=1e-12)

    @pytest.mark.parametrize("cls", BAND_LIMITED, ids=BAND_LIMITED_IDS)
    def test_one_confidence_qubit_is_rectangular(self, cls):
        """One confidence qubit keeps only frequency 0 (K = 0)."""
        np.testing.assert_allclose(cls(3, 1).taper(), np.full(16, 0.25), atol=1e-15)

    @pytest.mark.parametrize("num_precision", [0, 1, 2, 3])
    def test_bandwidth_optimized_matches_eigenvector_reference(self, num_precision):
        """Sec. VII B's eigenproblem: maximize in-band energy over the 2K + 1 kept frequencies.
        (An eigensolver resolves 1 - energy only to ~1e-16, so this uses 2 confidence qubits.)"""
        num_confidence = 2
        n, k_max = 2**(num_precision + num_confidence), half_bandwidth(num_confidence)
        frequencies = np.r_[0:k_max + 1, n - k_max:n]
        basis = np.exp(2j * np.pi * np.outer(np.arange(n), frequencies) / n) / np.sqrt(n)
        form = basis.conj().T @ sinc_kernel(n, k_max) @ basis
        reference = basis @ np.linalg.eigh(form)[1][:, -1]
        reference *= np.exp(-1j * np.angle(reference.sum()))
        prep = BandwidthOptimizedWindowState(num_precision, num_confidence)
        np.testing.assert_allclose(prep.taper(), reference, atol=1e-12)

    @pytest.mark.parametrize("num_precision, num_confidence", [(0, 2), (2, 2), (1, 3), (3, 3)])
    def test_average_failure_ordering(self, num_precision, num_confidence):
        """DPSS <= bandwidth-optimized <= truncated: each is optimal over a larger set."""
        dpss_, trunc, bwopt = (average_paper_failure(cls(num_precision, num_confidence))
                               for cls in TAPERED)
        tol = 1e-14  # round-off in 1 - w C w
        assert dpss_ <= bwopt + tol and bwopt <= trunc + tol
        assert trunc <= 4 * dpss_  # eq. 110 (here eps(0) ~ average)

    def test_average_failure_is_out_of_band_energy(self):
        """The paper's failure, averaged over phases, is 1 - w^T C w (what the tapers optimize)."""
        prep = TruncatedDPSSWindowState(2, 2)
        w, n = prep.taper(), 16
        # Midpoint rule; the integrand jumps only at theta = 1/(2N), a cell boundary here.
        theta = (np.arange(2000) + 0.5) / (2000 * n)
        average = np.mean([paper_failure(w, t, prep.half_bandwidth) for t in theta])
        assert average == pytest.approx(average_paper_failure(prep), rel=1e-6)

    def test_frequency_amplitudes_normalization(self):
        for cls in BAND_LIMITED:
            r = cls(2, 4).frequency_amplitudes()
            assert len(r) == half_bandwidth(4) + 1 and r[0] > 0
            assert r[0]**2 + 2 * np.sum(r[1:]**2) == pytest.approx(1.0)


# -------------------------------------------------------------------------------------------------
# Circuits
# -------------------------------------------------------------------------------------------------

class TestCircuits:

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("num_precision, num_confidence",
                             [(0, 1), (2, 1), (0, 2), (1, 2), (2, 3), (3, 3), (1, 4)])
    def test_circuit_prepares_taper(self, cls, num_precision, num_confidence):
        """Including the phase (the taper is real and positive-sum, so no global phase)."""
        prep = cls(num_precision, num_confidence)
        np.testing.assert_allclose(window(prep), prep.taper(), atol=1e-13)

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    def test_signature(self, cls):
        prep = cls(2, 3)
        (reg,) = prep.signature
        assert reg.name == 'qpe_reg' and reg.total_bits() == 5
        assert prep.m_bits == 5 and prep.half_bandwidth == 3

    def test_frequency_state_layout(self):
        """r[|k|] at index k mod 2(K + 1); index K + 1 is empty (here K = 3)."""
        prep = TruncatedDPSSWindowState(2, 3)
        r = prep.frequency_amplitudes()
        np.testing.assert_array_equal(prep.frequency_state(),
                                      [r[0], r[1], r[2], r[3], 0, r[3], r[2], r[1]])
        assert np.linalg.norm(prep.frequency_state()) == pytest.approx(1.0)

    def test_phase_ramp_example(self):
        assert TruncatedDPSSWindowState(1, 2).phase_ramp_exponents() == (-0.25, -0.875)

    @pytest.mark.parametrize("cls", BAND_LIMITED, ids=BAND_LIMITED_IDS)
    @pytest.mark.parametrize("num_precision, num_confidence", [(0, 1), (1, 2), (3, 3), (2, 4)])
    def test_phase_ramp_centres_frequencies(self, cls, num_precision, num_confidence):
        """The ZPowGates give each value of the confidence qubits, read as a two's-complement
        integer k, the phase e^{-i pi k (N-1)/N}."""
        exponents = cls(num_precision, num_confidence).phase_ramp_exponents()
        assert len(exponents) == num_confidence and all(-1 < t <= 1 for t in exponents)
        n = 2**(num_precision + num_confidence)
        for value in range(2**num_confidence):
            # Most significant bit first, as in phase_ramp_exponents.
            bits = [(value >> (num_confidence - 1 - j)) & 1 for j in range(num_confidence)]
            k = value - 2**num_confidence if bits[0] else value
            phase = np.exp(1j * np.pi * np.dot(exponents, bits))
            assert phase == pytest.approx(np.exp(-1j * np.pi * k * (n - 1) / n))

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("theta", [0.3141, 0.8])
    def test_flexible_qpe_outcomes_match_window(self, cls, theta):
        """In FlexibleQPE, the outcome distribution is that of the window (eq. 13)."""
        prep = cls(2, 3)
        mat = FlexibleQPE(ZPowGate(exponent=2 * theta), 5, prep).tensor_contract()
        probs = np.abs(mat[:, 1].reshape(32, 2)[:, 1])**2  # qpe_reg = |0>, target = |1>
        np.testing.assert_allclose(probs, outcome_probs(prep.taper(), theta), atol=1e-13)

    def test_flexible_qpe_tensor_contract_matches_decomposition(self):
        qpe = FlexibleQPE(ZPowGate(exponent=0.47), 3, TruncatedDPSSWindowState(1, 2))
        np.testing.assert_allclose(qpe.tensor_contract(),
                                   qpe.decompose_bloq().tensor_contract(), atol=1e-12)


# -------------------------------------------------------------------------------------------------
# Failure probabilities
# -------------------------------------------------------------------------------------------------

class TestFailureProbability:

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("phase_error, p_fail", [(0.3, 0.1), (0.1, 1e-3), (0.07, 1e-8)])
    def test_from_requirements_meets_requirements_in_worst_case(self, cls, phase_error, p_fail):
        """Over phases across a grid spacing, Pr[|estimate - phase| > delta] <= p_fail, where
        delta = 2^-(num_precision_qubits + 1) <= phase_error (so the requirement holds too)."""
        prep = cls.from_requirements(phase_error, p_fail)
        delta = 2**-(prep.num_precision_qubits + 1)
        assert delta <= phase_error
        w = window(prep)
        worst = max(qhat_failure(w, t, delta) for t in thetas(len(w)))
        assert worst <= p_fail

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("num_precision, num_confidence",
                             [(0, 3), (1, 3), (3, 3), (1, 4), (3, 4)])
    def test_worst_case_meets_sizing_bound(self, cls, num_precision, num_confidence):
        """Sizing uses factor * Thm 1, a bound on the average; check it holds for the worst case."""
        w = window(cls(num_precision, num_confidence))
        delta = 2**-(num_precision + 1)
        worst = max(qhat_failure(w, t, delta) for t in thetas(len(w)))
        bound = dpss_average_failure_bound(num_precision, num_confidence)
        assert worst <= cls.failure_bound_factor * bound

    @pytest.mark.parametrize("num_precision, num_confidence", [(0, 2), (2, 2), (1, 3), (3, 3)])
    def test_qhat_failure_at_most_paper_failure(self, num_precision, num_confidence):
        w = DPSSWindowState(num_precision, num_confidence).taper()
        delta, k_max = 2**-(num_precision + 1), half_bandwidth(num_confidence)
        for t in thetas(len(w), 9):
            assert qhat_failure(w, t, delta) <= paper_failure(w, t, k_max) + 1e-15

    @pytest.mark.parametrize("cls", BAND_LIMITED, ids=BAND_LIMITED_IDS)
    @pytest.mark.parametrize("num_precision, num_confidence",
                             [(0, 2), (2, 2), (4, 2), (1, 3), (3, 3)])
    def test_truncation_bound_pointwise(self, cls, num_precision, num_confidence):
        """Eq. 107 for the truncated DPSS: failure(theta) <= (sqrt(eps(0)) + sqrt(eps(theta)))^2,
        eps being the DPSS's failure. The bandwidth-optimized taper is not covered by eq. 107;
        here it meets the same bound."""
        k_max = half_bandwidth(num_confidence)
        w = window(cls(num_precision, num_confidence))
        w_dpss = DPSSWindowState(num_precision, num_confidence).taper()
        eps0 = paper_failure(w_dpss, 0.0, k_max)
        for t in thetas(len(w), 17):
            bound = (np.sqrt(eps0) + np.sqrt(paper_failure(w_dpss, t, k_max)))**2
            assert paper_failure(w, t, k_max) <= bound * (1 + 1e-9) + 1e-15

    @pytest.mark.parametrize("num_precision, num_confidence",
                             [(0, 3), (1, 3), (3, 3), (1, 4), (6, 4)])
    def test_theorem_1_bounds_dpss_average(self, num_precision, num_confidence):
        average = average_paper_failure(DPSSWindowState(num_precision, num_confidence))
        assert average <= dpss_average_failure_bound(num_precision, num_confidence)


# -------------------------------------------------------------------------------------------------
# Sizing
# -------------------------------------------------------------------------------------------------

class TestSizing:

    @pytest.mark.parametrize("num_confidence, k_max", [(1, 0), (2, 1), (3, 3), (4, 7), (6, 31)])
    def test_half_bandwidth(self, num_confidence, k_max):
        assert half_bandwidth(num_confidence) == k_max

    def test_theorem_1_example(self):
        # min(8 exp(-29 / ((2/pi^2) ln 1024)), 10 exp(-24 / ((2/pi^2) ln 1575))) = 8.64e-9
        assert dpss_average_failure_bound(3, 5) == pytest.approx(8.640019652e-9)

    def test_theorem_1_decreases_with_confidence_qubits(self):
        bounds = [dpss_average_failure_bound(3, num_confidence) for num_confidence in range(1, 9)]
        assert all(a > b for a, b in zip(bounds, bounds[1:]))

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("num_precision", [0, 3, 10])
    @pytest.mark.parametrize("p_fail", [0.3, 1e-3, 1e-12])
    def test_confidence_qubits_is_smallest_meeting_bound(self, cls, num_precision, p_fail):
        num_confidence = cls.confidence_qubits(num_precision, p_fail)
        factor = cls.failure_bound_factor
        assert factor * dpss_average_failure_bound(num_precision, num_confidence) <= p_fail
        assert (num_confidence == 1
                or factor * dpss_average_failure_bound(num_precision, num_confidence - 1) > p_fail)

    def test_confidence_qubits_example(self):
        assert [cls.confidence_qubits(3, 1e-8) for cls in TAPERED] == [5, 6, 6]
        assert [cls.failure_bound_factor for cls in TAPERED] == [1, 4, 4]

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("phase_error, num_precision", [
        (0.75, 0), (0.3, 1), (0.25, 1), (0.1, 3), (1 / 16, 3), (0.12499999999999999, 2),
    ])
    def test_from_requirements_precision(self, cls, phase_error, num_precision):
        """precision_bits(phase_error) - 1 precision qubits, so
        delta = 2^-(num_precision_qubits + 1) <= phase_error."""
        prep = cls.from_requirements(phase_error, 1e-4)
        assert prep.num_precision_qubits == num_precision == precision_bits(phase_error) - 1
        assert prep.num_confidence_qubits == cls.confidence_qubits(num_precision, 1e-4)

    def test_from_requirements_beats_rectangular(self):
        """Module docstring example: delta = 1/16 with p_fail 1e-8."""
        assert [cls.from_requirements(1 / 16, 1e-8).m_bits for cls in TAPERED] == [8, 9, 9]
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
    ], ids=["precision, confidence", "phase, confidence", "phase, precision", "all three",
            "numpy ints"])
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
    ], ids=["none", "precision only", "confidence only", "inconsistent", "confidence > phase",
            "precision > phase", "negative", "float", "bool"])
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

    def test_flexible_qpe_from_qubit_counts_with_qft_inv(self):
        u = ZPowGate(exponent=0.47)
        qft_inv = lambda bitsize: QFTTextBook(bitsize, with_reverse=False).adjoint()  # noqa: E731
        qpe = FlexibleQPE.from_qubit_counts(u, num_phase_qubits=5, num_confidence_qubits=2,
                                            ancilla_prep=TruncatedDPSSWindowState,
                                            qft_inv=qft_inv)
        assert qpe.ancilla_prep == TruncatedDPSSWindowState(3, 2)
        assert qpe.qft_inv == qft_inv(5)

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
    @pytest.mark.parametrize("num_precision, num_confidence",
                             [(0, 1), (2, 1), (1, 2), (2, 3), (1, 4)])
    def test_call_graph_matches_decomposition(self, cls, num_precision, num_confidence):
        prep = cls(num_precision, num_confidence)
        assert call_graph_by_type(prep) == gates_by_type(prep)

    @pytest.mark.parametrize("num_precision, num_confidence", [(2, 3), (1, 4)])
    def test_dpss_cost_is_generic_state_preparation(self, num_precision, num_confidence):
        """Symmetry of the taper zeroes some angles, but every rotation is still counted."""
        n = 2**(num_precision + num_confidence)
        assert t_complexity(DPSSWindowState(num_precision, num_confidence)) == TComplexity(
            rotations=n - 1, clifford=n - 2)

    @pytest.mark.parametrize("cls", BAND_LIMITED, ids=BAND_LIMITED_IDS)
    @pytest.mark.parametrize("num_precision, num_confidence",
                             [(0, 2), (1, 2), (2, 3), (3, 3), (1, 4)])
    def test_band_limited_t_complexity_is_sum_over_gates(self, cls, num_precision, num_confidence):
        """The kept-frequency angles are generic here, so gate-by-gate costs match exactly."""
        prep = cls(num_precision, num_confidence)
        cbloq = prep.decompose_bloq().flatten(
            lambda b: isinstance(b.bloq, ExactStatePreparation))
        total = TComplexity()
        for binst in cbloq.bloq_instances:
            if isinstance(binst.bloq, (Ry, CNOT, ZPowGate, QFTTextBook)):
                total += t_complexity(binst.bloq)
        assert t_complexity(prep) == total

    def test_band_limited_cost_example(self):
        """(3, 3): 7 Ry + 3 ramp rotations + 15 in the 6-qubit QFT; 6 + 3 CNOTs + 9 in the QFT."""
        assert t_complexity(QFTTextBook(6, with_reverse=True)) == TComplexity(rotations=15,
                                                                            clifford=9)
        for cls in BAND_LIMITED:
            assert t_complexity(cls(3, 3)) == TComplexity(rotations=25, clifford=18)

    @pytest.mark.parametrize("cls, method", [
        (DPSSWindowState, "taper"),
        (TruncatedDPSSWindowState, "frequency_amplitudes"),
        (BandwidthOptimizedWindowState, "frequency_amplitudes"),
    ], ids=TAPERED_IDS)
    def test_large_windows_are_costed_without_computing_taper(self, cls, method):
        prep = cls(40, 4)
        with mock.patch.object(cls, method, side_effect=AssertionError("taper computed")):
            assert get_cost_value(prep, QubitCount()) == 44
            cost = t_complexity(prep)
        if cls is DPSSWindowState:
            assert cost.rotations == 2**44 - 1
        else:
            # 15 Ry + 4 ramp rotations (none a Clifford or T here) + the 44-qubit QFT.
            qft = t_complexity(QFTTextBook(44, with_reverse=True))
            assert cost.rotations == 15 + 4 + qft.rotations

    @pytest.mark.parametrize("cls", TAPERED, ids=TAPERED_IDS)
    def test_flexible_qpe_qubit_count(self, cls):
        qpe = FlexibleQPE(ZPowGate(exponent=0.47), 5, cls(2, 3))
        assert get_cost_value(qpe, QubitCount()) == 5 + 1

    def test_band_limited_base_is_abstract(self):
        with pytest.raises(TypeError):
            BandLimitedWindowStateBase(2, 3)
