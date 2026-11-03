"""
Tests for the multi-run sizing of tapered QPE (`qhat.common.tapered_qpe_sizing`).

Reference values for the default problem come from the reference implementation of the
multi-run analysis by Yigit Subasi and Alan Bu. At very large N that implementation loses digits;
there the values are checked against 60-digit mpmath instead.
"""

import itertools
import math

import mpmath
import numpy as np
import pytest
from scipy.stats import binom

from qhat.common.qpe_window_state import precision_bits
from qhat.common.tapered_qpe_sizing import (
    COST_CEILING,
    fewest_runs,
    MAX_CONFIDENCE_QUBITS,
    MIN_CONFIDENCE_QUBITS,
    multi_run_failure,
    per_run_qpe_failure,
    size_tapered_qpe,
    split_phase_qubits,
)
from qhat.common.tapered_window_state import (
    BandwidthOptimizedWindowState,
    DPSSWindowState,
    TruncatedDPSSWindowState,
)

TAPERED = [DPSSWindowState, TruncatedDPSSWindowState, BandwidthOptimizedWindowState]
TAPERED_IDS = ["dpss", "truncated", "bandwidth-optimized"]

# Energy error 0.2 on a spectrum scaled by 1.7: 3 precision qubits.
PHASE_ERROR = 0.2 / 1.7
# Phase error for 10 precision qubits, as qhat's driver sets it (phase scale factor 1.01).
FINE_PHASE_ERROR = 1 / (1.01 * 2**10)


def direct_failure(n, k, under, over):
    """P(the k-th smallest of n runs is an under- or overestimate), summed term by term."""
    def at_least(j_min, q):
        return sum(math.comb(n, j) * q**j * (1 - q)**(n - j) for j in range(j_min, n + 1))
    # Underestimate: at least k runs underestimate. Overestimate: at least n - k + 1 overestimate.
    return at_least(k, under) + at_least(n - k + 1, over)


def mpmath_failure(window_cls, num_precision, num_confidence, n, k, overlap):
    """`multi_run_failure` in 60-digit arithmetic."""
    with mpmath.workdps(60):
        u = mpmath.mpf(per_run_qpe_failure(window_cls, num_precision, num_confidence))
        hit = mpmath.mpf(overlap) * (1 - u)

        def failure(under, not_over):
            return (mpmath.betainc(k, n - k + 1, 0, under, regularized=True)
                    + mpmath.betainc(k, n - k + 1, not_over, 1, regularized=True))
        return max(failure(0, hit), failure(u, u + hit))


class TestMultiRunFailure:
    @pytest.mark.parametrize("n, k", [(1, 1), (5, 2), (14, 2), (14, 9), (30, 30)])
    @pytest.mark.parametrize("overlap", [0.05, 0.5, 1.0])
    def test_matches_direct_sum(self, n, k, overlap):
        u = per_run_qpe_failure(DPSSWindowState, 1, 3)
        hit = overlap * (1 - u)
        expected = max(direct_failure(n, k, 0.0, 1 - hit),
                       direct_failure(n, k, u, 1 - u - hit))
        assert multi_run_failure(DPSSWindowState, 1, 3, n, k, overlap) == pytest.approx(
            expected, rel=1e-10, abs=1e-15)

    @pytest.mark.parametrize("n, k", [(7, 1), (14, 2), (25, 6)])
    def test_worst_split_is_an_endpoint(self, n, k):
        # The adversary may split the single-run failure 1 - hit into underestimates (at most u)
        # and overestimates however it likes; no split does worse than the two extremes.
        u = per_run_qpe_failure(DPSSWindowState, 0, 3)
        overlap = 0.4
        hit = overlap * (1 - u)
        bound = multi_run_failure(DPSSWindowState, 0, 3, n, k, overlap)
        for under in np.linspace(0, u, 21):
            assert direct_failure(n, k, under, 1 - hit - under) <= bound * (1 + 1e-12)

    def test_single_run_at_overlap_one_is_the_qpe_failure(self):
        u = per_run_qpe_failure(DPSSWindowState, 3, 4)
        assert multi_run_failure(DPSSWindowState, 3, 4, 1, 1, 1.0) == pytest.approx(u, rel=1e-12)

    def test_overlap_table(self):
        # The plan sized for overlap 0.5 (3 + 4 qubits, 2nd-smallest of 14 runs), if the true
        # overlap is different.
        overlaps = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        expected = [9.92e-1, 5.85e-1, 1.98e-1, 4.75e-2, 8.11e-3, 9.19e-4, 6.56e-5, 8.40e-6,
                    6.80e-6, 6.79e-6, 6.79e-6]
        got = [multi_run_failure(DPSSWindowState, 3, 4, 14, 2, p) for p in overlaps]
        assert got == pytest.approx(expected, rel=5e-3)

    @pytest.mark.parametrize("n, k, overlap", [(1383365, 1, 1e-5), (239448133, 2, 1e-7),
                                               (3367460167391, 6, 1e-11)])
    def test_large_n_matches_high_precision(self, n, k, overlap):
        expected = float(mpmath_failure(DPSSWindowState, 10, 6, n, k, overlap))
        assert multi_run_failure(DPSSWindowState, 10, 6, n, k, overlap) == pytest.approx(
            expected, rel=1e-9)


class TestSizeTaperedQPE:
    @pytest.mark.parametrize("overlap, runs, quantile, failure", [
        (1.0, 1, 1, 2.735e-4),
        (0.9, 5, 2, 4.644e-4),
        (0.5, 14, 2, 9.193e-4),
        (0.25, 33, 2, 9.361e-4),
        (0.1, 92, 2, 9.850e-4),
    ])
    def test_reference_values(self, overlap, runs, quantile, failure):
        plan = size_tapered_qpe(DPSSWindowState, PHASE_ERROR, 1e-3, overlap)
        assert (plan.num_precision_qubits, plan.num_confidence_qubits) == (3, 4)
        assert (plan.num_runs, plan.quantile) == (runs, quantile)
        assert plan.multi_run_failure == pytest.approx(failure, rel=1e-3)
        assert plan.single_run_failure == pytest.approx(
            1 - overlap * (1 - per_run_qpe_failure(DPSSWindowState, 3, 4)), rel=1e-12)
        assert plan.overlap == overlap
        assert plan.force_lowest_estimate is False

    def test_reference_preset_best_and_lowest_estimate(self):
        # Energy error 2^-10 on a unit spectrum padded by the error on both sides; 95% confidence.
        delta = 2**-10
        phase_error = delta / (1 + 2 * delta)
        best = size_tapered_qpe(DPSSWindowState, phase_error, 0.05, 0.5)
        assert (best.num_runs, best.quantile) == (8, 2)
        assert (best.num_precision_qubits, best.num_confidence_qubits) == (10, 4)
        assert 1 - best.multi_run_failure == pytest.approx(0.954553, abs=1e-6)
        lowest = size_tapered_qpe(DPSSWindowState, phase_error, 0.05, 0.5,
                                  force_lowest_estimate=True)
        assert (lowest.num_runs, lowest.quantile) == (5, 1)
        assert (lowest.num_precision_qubits, lowest.num_confidence_qubits) == (10, 5)
        assert lowest.force_lowest_estimate is True

    @pytest.mark.parametrize("window_cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("phase_error, probability_of_failure",
                             [(PHASE_ERROR, 1e-3), (0.03, 1e-8), (FINE_PHASE_ERROR, 0.2)])
    def test_overlap_one_matches_single_run_rule(self, window_cls, phase_error,
                                                 probability_of_failure):
        plan = size_tapered_qpe(window_cls, phase_error, probability_of_failure, 1.0)
        window = window_cls.from_requirements(phase_error, probability_of_failure)
        assert (plan.num_runs, plan.quantile) == (1, 1)
        assert plan.num_precision_qubits == window.num_precision_qubits
        assert plan.num_confidence_qubits == window.num_confidence_qubits

    def test_band_limited_needs_at_least_as_much(self):
        dpss = size_tapered_qpe(DPSSWindowState, PHASE_ERROR, 1e-6, 0.3)
        for window_cls in TAPERED[1:]:
            other = size_tapered_qpe(window_cls, PHASE_ERROR, 1e-6, 0.3)
            assert other.controlled_unitary_cost >= dpss.controlled_unitary_cost

    @pytest.mark.parametrize("overlap, probability_of_failure", [(0.5, 1e-3), (0.07, 1e-5)])
    @pytest.mark.parametrize("force", [False, True])
    def test_cheapest_over_confidence_qubits(self, overlap, probability_of_failure, force):
        plan = size_tapered_qpe(TruncatedDPSSWindowState, PHASE_ERROR, probability_of_failure,
                                overlap, force)
        cost = plan.controlled_unitary_cost
        for num_confidence in range(MIN_CONFIDENCE_QUBITS, MAX_CONFIDENCE_QUBITS + 1):
            other = fewest_runs(TruncatedDPSSWindowState, 3, num_confidence,
                                probability_of_failure, overlap, force)
            if other is not None:
                assert other.controlled_unitary_cost >= cost

    @pytest.mark.parametrize("overlap, probability_of_failure, runs, quantile", [
        (1e-5, 1e-6, 1383365, 1),
        (1e-5, 1e-9, 2393962, 2),
        (1e-7, 1e-9, 239448133, 2),
        (1e-11, 1e-9, 3367460167391, 6),
    ])
    def test_large_n(self, overlap, probability_of_failure, runs, quantile):
        plan = size_tapered_qpe(DPSSWindowState, FINE_PHASE_ERROR, probability_of_failure,
                                overlap)
        assert (plan.num_precision_qubits, plan.num_confidence_qubits) == (10, 6)
        assert (plan.num_runs, plan.quantile) == (runs, quantile)
        assert mpmath_failure(DPSSWindowState, 10, 6, runs, quantile,
                              overlap) <= probability_of_failure
        assert mpmath_failure(DPSSWindowState, 10, 6, runs - 1, quantile,
                              overlap) > probability_of_failure

    def test_infeasible_within_cost_ceiling(self):
        # Even K = 1 needs about 2e15 runs of at least 2^12 queries each.
        assert size_tapered_qpe(DPSSWindowState, FINE_PHASE_ERROR, 1e-9, 1e-14) is None
        assert size_tapered_qpe(DPSSWindowState, FINE_PHASE_ERROR, 1e-9, 1e-14, True) is None

    @pytest.mark.parametrize("overlap", [0.0, -0.1, 1.5])
    def test_rejects_overlap(self, overlap):
        with pytest.raises(ValueError, match="overlap"):
            size_tapered_qpe(DPSSWindowState, PHASE_ERROR, 1e-3, overlap)

    @pytest.mark.parametrize("probability_of_failure", [0.0, 1.0])
    def test_rejects_probability_of_failure(self, probability_of_failure):
        with pytest.raises(ValueError, match="probability of failure"):
            size_tapered_qpe(DPSSWindowState, PHASE_ERROR, probability_of_failure, 0.5)


class TestFewestRuns:
    @pytest.mark.parametrize("window_cls", TAPERED, ids=TAPERED_IDS)
    # Answers range from N = 6 to 57, K = 1 to 17; forcing K = 1 is impossible in some, and with 2
    # confidence qubits nothing works.
    @pytest.mark.parametrize("num_precision, num_confidence, overlap, probability_of_failure", [
        (1, 3, 0.6, 1e-2),
        (2, 4, 0.6, 1e-2),
        (0, 4, 0.8, 1e-4),
        (1, 4, 0.3, 1e-3),
        (2, 3, 0.5, 0.05),
        (3, 2, 0.9, 0.3),
    ])
    @pytest.mark.parametrize("force", [False, True])
    def test_matches_exhaustive_search(self, window_cls, num_precision, num_confidence, overlap,
                                       probability_of_failure, force):
        def failure(n, k):
            return multi_run_failure(window_cls, num_precision, num_confidence, n, k, overlap)

        expected = None
        for n in range(1, 120):
            quantiles = [1] if force else range(1, n + 1)
            if any(failure(n, k) <= probability_of_failure for k in quantiles):
                expected = (n, min(quantiles, key=lambda k: failure(n, k)))
                break
        plan = fewest_runs(window_cls, num_precision, num_confidence, probability_of_failure,
                           overlap, force)
        if expected is None:
            assert plan is None
            return
        assert plan.num_runs == expected[0]
        assert failure(plan.num_runs, plan.quantile) == pytest.approx(
            failure(*expected), rel=1e-9)
        assert plan.multi_run_failure <= probability_of_failure

    def test_no_runs_suffice(self):
        # One confidence qubit: the per-run failure bound is 1, so no run is ever trusted.
        assert fewest_runs(DPSSWindowState, 3, 1, 1e-3, 0.5) is None

    def test_respects_cost_ceiling(self):
        num_precision, num_confidence = 10, 6
        plan = fewest_runs(DPSSWindowState, num_precision, num_confidence, 1e-9, 1e-11)
        assert plan.controlled_unitary_cost <= COST_CEILING
        assert fewest_runs(DPSSWindowState, num_precision, num_confidence, 1e-9, 1e-14) is None


def outcome_probs(window, phase):
    """QPE outcome probabilities for a window state and a phase (turns): |QFT^-1 (w e^{2 pi i n
    phase})|^2."""
    n = np.arange(len(window))
    return np.abs(np.fft.fft(window * np.exp(2j * np.pi * n * phase)))**2 / len(window)


def simulated_multi_run_failure(window, phases, weights, phase_error, num_runs, quantile):
    """Exact failure of reporting the `quantile`-th smallest of `num_runs` runs, for an initial
    state with weights `weights[j]` on eigenphases `phases[j]`; phases[0] is the target (lowest).
    """
    num_outcomes = len(window)
    probs = sum(c * outcome_probs(window, phase) for c, phase in zip(weights, phases))
    estimates = np.arange(num_outcomes) / num_outcomes
    distance = np.abs(estimates - phases[0])
    accepted = np.minimum(distance, 1 - distance) <= phase_error
    under = ~accepted & (estimates < phases[0])
    p_under, p_over = probs[under].sum(), probs[~accepted & ~under].sum()
    # At least `quantile` underestimates, or at least num_runs - quantile + 1 overestimates.
    return (binom.sf(quantile - 1, num_runs, p_under)
            + binom.sf(num_runs - quantile, num_runs, p_over))


class TestModelAgainstQPE:
    """The multi-run bound against exact QPE outcome distributions for random spectra and initial
    states, with every eigenphase at least the phase error from the window edges."""

    @pytest.mark.parametrize("window_cls", TAPERED, ids=TAPERED_IDS)
    def test_bound_holds_for_random_spectra(self, window_cls):
        rng = np.random.default_rng(7)
        for _ in range(40):
            num_precision = int(rng.integers(1, 4))
            num_confidence = int(rng.integers(2, 5))
            overlap = float(rng.uniform(0.05, 1.0))
            probability_of_failure = 10**rng.uniform(-6, -1)
            plan = fewest_runs(window_cls, num_precision, num_confidence, probability_of_failure,
                               overlap)
            if plan is None:
                continue
            phase_error = 2.0**-(num_precision + 1)
            num_excited = int(rng.integers(0, 6))
            phases = np.sort(rng.uniform(phase_error, 1 - phase_error, num_excited + 1))
            weights = np.concatenate(
                [[overlap], rng.dirichlet(np.ones(num_excited)) * (1 - overlap)
                 if num_excited else []])
            if num_excited == 0:
                weights = np.array([1.0])
            window = window_cls(num_precision, num_confidence).taper()
            failure = simulated_multi_run_failure(window, phases, weights, phase_error,
                                                  plan.num_runs, plan.quantile)
            assert failure <= plan.multi_run_failure * (1 + 1e-9)

    def test_bound_needs_margin_below_phase_one(self):
        # An excited state within the phase error of phase 1: its correct estimates wrap to
        # near 0 and count as underestimates, so the bound fails by a factor of about 27. This
        # is the window-edge margin the driver does not yet guarantee.
        num_precision, num_confidence, overlap = 3, 4, 0.5
        phase_error = 2.0**-(num_precision + 1)
        plan = fewest_runs(DPSSWindowState, num_precision, num_confidence, 1e-3, overlap)
        window = DPSSWindowState(num_precision, num_confidence).taper()
        weights = [overlap, 1 - overlap]
        inside = simulated_multi_run_failure(window, [0.3, 0.6], weights, phase_error,
                                             plan.num_runs, plan.quantile)
        near_edge = simulated_multi_run_failure(window, [0.3, 1 - 0.3 * phase_error], weights,
                                                phase_error, plan.num_runs, plan.quantile)
        assert inside <= plan.multi_run_failure
        assert near_edge > 20 * plan.multi_run_failure


class TestSplitPhaseQubits:
    """The cheapest split of a given total: the fewest precision qubits, the rest confidence."""

    @pytest.mark.parametrize("phase_error, num_phase_qubits, split",
                             [(PHASE_ERROR, 9, (3, 6)), (PHASE_ERROR, 4, (3, 1)),
                              (FINE_PHASE_ERROR, 14, (10, 4)), (0.3, 5, (1, 4))])
    def test_values(self, phase_error, num_phase_qubits, split):
        assert split_phase_qubits(phase_error, num_phase_qubits) == split

    @pytest.mark.parametrize("window_cls", TAPERED, ids=TAPERED_IDS)
    @pytest.mark.parametrize("phase_error", [PHASE_ERROR, 0.03, FINE_PHASE_ERROR])
    @pytest.mark.parametrize("extra", [2, 4, 7])
    @pytest.mark.parametrize("overlap, probability_of_failure", [(1.0, 1e-6), (0.5, 1e-2),
                                                                 (0.1, 1e-6)])
    def test_cheapest_split(self, window_cls, phase_error, extra, overlap,
                            probability_of_failure):
        # Every split meeting the phase error costs the same per run, so the cheapest needs the
        # fewest runs; none needs fewer than the minimal-precision split.
        num_precision, num_confidence = split_phase_qubits(
            phase_error, precision_bits(phase_error) - 1 + extra)
        total = num_precision + num_confidence
        best = fewest_runs(window_cls, num_precision, num_confidence, probability_of_failure,
                           overlap)
        for more_precision in range(num_precision + 1, total):
            other = fewest_runs(window_cls, more_precision, total - more_precision,
                                probability_of_failure, overlap)
            if other is not None:
                assert best is not None and best.num_runs <= other.num_runs

    def test_too_few_qubits(self):
        # 0.2/1.7 turns needs 3 precision qubits, so 3 in total leaves no confidence qubit.
        with pytest.raises(ValueError, match="3 phase qubits cannot meet"):
            split_phase_qubits(PHASE_ERROR, 3)
