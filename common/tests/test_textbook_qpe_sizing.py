"""
Tests for the multi-run sizing of textbook QPE (`qhat.common.textbook_qpe_sizing`).

The per-run bound is checked against the exact outcome distribution of textbook QPE (summed over
all 2^m outcomes) and against mpmath; plans are checked against brute force.
"""

import math

import mpmath
import numpy as np
import pytest
from scipy.stats import binom

from qhat.common.multi_run_qpe_sizing import COST_CEILING, worst_failure
from qhat.common.qpe_window_state import precision_bits, RectangularWindowState
from qhat.common.textbook_qpe_sizing import (
    fewest_runs,
    multi_run_failure,
    per_run_qpe_failure,
    size_textbook_qpe,
)

# Energy error 0.2 on a spectrum scaled by 1.7.
PHASE_ERROR = 0.2 / 1.7
# Phase error as qhat's driver sets it for 10 bits (phase scale factor 1.01).
FINE_PHASE_ERROR = 1 / (1.01 * 2**10)


def exact_failure(phase_error, num_phase_qubits, offset):
    """P(estimate more than `phase_error` turns from the phase, circularly), for a phase `offset`
    grid steps above a grid point, summed over all 2^m outcomes of textbook QPE."""
    m = 2**num_phase_qubits
    d = np.arange(m)
    # distance in grid steps from the phase to each outcome, wrapped to [-m/2, m/2)
    x = (d - offset + m / 2) % m - m / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.sin(np.pi * offset)**2 / (m**2 * np.sin(np.pi * (d - offset) / m)**2)
    p = np.where(np.isclose(d, offset), 1.0, p)
    return float(p[np.abs(x) > phase_error * m].sum())


def worst_exact_failure(phase_error, num_phase_qubits, points=2001):
    """Largest `exact_failure` over offsets. The peaks lie where the acceptance window's edge
    crosses a grid point, so those offsets are included."""
    h = phase_error * 2**num_phase_qubits
    edges = [(h % 1) + s for s in (-1e-9, 1e-9)] + [(-h % 1) + s for s in (-1e-9, 1e-9)]
    offsets = np.concatenate([np.linspace(1e-6, 1 - 1e-6, points), np.clip(edges, 1e-6, 1 - 1e-6)])
    return max(exact_failure(phase_error, num_phase_qubits, f) for f in offsets)


class TestPerRunBound:
    @pytest.mark.parametrize("phase_error", [PHASE_ERROR, 0.07, 0.3, 1 / (1.01 * 2**5)])
    @pytest.mark.parametrize("num_phase_qubits", [2, 3, 5, 7, 9])
    def test_bounds_exact_failure_at_every_offset(self, phase_error, num_phase_qubits):
        u = per_run_qpe_failure(phase_error, num_phase_qubits)
        assert worst_exact_failure(phase_error, num_phase_qubits, points=101) <= u * (1 + 1e-12)

    # True ratios at these sizes: 1.123, 1.061 and 1.048.
    @pytest.mark.parametrize("num_phase_qubits, ratio", [(6, 1.13), (8, 1.07), (10, 1.06)])
    def test_close_to_exact_worst_case(self, num_phase_qubits, ratio):
        u = per_run_qpe_failure(PHASE_ERROR, num_phase_qubits)
        worst = worst_exact_failure(PHASE_ERROR, num_phase_qubits)
        assert worst <= u <= ratio * worst

    @pytest.mark.parametrize("h", [0.75, 1.0, 3.7, 30.1, 481.9, 1.5e6])
    def test_matches_high_precision(self, h):
        num_phase_qubits = 30
        phase_error = h / 2**num_phase_qubits
        with mpmath.workdps(40):
            expected = float(2 / mpmath.pi**2
                             * mpmath.psi(1, mpmath.mpf(phase_error) * 2**num_phase_qubits))
        u = per_run_qpe_failure(phase_error, num_phase_qubits)
        assert expected <= u <= expected * (1 + 1e-11)

    def test_tighter_than_nielsen_and_chuang(self):
        n = precision_bits(PHASE_ERROR)
        for confidence in range(3, 12):
            nielsen_chuang = 1 / (2 * (2**confidence - 2))
            assert per_run_qpe_failure(PHASE_ERROR, n + confidence) < nielsen_chuang / 2

    def test_capped_at_one(self):
        assert per_run_qpe_failure(0.2, 1) == 1.0

    @pytest.mark.parametrize("phase_error", [0.0, -0.1, 1.0, 1.3])
    def test_rejects_phase_error(self, phase_error):
        with pytest.raises(ValueError, match="phase error"):
            per_run_qpe_failure(phase_error, 5)

    @pytest.mark.parametrize("phase_error", [0.5, 0.7, 0.99])
    def test_large_phase_error(self, phase_error):
        # From W = 1/2 on nothing is rejected. The driver reaches W near 1 when energy_error is
        # the whole spectrum width; the old single-run rule accepted it, so the sizing must too.
        assert worst_exact_failure(phase_error, 3, points=101) == 0.0
        assert per_run_qpe_failure(phase_error, 3) < 0.06
        plan = size_textbook_qpe(phase_error, 1e-3, 0.5)
        assert plan.multi_run_failure <= 1e-3


class TestSizeTextbookQPE:
    @pytest.mark.parametrize("phase_error, probability_of_failure, overlap, qubits, runs, quantile",
                             [
                                 (PHASE_ERROR, 1e-3, 1.0, 4, 13, 7),
                                 (PHASE_ERROR, 1e-3, 0.5, 5, 38, 9),
                                 (PHASE_ERROR, 1e-3, 0.1, 6, 450, 26),
                                 (FINE_PHASE_ERROR, 1e-6, 1.0, 11, 29, 15),
                                 (FINE_PHASE_ERROR, 1e-6, 0.5, 12, 82, 18),
                                 (FINE_PHASE_ERROR, 1e-6, 0.1, 13, 983, 55),
                             ])
    def test_reference_values(self, phase_error, probability_of_failure, overlap, qubits, runs,
                              quantile):
        plan = size_textbook_qpe(phase_error, probability_of_failure, overlap)
        assert (plan.num_phase_qubits, plan.num_runs, plan.quantile) == (qubits, runs, quantile)
        u = per_run_qpe_failure(phase_error, qubits)
        assert plan.single_run_failure == pytest.approx(1 - overlap * (1 - u), rel=1e-12)
        assert plan.multi_run_failure <= probability_of_failure
        assert plan.overlap == overlap
        assert plan.force_lowest_estimate is False

    @pytest.mark.parametrize("phase_error, probability_of_failure, overlap",
                             [(PHASE_ERROR, 1e-3, 0.5), (FINE_PHASE_ERROR, 1e-6, 0.1),
                              (0.03, 1e-2, 0.8)])
    @pytest.mark.parametrize("force", [False, True])
    def test_cheapest_over_phase_qubits(self, phase_error, probability_of_failure, overlap, force):
        plan = size_textbook_qpe(phase_error, probability_of_failure, overlap, force)
        cost = plan.num_runs * 2**plan.num_phase_qubits
        for num_phase_qubits in range(1, 40):
            other = fewest_runs(phase_error, num_phase_qubits, probability_of_failure, overlap,
                                force)
            if other is not None:
                assert other.num_runs * 2**num_phase_qubits >= cost

    def test_fewest_runs_matches_brute_force(self):
        # The plans above are the fewest runs at their qubit counts.
        for overlap, qubits, runs in [(1.0, 4, 13), (0.5, 5, 38), (0.1, 6, 450)]:
            u = per_run_qpe_failure(PHASE_ERROR, qubits)
            hit = overlap * (1 - u)
            assert min(worst_failure(runs, k, u, hit) for k in range(1, runs + 1)) <= 1e-3
            assert min(worst_failure(runs - 1, k, u, hit) for k in range(1, runs)) > 1e-3

    def test_lowest_estimate_needs_more_qubits(self):
        best = size_textbook_qpe(FINE_PHASE_ERROR, 1e-6, 0.5)
        lowest = size_textbook_qpe(FINE_PHASE_ERROR, 1e-6, 0.5, force_lowest_estimate=True)
        assert lowest.quantile == 1
        assert lowest.force_lowest_estimate is True
        assert lowest.num_phase_qubits > best.num_phase_qubits
        assert lowest.multi_run_failure <= 1e-6

    def test_cheaper_than_single_run_nielsen_and_chuang(self):
        # At overlap 1 the multi-run plan beats today's single-run rule.
        plan = size_textbook_qpe(FINE_PHASE_ERROR, 1e-6, 1.0)
        single = RectangularWindowState.from_requirements(FINE_PHASE_ERROR, 1e-6)
        assert plan.num_runs * 2**plan.num_phase_qubits < 2**single.m_bits

    def test_qubits_below_nielsen_and_chuang_minimum_allowed(self):
        # Rare (10 of 882 cases surveyed, all with a loose target at overlap 1): one run with one
        # qubit below precision_bits(W) is cheapest.
        plan = size_textbook_qpe(PHASE_ERROR, 0.4, 1.0)
        assert plan.num_phase_qubits == precision_bits(PHASE_ERROR) - 1
        assert (plan.num_runs, plan.quantile) == (1, 1)
        assert plan.multi_run_failure <= 0.4

    def test_large_n_is_fast(self):
        # Took 88 s when the K search stepped one K at a time.
        plan = size_textbook_qpe(FINE_PHASE_ERROR, 1e-9, 1e-3)
        assert (plan.num_phase_qubits, plan.num_runs, plan.quantile) == (20, 115883, 58)
        assert plan.num_runs * 2**plan.num_phase_qubits <= COST_CEILING

    def test_infeasible_within_cost_ceiling(self):
        assert size_textbook_qpe(FINE_PHASE_ERROR, 1e-9, 1e-14) is None

    @pytest.mark.parametrize("overlap", [0.0, 1.5])
    def test_rejects_overlap(self, overlap):
        with pytest.raises(ValueError, match="overlap"):
            size_textbook_qpe(PHASE_ERROR, 1e-3, overlap)

    @pytest.mark.parametrize("probability_of_failure", [0.0, 1.0])
    def test_rejects_probability_of_failure(self, probability_of_failure):
        with pytest.raises(ValueError, match="probability of failure"):
            size_textbook_qpe(PHASE_ERROR, probability_of_failure, 0.5)


class TestFewestRuns:
    def test_given_qubits(self):
        plan = fewest_runs(PHASE_ERROR, 7, 1e-3, 0.4)
        assert plan.num_phase_qubits == 7
        assert plan.multi_run_failure == pytest.approx(
            multi_run_failure(PHASE_ERROR, 7, plan.num_runs, plan.quantile, 0.4), rel=1e-12)
        assert plan.multi_run_failure <= 1e-3
        u = per_run_qpe_failure(PHASE_ERROR, 7)
        hit = 0.4 * (1 - u)
        assert all(worst_failure(plan.num_runs - 1, k, u, hit) > 1e-3
                   for k in range(1, plan.num_runs))

    def test_too_few_qubits(self):
        # 2 qubits: u = 1, so no run is ever trusted.
        assert fewest_runs(PHASE_ERROR, 2, 1e-3, 0.5) is None

    def test_overlap_table_value(self):
        # The plan sized for overlap 0.5, if the true overlap is 1: failure drops.
        plan = size_textbook_qpe(PHASE_ERROR, 1e-3, 0.5)
        at_one = multi_run_failure(PHASE_ERROR, plan.num_phase_qubits, plan.num_runs,
                                   plan.quantile, 1.0)
        assert at_one < plan.multi_run_failure
        assert math.isclose(
            multi_run_failure(PHASE_ERROR, plan.num_phase_qubits, plan.num_runs, plan.quantile,
                              0.5), plan.multi_run_failure, rel_tol=1e-12)


class TestModelAgainstQPE:
    """The multi-run bound against exact textbook QPE for random spectra and initial states, with
    every eigenphase at least the phase error from the window edges."""

    def test_bound_holds_for_random_spectra(self):
        rng = np.random.default_rng(11)
        checked = 0
        for _ in range(60):
            phase_error = float(rng.uniform(0.02, 0.2))
            overlap = float(rng.uniform(0.05, 1.0))
            probability_of_failure = 10**rng.uniform(-6, -1)
            plan = size_textbook_qpe(phase_error, probability_of_failure, overlap)
            if plan is None or plan.num_phase_qubits > 12:
                continue
            num_excited = int(rng.integers(0, 6))
            phases = np.sort(rng.uniform(phase_error, 1 - phase_error, num_excited + 1))
            weights = [1.0] if num_excited == 0 else np.concatenate(
                [[overlap], rng.dirichlet(np.ones(num_excited)) * (1 - overlap)])
            m = 2**plan.num_phase_qubits
            offsets = np.arange(m)
            probs = 0
            for c, phase in zip(weights, phases):
                probs = probs + c * np.array(
                    [exact_failure_outcome(plan.num_phase_qubits, phase, d) for d in offsets])
            estimates = offsets / m
            distance = np.abs(estimates - phases[0])
            accepted = np.minimum(distance, 1 - distance) <= phase_error
            under = ~accepted & (estimates < phases[0])
            failure = (binom.sf(plan.quantile - 1, plan.num_runs, probs[under].sum())
                       + binom.sf(plan.num_runs - plan.quantile, plan.num_runs,
                                  probs[~accepted & ~under].sum()))
            assert failure <= plan.multi_run_failure * (1 + 1e-9)
            checked += 1
        assert checked >= 30


def exact_failure_outcome(num_phase_qubits, phase, outcome):
    """Probability that textbook QPE with `num_phase_qubits` returns `outcome` for `phase`."""
    m = 2**num_phase_qubits
    x = phase * m - outcome
    if abs(np.sin(np.pi * x / m)) < 1e-15:
        return 1.0
    return float(np.sin(np.pi * x)**2 / (m**2 * np.sin(np.pi * x / m)**2))
