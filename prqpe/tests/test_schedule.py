"""Property tests for the robust-phase-estimation schedule."""

import math

import pytest

from qhat.prqpe import schedule
from qhat.prqpe.schedule import (
    RPESchedule,
    num_rounds,
    num_rounds_randomized,
    resolve_depth_tradeoff,
    round_circuit_toffolis,
    rpe_schedule,
)


# ----------------------------------------------------------------------------
# Depth-reduction factor xi.
# ----------------------------------------------------------------------------

def test_xi_none_means_no_reduction():
    assert resolve_depth_tradeoff(overlap=1.0, xi=None) == 1.0


def test_xi_passes_through_for_exact_eigenstate():
    assert resolve_depth_tradeoff(overlap=1.0, xi=0.1) == 0.1


def test_xi_out_of_range_rejected():
    with pytest.raises(ValueError):
        resolve_depth_tradeoff(overlap=1.0, xi=1.5)


def test_xi_below_overlap_floor_rejected():
    # For overlap < 1, xi must exceed arcsin((1-overlap)/overlap).
    with pytest.raises(ValueError):
        resolve_depth_tradeoff(overlap=0.9, xi=1e-4)


def test_low_overlap_rejected():
    with pytest.raises(ValueError):
        resolve_depth_tradeoff(overlap=0.4, xi=0.5)


# ----------------------------------------------------------------------------
# Round count.
# ----------------------------------------------------------------------------

def test_num_rounds_is_ceil_log2():
    eps_qpe, delta = 1.3e-3, 0.02
    expected = math.ceil(math.log2(0.08 * math.pi / (eps_qpe * delta)))
    assert num_rounds(eps_qpe, delta) == expected


def test_xi_reduces_round_count():
    eps_qpe, delta = 1.3e-3, 0.02
    assert num_rounds(eps_qpe, delta, xi=0.1) <= num_rounds(eps_qpe, delta, xi=1.0)


def test_randomized_round_count_matches_paper_femoco():
    # FeMoco fully randomized at xi = 0.1: M = ceil(log2(xi*lambda/eps)) = 15.
    assert num_rounds_randomized(405.0, 0.0016, xi=0.1) == 15


# ----------------------------------------------------------------------------
# Sample schedule.
# ----------------------------------------------------------------------------

def test_samples_decrease_toward_deepest_round():
    sched = RPESchedule(M=10, inflation=1.0)
    counts = [sched.samples(m) for m in range(sched.M + 1)]
    assert all(a >= b for a, b in zip(counts, counts[1:]))
    assert sched.samples(sched.M) == 11  # N_M at inflation 1


def test_num_circuits_is_twice_the_sample_sum():
    sched = RPESchedule(M=8, inflation=math.e)
    assert sched.num_circuits == 2 * sum(sched.samples(m) for m in range(9))


def test_rte_inflation_is_e():
    # RTE optimal kappa = 2 inflates samples by e^(2/2) = e.
    sched = rpe_schedule(M=5, kappa=2.0, overlap=1.0)
    assert sched.inflation == pytest.approx(math.e)


def test_deterministic_has_no_inflation():
    sched = rpe_schedule(M=5, kappa=math.inf, overlap=1.0)
    assert sched.inflation == pytest.approx(1.0)


def test_lower_overlap_inflates_samples():
    high = rpe_schedule(M=5, kappa=2.0, overlap=1.0)
    low = rpe_schedule(M=5, kappa=2.0, overlap=0.8)
    assert low.inflation > high.inflation


# ----------------------------------------------------------------------------
# Per-circuit cost.
# ----------------------------------------------------------------------------

def test_deterministic_depth_uses_halving():
    # The deterministic part runs 2^(m-1) halved Trotter steps.
    circuit = dict(G_det=23.0, G_rand=2.5, L_D=100, lam_R=0.0,
                   delta=0.02, kappa=math.inf)
    assert round_circuit_toffolis(3, **circuit) == pytest.approx(
        23.0 * 2 * 100 * 2.0 ** 2)


def test_deepest_circuit_is_round_M():
    sched = RPESchedule(M=6, inflation=math.e)
    circuit = dict(G_det=23.0, G_rand=2.5, L_D=500, lam_R=30.0,
                   delta=0.02, kappa=4.0)
    deepest = schedule.max_toffolis_per_circuit(sched, **circuit)
    assert deepest == round_circuit_toffolis(6, **circuit)
    assert deepest >= round_circuit_toffolis(5, **circuit)
