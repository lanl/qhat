"""Paper-reproduction anchors and closed-form property tests for the cost core.

The reference points are the worked FeMoco numbers and the equations they come
from.  The FeMoco system has weight ``lambda = 405`` and is estimated to
precision ``eps = 0.0016`` Hartree with commuting-group length ``K = 10``.
"""

import math

import numpy as np
import pytest

from qhat.prqpe import costs
from qhat.prqpe.costs import (
    minimize_delta_kappa,
    optimal_kappa,
    partial_toffoli_bound,
    phase_state_precision,
    qpe_trotter_split,
    rotation_toffolis,
    trotter_step_size,
)

EPS = 0.0016
LAM_FEMOCO = 405.0
K = 10


# Closed-form totals for the two degenerate splits, used here as cross-check
# oracles.  They are exactly the L_D=0 / lam_R=0 cases of minimize_delta_kappa.
def _randomized_total(lam, eps, G_rand):
    return minimize_delta_kappa(0, lam, eps, 0.1, G_det=0.0, G_rand=G_rand).toffoli_total


def _deterministic_total(L, eps, C_gs, G_det):
    return minimize_delta_kappa(L, 0.0, eps, C_gs, G_det=G_det, G_rand=2.5).toffoli_total


# ----------------------------------------------------------------------------
# FeMoco paper anchors (weight-independent: derivable from lambda and eps).
# ----------------------------------------------------------------------------

def test_femoco_dyadic_exponent():
    # J = ceil(log2(lambda / (pi*eps))).
    assert phase_state_precision(LAM_FEMOCO, EPS) == 17


def test_femoco_toffolis_per_rotation():
    # 1 + (J-2)/K = 1 + 15/10 = 2.5.
    assert rotation_toffolis(17, K) == pytest.approx(2.5)


def test_dyadic_exponent_floored_at_one_bit():
    # A randomized weight below the precision scale still needs one bit, keeping
    # the per-rotation cost positive.
    assert phase_state_precision(1e-9, EPS) == 1
    assert rotation_toffolis(phase_state_precision(1e-9, EPS), K) > 0.0


def test_femoco_phasing_ancillas():
    # K + 2J - 2 = 10 + 34 - 2 = 42.
    assert costs.phasing_ancillas(17, K) == 42


def test_femoco_logical_qubits():
    # 152 system qubits + 1 Hadamard-test ancilla + 42 phasing ancillas.
    assert costs.logical_qubits(152, 17, K) == 195


def test_femoco_fully_randomized_total():
    # The fully randomized arm depends only on lambda; its value sets the scale
    # the partially randomized estimate improves on (paper: ~2.8e12 Toffoli,
    # the upper method in the FeMoco comparison figure).
    G_rand = rotation_toffolis(17, K)
    total = _randomized_total(LAM_FEMOCO, EPS, G_rand)
    assert total == pytest.approx(2.83e12, rel=0.05)


def test_fully_randomized_matches_rte_rotation_bookkeeping():
    # Independent cross-check: the App. E.2 randomized bookkeeping counts
    # r = 16.3 lambda^2/eps^2 Pauli rotations, so G_rand*r Toffolis.  This must
    # agree with the schedule-summed fully randomized total to ~10%.
    G_rand = rotation_toffolis(17, K)
    rotation_count = 16.3 * LAM_FEMOCO ** 2 / EPS ** 2
    bookkeeping = G_rand * rotation_count
    closed_form = _randomized_total(LAM_FEMOCO, EPS, G_rand)
    assert closed_form == pytest.approx(bookkeeping, rel=0.1)


# ----------------------------------------------------------------------------
# Error-budget and Trotter-step consistency.
# ----------------------------------------------------------------------------

def test_qpe_trotter_split_adds_in_quadrature():
    eps_qpe, eps_trot = qpe_trotter_split(EPS, order=2)
    assert eps_qpe ** 2 + eps_trot ** 2 == pytest.approx(EPS ** 2)
    assert eps_qpe == pytest.approx(EPS * math.sqrt(2.0 / 3.0))


@pytest.mark.parametrize("C_gs", [1e-3, 0.1, 1.4])
def test_trotter_step_is_self_consistent(C_gs):
    # The Trotter step must satisfy C_gs*delta^p = eps_trotter exactly, so that
    # eps^2 = eps_qpe^2 + eps_trotter^2 holds at the optimum.
    delta = trotter_step_size(EPS, C_gs, order=2)
    _, eps_trot = qpe_trotter_split(EPS, order=2)
    assert C_gs * delta ** 2 == pytest.approx(eps_trot)


# ----------------------------------------------------------------------------
# kappa optimization.
# ----------------------------------------------------------------------------

@pytest.mark.parametrize("A,B", [(1.0, 1.0), (10.0, 1.0), (0.5, 3.0)])
def test_optimal_kappa_solves_quadratic(A, B):
    kappa = optimal_kappa(A, B)
    # dG/dkappa = 0  <=>  B*kappa^2 - 2*B*kappa - 2*A = 0.
    assert B * kappa ** 2 - 2.0 * B * kappa - 2.0 * A == pytest.approx(0.0, abs=1e-9)


def test_optimal_kappa_no_randomized_part_is_infinite():
    assert math.isinf(optimal_kappa(1.0, 0.0))


def test_optimal_kappa_no_deterministic_part_is_two():
    # A = 0 recovers the RTE-alone optimum r = 2 t^2.
    assert optimal_kappa(0.0, 5.0) == pytest.approx(2.0)


# ----------------------------------------------------------------------------
# Split endpoints: the optimizer reduces to the closed-form arms.
# ----------------------------------------------------------------------------

def test_fully_randomized_endpoint():
    G_rand = rotation_toffolis(17, K)
    opt = minimize_delta_kappa(0, LAM_FEMOCO, EPS, 0.1, G_det=0.0, G_rand=G_rand)
    assert opt.kappa == pytest.approx(2.0)
    assert opt.eps_qpe == pytest.approx(EPS)
    assert opt.toffoli_total == pytest.approx(
        _randomized_total(LAM_FEMOCO, EPS, G_rand))


def test_fully_deterministic_endpoint():
    L, G_det = 5000, 23.0
    opt = minimize_delta_kappa(L, 0.0, EPS, 0.1, G_det=G_det, G_rand=2.5)
    assert math.isinf(opt.kappa)
    assert opt.eps_qpe == pytest.approx(EPS * math.sqrt(2.0 / 3.0))
    assert opt.toffoli_total == pytest.approx(
        _deterministic_total(L, EPS, 0.1, G_det))


def test_bound_matches_optimizer_at_returned_point():
    # Re-evaluating the closed-form bound at the optimizer's (delta, kappa)
    # reproduces its reported total.
    opt = minimize_delta_kappa(2000, 40.0, EPS, 0.1, G_det=23.0, G_rand=2.5)
    recomputed = partial_toffoli_bound(
        2000, 40.0, opt.delta, opt.eps_qpe, opt.kappa, 23.0, 2.5)
    assert recomputed == pytest.approx(opt.toffoli_total)


# ----------------------------------------------------------------------------
# Asymptotic scaling with target precision.
# ----------------------------------------------------------------------------

def test_fully_randomized_scales_as_inverse_eps_squared():
    G_rand = 2.5
    coarse = _randomized_total(LAM_FEMOCO, EPS, G_rand)
    fine = _randomized_total(LAM_FEMOCO, EPS / 2.0, G_rand)
    assert fine / coarse == pytest.approx(4.0)


def test_fully_deterministic_scales_as_eps_to_the_minus_three_halves():
    L, G_det = 5000, 23.0
    coarse = _deterministic_total(L, EPS, 0.1, G_det)
    fine = _deterministic_total(L, EPS / 2.0, 0.1, G_det)
    assert fine / coarse == pytest.approx(2.0 ** 1.5)
