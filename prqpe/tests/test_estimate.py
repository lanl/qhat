"""Integration and property tests for ``estimate_resources`` through the public API."""

import math

import numpy as np
import pytest

from qhat.prqpe import estimate_resources
from qhat.prqpe.costs import (
    minimize_delta_kappa,
    phase_state_precision,
    rotation_toffolis,
)

EPS = 0.0016
K = 10


# Closed-form totals for the degenerate splits (the L_D=0 / lam_R=0 cases of
# minimize_delta_kappa), used here as cross-check oracles for the endpoints.
def _randomized_total(lam, eps, G_rand):
    return minimize_delta_kappa(0, lam, eps, 0.1, G_det=0.0, G_rand=G_rand).toffoli_total


def _deterministic_total(L, eps, C_gs, G_det):
    return minimize_delta_kappa(L, 0.0, eps, C_gs, G_det=G_det, G_rand=2.5).toffoli_total


def representative_weights(lam: float, n_terms: int, seed: int = 0) -> np.ndarray:
    """A heavy-tailed positive weight array summing to ``lam``.

    Stands in for a molecular Pauli spectrum: a few large coefficients and a
    long tail of small ones.  Not any specific molecule.
    """
    rng = np.random.default_rng(seed)
    w = rng.pareto(1.5, size=n_terms) + 1e-3
    return np.sort(w)[::-1] * (lam / w.sum())


# ----------------------------------------------------------------------------
# Endpoint methods reduce to the closed-form arms.
# ----------------------------------------------------------------------------

def test_randomized_method_matches_closed_form():
    w = representative_weights(405.0, 3000)
    est = estimate_resources(w, 152, C_gs=0.1, method="randomized")
    G_rand = rotation_toffolis(phase_state_precision(405.0, EPS), K)
    assert est.toffoli_count == pytest.approx(
        _randomized_total(405.0, EPS, G_rand))
    assert est.metadata["J"] == 17
    assert est.metadata["kappa"] == pytest.approx(2.0)
    assert est.logical_qubits == 195  # 152 + 1 + 42
    assert len(est.pareto) == 1


def test_deterministic_method_has_no_randomized_part():
    w = representative_weights(405.0, 3000)
    est = estimate_resources(w, 152, C_gs=0.1, method="deterministic")
    assert est.metadata["lambda_R"] == 0.0
    assert math.isinf(est.metadata["kappa"])
    assert est.logical_qubits == 153  # system + Hadamard-test ancilla only
    assert est.toffoli_count == pytest.approx(
        _deterministic_total(len(w), EPS, 0.1, est.metadata["G_det"]))


# ----------------------------------------------------------------------------
# The partially randomized sweep.
# ----------------------------------------------------------------------------

def test_partial_sweep_returns_full_pareto_curve():
    w = representative_weights(405.0, 4000)
    est = estimate_resources(w, 152, C_gs=0.1, method="partially_randomized")
    assert len(est.pareto) > 10
    L_D = [p.L_D for p in est.pareto]
    assert L_D == sorted(L_D)
    assert L_D[0] == 0 and L_D[-1] == len(w)


def test_partial_optimum_beats_both_endpoints():
    w = representative_weights(405.0, 4000)
    est = estimate_resources(w, 152, C_gs=0.1, method="partially_randomized")
    randomized = est.pareto[0].toffoli_total       # L_D = 0
    deterministic = est.pareto[-1].toffoli_total    # L_D = L
    assert est.toffoli_count <= randomized
    assert est.toffoli_count <= deterministic


def test_partial_endpoints_match_closed_forms():
    w = representative_weights(405.0, 4000)
    est = estimate_resources(w, 152, C_gs=0.1, method="partially_randomized")
    G_rand = rotation_toffolis(phase_state_precision(405.0, EPS), K)
    assert est.pareto[0].toffoli_total == pytest.approx(
        _randomized_total(405.0, EPS, G_rand))
    assert est.pareto[-1].lambda_R == 0.0


def test_sweep_is_unimodal():
    # Cost as a function of L_D falls to a single interior minimum and rises.
    w = representative_weights(405.0, 4000)
    est = estimate_resources(w, 152, C_gs=0.1, method="partially_randomized")
    totals = np.array([p.toffoli_total for p in est.pareto])
    signs = np.sign(np.diff(totals))
    signs = signs[signs != 0]
    sign_changes = int(np.sum(np.diff(signs) != 0))
    assert sign_changes <= 1  # down-then-up: one change at the minimum


def test_partial_improves_on_fully_randomized_at_scale():
    # A FeMoco-scale problem: the optimum strictly improves on full randomization.
    w = representative_weights(405.0, 4000)
    est = estimate_resources(w, 152, C_gs=0.1, method="partially_randomized")
    G_rand = rotation_toffolis(17, K)
    assert est.toffoli_count < _randomized_total(405.0, EPS, G_rand)
    assert est.metadata["L_D"] not in (0, len(w))  # genuinely partial


# ----------------------------------------------------------------------------
# Reporting and validation.
# ----------------------------------------------------------------------------

def test_metadata_is_complete():
    w = representative_weights(405.0, 2000)
    est = estimate_resources(w, 152, C_gs=0.1)
    for key in ("L_D", "lambda_R", "delta", "M", "kappa", "xi",
                "eps_qpe", "J", "K", "G_det", "G_rand"):
        assert key in est.metadata


def test_qdrift_is_cheaper_than_rte():
    w = representative_weights(405.0, 3000)
    rte = estimate_resources(w, 152, C_gs=0.1, method="randomized", randomizer="rte")
    qdrift = estimate_resources(w, 152, C_gs=0.1, method="randomized",
                                randomizer="qdrift")
    assert qdrift.toffoli_count == pytest.approx(0.5 * rte.toffoli_count)


def test_unknown_representation_rejected():
    w = representative_weights(405.0, 100)
    with pytest.raises(ValueError):
        estimate_resources(w, 152, C_gs=0.1, representation="double_factorized")


def test_unknown_method_rejected():
    w = representative_weights(405.0, 100)
    with pytest.raises(ValueError):
        estimate_resources(w, 152, C_gs=0.1, method="nonsense")


def test_empty_weights_rejected():
    with pytest.raises(ValueError):
        estimate_resources([], 152, C_gs=0.1)
