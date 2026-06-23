"""Closed-form cost primitives for (partially) randomized phase estimation.

Every quantity in this module is an analytic resource formula for ground-state
energy estimation with second-order product formulas that are *deterministic*,
*randomized*, or a *partial* mix of the two.

The cost model has three layers, each a small group of functions below:

1.  Product-formula error control -- splitting the precision budget between
    phase estimation and Trotter error, and the resulting Trotter step size.
2.  Per-operation gate costs -- Toffolis per Pauli rotation (Hamming-weight
    phasing with a shared catalyst phase state) and per arbitrary-angle
    rotation synthesis.
3.  The partially randomized total -- a closed-form gate bound at a fixed
    deterministic/randomized split, jointly optimized over the Trotter step
    ``delta`` and the randomization scale ``kappa``.

All formulas take the *unnormalized* Hamiltonian 1-norm ``lam`` and the target
precision ``eps`` in the same (Hartree) units.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# ----------------------------------------------------------------------------
# Empirical robust-phase-estimation sample-schedule constants.
#
# Robust phase estimation runs rounds m = 0, ..., M and takes
# N_m = N_M + D*(M - m) samples per quadrature in round m.  The numerically
# calibrated values for an exact eigenstate are N_M = 11 in the deepest round
# and slope D = 4.  These two constants set the closed-form prefactors of the
# total-gate bound below.
# ----------------------------------------------------------------------------

RPE_LAST_ROUND_SAMPLES = 11.0  # N_M: samples per quadrature in the deepest round
RPE_SLOPE = 4.0  # D: extra samples per quadrature per round below the deepest

#: Empirical depth constant: the deepest round reaches t_max with
#: eps_qpe ~= C_MAX_TIME * pi / t_max, so M = ceil(log2(C_MAX_TIME*pi/(eps*delta))).
C_MAX_TIME = 0.08

#: Depth constant in the closed-form gate bound: 2^M = C_MAX_GATE*pi/(delta*eps_qpe).
C_MAX_GATE = 0.1

#: Two T gates are counted as one Toffoli gate.
T_GATES_PER_TOFFOLI = 2.0


def _gate_bound_prefactors(N_M: float = RPE_LAST_ROUND_SAMPLES,
                           D: float = RPE_SLOPE) -> tuple[float, float]:
    """Closed-form prefactors (deterministic, randomized) of the total bound.

    Summing the per-round sample schedule against the per-round gate counts
    gives ``2*(N_M + D)`` for the deterministic term and ``8*(N_M/3 + D/9)``
    for the randomized term.  For (N_M, D) = (11, 4) these are 30 and 296/9.
    """
    return 2.0 * (N_M + D), 8.0 * (N_M / 3.0 + D / 9.0)


# ----------------------------------------------------------------------------
# 1. Product-formula error control.
# ----------------------------------------------------------------------------

def qpe_trotter_split(eps: float, order: int = 2) -> tuple[float, float]:
    """Split the precision budget as eps^2 = eps_qpe^2 + eps_trotter^2.

    The phase-estimation error is unbiased and independent of the Trotter
    error, so the two add in quadrature.  Minimizing the total gate count
    under this constraint puts ``eps_qpe = eps*sqrt(p/(p+1))`` and hence
    ``eps_trotter = eps/sqrt(p+1)``.
    """
    eps_qpe = eps * math.sqrt(order / (order + 1.0))
    eps_trotter = eps / math.sqrt(order + 1.0)
    return eps_qpe, eps_trotter


def trotter_step_size(eps: float, C_gs: float, order: int = 2) -> float:
    """Cost-optimal Trotter step ``delta = (eps / (C_gs*sqrt(p+1)))^(1/p)``.

    The ground-state Trotter bias obeys ``|E0 - E0~| = C_gs * delta^p``, so the
    Trotter share ``eps_trotter = eps/sqrt(p+1)`` of the budget fixes
    ``delta^p = eps / (C_gs*sqrt(p+1))``.  Equivalently
    ``delta = (eps/C_gs)^(1/p) * (1/(p+1))^(1/(2p))``; this is the step that is
    self-consistent with ``eps^2 = eps_qpe^2 + eps_trotter^2`` and with the
    per-rotation synthesis allocation in :func:`synthesis_precision_per_rotation`.
    """
    return (eps / (C_gs * math.sqrt(order + 1.0))) ** (1.0 / order)


# ----------------------------------------------------------------------------
# 2. Per-operation gate costs.
# ----------------------------------------------------------------------------

def phase_state_precision(lam: float, eps: float) -> int:
    """Dyadic-angle exponent ``J = ceil(log2(lam / (pi*eps)))`` (at least one bit).

    The randomized product formula is free to choose its step so that every
    Pauli rotation is over the same angle ``pi*2^-J``, implemented from a shared
    ``J``-qubit catalyst phase state.  ``J`` is fixed by requiring the step to be
    ``O(eps/lam)``.  When the weight is below the precision scale (``lam < pi*eps``)
    a single bit suffices and the randomized contribution is negligible.
    """
    return max(1, int(math.ceil(math.log2(lam / (math.pi * eps)))))


def rotation_toffolis(J: int, K: int) -> float:
    """Toffolis per Pauli rotation with Hamming-weight phasing: ``1 + (J-2)/K``.

    A run of ``K`` mutually commuting Paulis is mapped to single-qubit Z
    rotations and phased through the catalyst state: ``K-1`` Toffolis compute and
    uncompute the Hamming weight and ``J-1`` Toffolis apply the phase, i.e.
    ``((K-1)+(J-1))/K = 1 + (J-2)/K`` Toffolis per rotation.
    """
    return 1.0 + (J - 2.0) / K


def phasing_ancillas(J: int, K: int) -> int:
    """Ancilla qubits for Hamming-weight phasing: ``K + 2J - 2``.

    A ``J``-qubit catalyst phase state, ``K-1`` Hamming-weight register qubits,
    and ``J-1`` adder ancillas: ``J + (K-1) + (J-1) = K + 2J - 2``.
    """
    return K + 2 * J - 2


def logical_qubits(n_system: int, J: int, K: int) -> int:
    """Total logical qubits: system + one phase-estimation ancilla + phasing.

    ``n_system + 1 + (K + 2J - 2)``.  The single extra ancilla carries the
    Hadamard test; the remaining ``K + 2J - 2`` are the Hamming-weight phasing
    register, adder ancillas, and catalyst phase state of :func:`phasing_ancillas`.
    """
    return n_system + 1 + phasing_ancillas(J, K)


def synthesis_toffolis(eps_prime: float) -> float:
    """Toffolis to synthesize one arbitrary-angle rotation to error ``eps'``.

    A rotation costs ``1.14*log2(1/eps') + 9.2`` T gates, i.e. half that many
    Toffolis at two T per Toffoli.
    """
    t_gates = 1.14 * math.log2(1.0 / eps_prime) + 9.2
    return t_gates / T_GATES_PER_TOFFOLI


def synthesis_precision_per_rotation(delta: float, eps_synth: float,
                                     n_stage: int, n_rot: int) -> float:
    """Per-rotation synthesis error ``eps' = delta*eps_synth/(n_stage*n_rot)``.

    Imperfect synthesis of every rotation perturbs the effective Hamiltonian by
    at most ``n_stage*n_rot*eps'/delta``; keeping the induced ground-state error
    below ``eps_synth`` fixes ``eps'`` as above.  ``n_rot`` is the number of
    rotations per Trotter stage (the deterministic term count for a Pauli
    Hamiltonian).
    """
    return delta * eps_synth / (n_stage * n_rot)


# ----------------------------------------------------------------------------
# 3. The partially randomized total.
# ----------------------------------------------------------------------------

@dataclass(frozen=True)
class SplitOptimum:
    """Cost-optimal operating point at a fixed deterministic/randomized split."""

    toffoli_total: float
    delta: float
    kappa: float
    eps_qpe: float


def optimal_kappa(A: float, B: float) -> float:
    """Minimizer of ``G(kappa) = (A + B*kappa) * e^(2/kappa)`` over ``kappa > 0``.

    ``kappa`` scales the randomized rotation budget ``r = kappa*lam_R^2*delta^2*4^m``
    and sets the damping ``B_damp <= e^(1/kappa)`` that inflates every sample
    count by ``e^(2/kappa)``.  Setting ``dG/dkappa = 0`` gives the quadratic
    ``B*kappa^2 - 2*B*kappa - 2*A = 0`` with positive root
    ``kappa* = 1 + sqrt(1 + 2A/B)``.  ``A`` is the deterministic prefactor
    (no randomized part -> ``kappa* -> inf``); ``A = 0`` recovers the fully
    randomized optimum ``kappa* = 2``.
    """
    if B == 0.0:
        return math.inf
    return 1.0 + math.sqrt(1.0 + 2.0 * A / B)


def partial_toffoli_bound(L_D: int, lam_R: float, delta: float, eps_qpe: float,
                          kappa: float, G_det: float, G_rand: float,
                          n_stage: int = 2, randomizer_factor: float = 1.0,
                          N_M: float = RPE_LAST_ROUND_SAMPLES,
                          D: float = RPE_SLOPE) -> float:
    """Closed-form total Toffoli bound at a fixed split and operating point.

    Summing ``sum_m 2*N_m*(G_det*n_stage*L_D*2^(m-1) + G_rand*kappa*lam_R^2*delta^2*4^m)``
    over the robust-phase-estimation rounds, with ``N_m = e^(2/kappa)*(N_M + D*(M-m))``
    and ``2^M = C_MAX_GATE*pi/(delta*eps_qpe)``, gives

        G = 2*(N_M+D) * G_det*n_stage*L_D * e^(2/kappa) * C_MAX_GATE*pi/(delta*eps_qpe)
          + 8*(N_M/3+D/9) * G_rand*kappa*e^(2/kappa) * (C_MAX_GATE*pi*lam_R)^2/eps_qpe^2.

    The deterministic term already includes the factor-2 reduction from halving
    the deterministic evolution time (``2^(m-1)`` rather than ``2^m``).  ``L_D = 0``
    drops the deterministic term; ``lam_R = 0`` (``kappa = inf``) drops the
    randomized term.  ``randomizer_factor`` scales the randomized rotation count
    (``1`` for RTE, ``1/2`` for qDRIFT).
    """
    det_pre, rand_pre = _gate_bound_prefactors(N_M, D)
    inflation = 1.0 if math.isinf(kappa) else math.exp(2.0 / kappa)

    det = 0.0
    if L_D > 0 and delta > 0.0:
        det = (det_pre * G_det * n_stage * L_D * inflation
               * C_MAX_GATE * math.pi / (delta * eps_qpe))
    rand = 0.0
    if lam_R > 0.0:
        rand = (rand_pre * randomizer_factor * G_rand * kappa * inflation
                * (C_MAX_GATE * math.pi * lam_R) ** 2 / eps_qpe ** 2)
    return det + rand


def minimize_delta_kappa(L_D: int, lam_R: float, eps: float, C_gs: float,
                         G_det: float, G_rand: float, *, order: int = 2,
                         n_stage: int = 2, randomizer_factor: float = 1.0,
                         n_grid: int = 400,
                         N_M: float = RPE_LAST_ROUND_SAMPLES,
                         D: float = RPE_SLOPE) -> SplitOptimum:
    """Jointly optimize the total bound over ``(delta, kappa)`` at fixed ``L_D``.

    Substituting ``eps_qpe = sqrt(eps^2 - C_gs^2*delta^4)`` couples the Trotter
    step to the phase-estimation budget.  Parameterize ``x = C_gs*delta^2/eps``
    in ``(0, 1)``, so ``delta = sqrt(x*eps/C_gs)`` and ``eps_qpe = eps*sqrt(1-x^2)``;
    at each grid point ``kappa`` is the closed-form root of :func:`optimal_kappa`,
    and the global minimum over the grid is returned.

    The two degenerate splits are exact closed forms (no grid):

    * ``L_D = 0`` (fully randomized): ``kappa = 2``, ``eps_qpe = eps``,
      ``delta`` immaterial.
    * ``lam_R = 0`` (fully deterministic): ``kappa = inf`` and the optimum sits
      at ``x = 1/sqrt(3)``, i.e. ``eps_qpe = eps*sqrt(2/3)``.
    """
    if order != 2:
        raise NotImplementedError(
            "the partially randomized cost model is implemented for order 2")

    if L_D == 0:
        total = partial_toffoli_bound(0, lam_R, 0.0, eps, 2.0, G_det, G_rand,
                                      n_stage=n_stage,
                                      randomizer_factor=randomizer_factor,
                                      N_M=N_M, D=D)
        return SplitOptimum(total, delta=0.0, kappa=2.0, eps_qpe=eps)

    if lam_R == 0.0:
        x = 1.0 / math.sqrt(3.0)
        delta = math.sqrt(x * eps / C_gs)
        eps_qpe = eps * math.sqrt(1.0 - x ** 2)
        total = partial_toffoli_bound(L_D, 0.0, delta, eps_qpe, math.inf,
                                      G_det, G_rand, n_stage=n_stage,
                                      randomizer_factor=randomizer_factor,
                                      N_M=N_M, D=D)
        return SplitOptimum(total, delta=delta, kappa=math.inf, eps_qpe=eps_qpe)

    # Grid in x, denser near 0, with the analytic fully-deterministic point
    # x = 1/sqrt(3) included so the interior optimum is never missed.
    x = np.linspace(1e-4, 1.0 - 1e-6, n_grid) ** 2
    x = np.sort(np.append(x, 1.0 / math.sqrt(3.0)))
    delta = np.sqrt(x * eps / C_gs)
    eps_qpe = eps * np.sqrt(1.0 - x ** 2)

    det_pre, rand_pre = _gate_bound_prefactors(N_M, D)
    A = det_pre * G_det * n_stage * L_D * C_MAX_GATE * math.pi / (delta * eps_qpe)
    B = (rand_pre * randomizer_factor * G_rand
         * (C_MAX_GATE * math.pi * lam_R) ** 2 / eps_qpe ** 2)
    kappa = 1.0 + np.sqrt(1.0 + 2.0 * A / B)
    total = (A + B * kappa) * np.exp(2.0 / kappa)

    total = np.where(np.isfinite(total), total, np.inf)
    i = int(np.argmin(total))
    return SplitOptimum(float(total[i]), float(delta[i]), float(kappa[i]),
                        float(eps_qpe[i]))


def randomized_deepest_circuit_toffolis(total_randomized_toffolis: float,
                                        N_M: float = RPE_LAST_ROUND_SAMPLES,
                                        D: float = RPE_SLOPE) -> float:
    """Toffolis in the deepest randomized circuit, from the randomized total.

    The schedule sums to ``r_tot = 8*(N_M/3 + D/9) * r_max``, so the deepest
    single circuit carries ``total / (8*(N_M/3 + D/9))`` of the Toffolis.
    """
    _, rand_pre = _gate_bound_prefactors(N_M, D)
    return total_randomized_toffolis / rand_pre
