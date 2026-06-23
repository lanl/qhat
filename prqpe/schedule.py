"""Robust-phase-estimation round and sample schedule.

Phase estimation runs rounds ``m = 0, ..., M`` at evolution time ``t = 2^m``,
taking ``N_m`` samples per quadrature in round ``m``.  This module owns the
round count ``M``, the sample schedule ``N_m = N_M + D*(M - m)`` (inflated by the
randomized-formula damping), the resulting circuit count, and the per-circuit
Toffoli cost from which the deepest circuit is read off.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .costs import C_MAX_TIME, RPE_LAST_ROUND_SAMPLES, RPE_SLOPE

#: Lowest ground-state overlap for which robust phase estimation is guaranteed.
MIN_OVERLAP = 4.0 - 2.0 * math.sqrt(3.0)  # ~= 0.536


def resolve_depth_tradeoff(overlap: float = 1.0, xi: float | None = None) -> float:
    """Return a valid depth-reduction factor ``xi`` in ``(0, 1]``.

    A smaller ``xi`` shortens the deepest circuit (``t_max ~ xi/eps``) at the
    cost of more circuits.  Robust phase estimation admits depth reduction only
    when the guiding-state overlap is high: ``xi`` must exceed
    ``arcsin((1-overlap)/overlap)``, which is ``0`` for an exact eigenstate
    (``overlap = 1``).  ``xi = None`` means no reduction and returns ``1.0``.
    """
    if xi is None:
        return 1.0
    if not 0.0 < xi <= 1.0:
        raise ValueError(f"xi must lie in (0, 1], got {xi}")
    if overlap <= MIN_OVERLAP:
        raise ValueError(
            f"overlap must exceed {MIN_OVERLAP:.3f} for robust phase estimation")
    if overlap < 1.0:
        xi_min = math.asin((1.0 - overlap) / overlap)
        if xi <= xi_min:
            raise ValueError(
                f"xi must exceed arcsin((1-overlap)/overlap) = {xi_min:.4f}")
    return xi


def _overlap_inflation(overlap: float) -> float:
    """Sample-count inflation from imperfect overlap, normalized to ``1`` at ``eta=1``.

    The per-round sample count scales as ``beta^-2`` with the signal headroom
    ``beta = overlap*(1 + sin(pi/3)) - 1``.  The empirical schedule constants are
    calibrated at ``overlap = 1`` (``beta = sin(pi/3)``), so the relative
    inflation for general overlap is ``(beta(1)/beta(overlap))^2``.
    """
    if overlap >= 1.0:
        return 1.0
    if overlap <= MIN_OVERLAP:
        raise ValueError(
            f"overlap must exceed {MIN_OVERLAP:.3f} for robust phase estimation")
    headroom = overlap * (1.0 + math.sin(math.pi / 3.0)) - 1.0
    return (math.sin(math.pi / 3.0) / headroom) ** 2


def num_rounds(eps_qpe: float, delta: float, xi: float = 1.0) -> int:
    """Round count for a Trotter-stepped arm: ``M = ceil(log2(C*pi*xi/(eps_qpe*delta)))``.

    The deepest round reaches ``t_max = 2^M`` with ``eps_qpe ~= C_MAX_TIME*pi/t_max``
    in units of the Trotter step ``delta``; the depth-reduction factor ``xi``
    shrinks ``M`` when the guiding-state overlap is high.
    """
    return int(math.ceil(math.log2(C_MAX_TIME * math.pi * xi / (eps_qpe * delta))))


def num_rounds_randomized(lam_R: float, eps_qpe: float, xi: float = 1.0) -> int:
    """Round count for the fully randomized arm: ``M = ceil(log2(xi*lam_R/eps_qpe))``.

    With no Trotter step the deepest round estimates the normalized Hamiltonian
    ``H/lam_R`` to precision ``eps_qpe/lam_R`` at time ``t_max = 2^M``.
    """
    return int(math.ceil(math.log2(xi * lam_R / eps_qpe)))


@dataclass(frozen=True)
class RPESchedule:
    """Sample schedule ``N_m = ceil(inflation*(N_M + D*(M - m)))`` over rounds ``0..M``."""

    M: int
    inflation: float
    N_M: float = RPE_LAST_ROUND_SAMPLES
    D: float = RPE_SLOPE

    def samples(self, m: int) -> int:
        """Samples per quadrature in round ``m`` (decreasing toward the deepest)."""
        return int(math.ceil(self.inflation * (self.N_M + self.D * (self.M - m))))

    @property
    def num_circuits(self) -> int:
        """Total Hadamard-test circuits ``sum_m 2*N_m`` (two quadratures per round)."""
        return int(2 * sum(self.samples(m) for m in range(self.M + 1)))


def rpe_schedule(M: int, kappa: float, overlap: float = 1.0) -> RPESchedule:
    """Build the sample schedule for ``M`` rounds at randomization scale ``kappa``.

    The randomized product formula damps the signal by ``B <= e^(1/kappa)``, so
    every sample count is inflated by ``B^2 = e^(2/kappa)``; a fully
    deterministic arm (``kappa = inf``) has no inflation.  Imperfect overlap
    adds a further ``(beta(1)/beta(overlap))^2`` factor.
    """
    damping = 1.0 if math.isinf(kappa) else math.exp(2.0 / kappa)
    return RPESchedule(M=M, inflation=damping * _overlap_inflation(overlap))


def round_circuit_toffolis(m: int, *, G_det: float, G_rand: float, L_D: int,
                           lam_R: float, delta: float, kappa: float,
                           n_stage: int = 2) -> float:
    """Toffolis in the single depth-``m`` circuit (no sample factor).

    The deterministic part applies ``n_stage*L_D`` rotations over ``2^(m-1)``
    halved Trotter steps; the randomized part draws ``kappa*lam_R^2*delta^2*4^m``
    Pauli rotations.  Each kind costs ``G_det`` and ``G_rand`` Toffolis.
    """
    det = G_det * n_stage * L_D * 2.0 ** (m - 1) if L_D > 0 else 0.0
    rand = (G_rand * kappa * lam_R ** 2 * delta ** 2 * 4.0 ** m
            if lam_R > 0.0 and delta > 0.0 else 0.0)
    return det + rand


def max_toffolis_per_circuit(schedule: RPESchedule, **circuit) -> float:
    """Toffolis in the deepest circuit, round ``m = M``."""
    return round_circuit_toffolis(schedule.M, **circuit)
