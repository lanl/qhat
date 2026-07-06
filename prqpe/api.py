"""Public entry point: estimate phase-estimation resources from term weights.

``estimate_resources`` takes the cost-relevant *data* about a Hamiltonian -- its
sorted Pauli weights ``|h_l|``, the system qubit count, and the ground-state
Trotter constant ``C_gs`` -- and returns a :class:`ResourceEstimate`.  Weight
reduction, the fermion-to-qubit mapping, and the ``C_gs`` measurement all happen
upstream.

The estimator sweeps the deterministic/randomized split ``L_D`` (the number of
largest-weight terms compiled deterministically), jointly optimizes the Trotter
step and randomization scale at each split, and keeps the whole sweep as a
Pareto curve alongside the cheapest point.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from . import costs, schedule

Method = Literal["deterministic", "randomized", "partially_randomized"]
Randomizer = Literal["rte", "qdrift"]

#: Default commuting-group length for Hamming-weight phasing.
DEFAULT_COMMUTING_GROUP_SIZE = 10

#: qDRIFT draws half as many rotations as RTE for the same accuracy.
_RANDOMIZER_FACTOR = {"rte": 1.0, "qdrift": 0.5}


# ----------------------------------------------------------------------------
# Result containers.
# ----------------------------------------------------------------------------

@dataclass(frozen=True)
class SplitEstimate:
    """Cost of one deterministic/randomized split -- one point on the sweep."""

    L_D: int
    lambda_R: float
    delta: float
    kappa: float
    eps_qpe: float
    toffoli_total: float
    max_toffoli_per_circuit: float
    num_circuits: int
    logical_qubits: int


@dataclass(frozen=True)
class ResourceEstimate:
    """The cheapest split, plus the full Pareto curve over ``L_D``."""

    toffoli_count: float
    logical_qubits: int
    max_toffoli_per_circuit: float
    num_circuits: int  # one guiding-state preparation per circuit
    method: str
    metadata: dict = field(default_factory=dict)
    pareto: list[SplitEstimate] = field(default_factory=list)


# ----------------------------------------------------------------------------
# Split selection.
# ----------------------------------------------------------------------------

def _candidate_splits(L: int, n_points: int = 200) -> np.ndarray:
    """Values of ``L_D`` to evaluate: log-spaced over ``[1, L]`` plus ``0`` and ``L``.

    The terms are sorted by descending weight, so ``L_D`` interpolates between
    the fully randomized (``L_D = 0``) and fully deterministic (``L_D = L``)
    endpoints.  A geometric spacing resolves the interior minimum without
    evaluating every integer split.
    """
    if L <= 1:
        return np.array([0, L][: L + 1], dtype=np.int64)
    interior = np.round(np.geomspace(1, L, n_points)).astype(np.int64)
    return np.unique(np.concatenate([[0], interior, [L]]))


# ----------------------------------------------------------------------------
# Per-split cost.
# ----------------------------------------------------------------------------

def _estimate_for_split(L_D: int, lam_R: float, n_system_qubits: int, lam: float,
                        C_gs: float, *, eps: float, order: int, K: int,
                        overlap: float, xi: float, randomizer: Randomizer,
                        eps_synth: float, n_stage: int = 2) -> tuple[SplitEstimate, dict]:
    """Optimize the cost at one split and return ``(SplitEstimate, metadata)``.

    The randomized part has weight ``lam_R`` and dyadic-angle exponent ``J``; the
    deterministic part has ``L_D`` terms whose synthesis cost ``G_det`` depends on
    the chosen Trotter step ``delta``, which in turn depends on ``G_det`` through
    the joint optimum -- a small fixed point resolved in a couple of iterations.
    """
    randomizer_factor = _RANDOMIZER_FACTOR[randomizer]

    # Randomized part: dyadic-angle exponent and per-rotation Toffoli cost.
    J = costs.phase_state_precision(lam_R if lam_R > 0.0 else lam, eps)
    G_rand = costs.rotation_toffolis(J, K)

    # Deterministic part: per-rotation synthesis cost.
    if L_D > 0:
        delta_seed = costs.trotter_step_size(eps, C_gs, order)
        eps_prime = costs.synthesis_precision_per_rotation(
            delta_seed, eps_synth, n_stage, L_D)
        G_det = costs.synthesis_toffolis(eps_prime)
        for _ in range(2):
            opt = costs.minimize_delta_kappa(
                L_D, lam_R, eps, C_gs, G_det, G_rand, order=order,
                n_stage=n_stage, randomizer_factor=randomizer_factor)
            if opt.delta > 0.0:
                eps_prime = costs.synthesis_precision_per_rotation(
                    opt.delta, eps_synth, n_stage, L_D)
                G_det = costs.synthesis_toffolis(eps_prime)
        opt = costs.minimize_delta_kappa(
            L_D, lam_R, eps, C_gs, G_det, G_rand, order=order,
            n_stage=n_stage, randomizer_factor=randomizer_factor)
    else:
        G_det = 0.0
        opt = costs.minimize_delta_kappa(
            0, lam_R, eps, C_gs, G_det, G_rand, order=order,
            n_stage=n_stage, randomizer_factor=randomizer_factor)

    # Round count and circuit schedule (at least one round for any input).
    if L_D == 0:
        M = schedule.num_rounds_randomized(lam_R, opt.eps_qpe, xi)
    else:
        M = schedule.num_rounds(opt.eps_qpe, opt.delta, xi)
    M = max(M, 0)
    sched = schedule.rpe_schedule(M, opt.kappa, overlap)

    # Deepest single circuit.  With no Trotter step the randomized arm reads it
    # off the closed-form total; otherwise it is the depth-M circuit explicitly.
    if L_D == 0:
        max_per_circuit = costs.randomized_deepest_circuit_toffolis(
            opt.toffoli_total)
    else:
        max_per_circuit = schedule.max_toffolis_per_circuit(
            sched, G_det=G_det, G_rand=G_rand * randomizer_factor, L_D=L_D,
            lam_R=lam_R, delta=opt.delta, kappa=opt.kappa, n_stage=n_stage)

    # Qubits: a randomized part needs the phasing register; a purely
    # deterministic arm only needs the system and the Hadamard-test ancilla.
    if lam_R > 0.0:
        qubits = costs.logical_qubits(n_system_qubits, J, K)
    else:
        qubits = n_system_qubits + 1

    estimate = SplitEstimate(
        L_D=int(L_D), lambda_R=float(lam_R), delta=opt.delta, kappa=opt.kappa,
        eps_qpe=opt.eps_qpe, toffoli_total=opt.toffoli_total,
        max_toffoli_per_circuit=float(max_per_circuit),
        num_circuits=sched.num_circuits, logical_qubits=int(qubits))
    meta = dict(M=M, J=J, K=K, G_det=G_det, G_rand=G_rand)
    return estimate, meta


# ----------------------------------------------------------------------------
# Entry point.
# ----------------------------------------------------------------------------

def estimate_resources(weights, n_system_qubits: int, C_gs: float, *,
                       target_precision: float = 0.0016,
                       method: Method = "partially_randomized",
                       representation: str = "pauli",
                       order: int = 2, overlap: float = 1.0,
                       xi: float | None = None,
                       randomizer: Randomizer = "rte",
                       synthesis_precision: float = 1e-4,
                       commuting_group_size: int | None = None) -> ResourceEstimate:
    """Estimate Toffoli and qubit costs for ground-state phase estimation.

    Parameters
    ----------
    weights : array of Pauli coefficients ``|h_l|`` (sorted internally, descending).
    n_system_qubits : qubits carrying the Hamiltonian.
    C_gs : ground-state Trotter constant, ``|E0 - E0~| = C_gs * delta^p``.
    target_precision : energy precision ``eps`` (Hartree).
    method : ``deterministic`` (L_D = L), ``randomized`` (L_D = 0), or
        ``partially_randomized`` (the full ``L_D`` sweep).
    representation : only ``pauli`` is implemented.
    overlap, xi : guiding-state overlap and depth-reduction factor.
    randomizer : ``rte`` (unbiased) or ``qdrift`` (half as many rotations).
    synthesis_precision : ground-state error budget for rotation synthesis.
    commuting_group_size : ``K`` for Hamming-weight phasing (default 10).

    Returns
    -------
    ResourceEstimate
        The cheapest split, with the full Pareto curve over ``L_D`` in ``pareto``.
    """
    if representation != "pauli":
        raise ValueError(
            f"unknown representation {representation!r}; only 'pauli' is implemented")
    if randomizer not in _RANDOMIZER_FACTOR:
        raise ValueError(f"unknown randomizer {randomizer!r}")

    w = np.sort(np.abs(np.asarray(weights, dtype=float)))[::-1]
    L = int(w.size)
    if L == 0:
        raise ValueError("weights is empty")
    lam = float(w.sum())
    K = commuting_group_size or DEFAULT_COMMUTING_GROUP_SIZE
    xi_value = schedule.resolve_depth_tradeoff(overlap, xi)

    if method == "deterministic":
        ld_values = np.array([L], dtype=np.int64)
    elif method == "randomized":
        ld_values = np.array([0], dtype=np.int64)
    elif method == "partially_randomized":
        ld_values = _candidate_splits(L)
    else:
        raise ValueError(f"unknown method {method!r}")

    # Tail weight: lambda_R(L_D) = sum of the weights not treated deterministically.
    tail = np.concatenate([np.cumsum(w[::-1])[::-1], [0.0]])

    results = [
        _estimate_for_split(
            int(L_D), float(tail[L_D]), n_system_qubits, lam, C_gs,
            eps=target_precision, order=order, K=K, overlap=overlap,
            xi=xi_value, randomizer=randomizer, eps_synth=synthesis_precision)
        for L_D in ld_values
    ]

    pareto = [est for est, _ in results]
    best_est, best_meta = min(results, key=lambda r: r[0].toffoli_total)

    metadata = dict(
        L_D=best_est.L_D, lambda_R=best_est.lambda_R, delta=best_est.delta,
        M=best_meta["M"], kappa=best_est.kappa, xi=xi_value,
        eps_qpe=best_est.eps_qpe, J=best_meta["J"], K=best_meta["K"],
        G_det=best_meta["G_det"], G_rand=best_meta["G_rand"])

    return ResourceEstimate(
        toffoli_count=best_est.toffoli_total,
        logical_qubits=best_est.logical_qubits,
        max_toffoli_per_circuit=best_est.max_toffoli_per_circuit,
        num_circuits=best_est.num_circuits,
        method=method, metadata=metadata, pareto=pareto)
