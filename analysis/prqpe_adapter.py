"""QHAT -> prqpe adapter for analytic phase-estimation resource estimation.

This module marshals a QHAT ``Hamiltonian`` into the data the analytic ``prqpe``
estimator consumes (descending-sorted Pauli weights ``|h_l|``, the system qubit
count, and the ground-state Trotter constant ``C_gs``) and shapes the resulting
:class:`~qhat.prqpe.api.ResourceEstimate` into a flat, TOML-serializable dict.

The estimator itself implements the partially-randomized quantum phase
estimation cost model of Gunther et al., "Phase estimation with partially
randomized time evolution" (arXiv:2503.05647); this file is only the
QHAT-side glue and does not modify the prqpe core.
"""

from __future__ import annotations

import logging

import numpy as np
from openfermion import QubitOperator

from qhat.prqpe import estimate_resources
from qhat.prqpe.representation import build_representation, sorted_weights
from qhat.prqpe.cgs import estimate_cgs, cgs_rule

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------


def _qubit_operator_from_qhat(hamiltonian) -> QubitOperator:
    """Rebuild an openfermion ``QubitOperator`` from a QHAT Hamiltonian.

    The identity term (empty Pauli tuple ``()``) carries no resource cost and is
    skipped here; ``build_representation`` also drops it downstream.
    """
    qop = QubitOperator()
    for pauli_tuple, coef in hamiltonian.get_all_pauli_strings(return_as="tuples").items():
        if pauli_tuple == ():  # identity term -> classical constant, no cost
            continue
        qop += QubitOperator(pauli_tuple, coef)
    return qop


def _pick(obj, name, default):
    """``getattr`` with a fallback that also treats ``None`` as "use default"."""
    v = getattr(obj, name, None)
    return default if v is None else v


def _resolve_C_gs(pauli_ham, n_qubits, lam, config_analysis) -> float:
    """Resolve the ground-state Trotter constant ``C_gs`` (tiered strategy)."""
    explicit = getattr(config_analysis, "prqpe_C_gs", None)
    if explicit is not None:
        return float(explicit)

    rule = getattr(config_analysis, "prqpe_cgs_rule", None)
    if rule is not None:
        A, b = rule
        return float(cgs_rule(lam, A, b))

    max_qubits = getattr(config_analysis, "prqpe_cgs_max_qubits", 14)
    if n_qubits <= max_qubits:
        logger.info("Estimating C_gs by diagonalization (%d qubits)." % n_qubits)
        return float(estimate_cgs(pauli_ham))

    raise ValueError(
        f"C_gs required for {n_qubits} qubits (> prqpe_cgs_max_qubits={max_qubits}): "
        "the system is too large to diagonalize. Set analysis.prqpe_C_gs to an "
        "explicit value or analysis.prqpe_cgs_rule=(A, b) for the C_gs = A * lambda**b "
        "extrapolation."
    )


def _to_toml_primitive(v):
    """Coerce a metadata value to a TOML-serializable python primitive."""
    if isinstance(v, np.floating):
        return float(v)
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, bool):
        return v
    if isinstance(v, int):
        return v
    if isinstance(v, float):
        return v
    if isinstance(v, str):
        return v
    if isinstance(v, (np.complexfloating, complex)):
        return str(v)
    if isinstance(v, np.bool_):
        return bool(v)
    try:
        return float(v)
    except (TypeError, ValueError):
        return str(v)


# -------------------------------------------------------------------------------------------------


def resource_estimation_prqpe(config_analysis, algorithm, hamiltonian=None) -> dict:
    """Estimate phase-estimation resources for ``hamiltonian`` with prqpe."""
    if hamiltonian is None:
        raise ValueError("PRQPE resource estimation requires the Hamiltonian.")

    pauli_ham = build_representation(_qubit_operator_from_qhat(hamiltonian))
    weights = sorted_weights(pauli_ham)
    lam = float(weights.sum())
    n_qubits = hamiltonian.num_qubits()
    C_gs = _resolve_C_gs(pauli_ham, n_qubits, lam, config_analysis)

    est = estimate_resources(
        weights, n_qubits, C_gs,
        target_precision=_pick(config_analysis, "prqpe_target_precision", 0.0016),
        method=_pick(algorithm, "method", "partially_randomized"),
        overlap=_pick(algorithm, "overlap", 1.0),
        xi=getattr(algorithm, "xi", None),
        randomizer=_pick(algorithm, "randomizer", "rte"),
        commuting_group_size=getattr(algorithm, "commuting_group_size", None),
    )

    return {
        "estimator": "prqpe",
        "method": est.method,
        "toffoli_count": float(est.toffoli_count),
        "logical_qubits": int(est.logical_qubits),
        "max_toffoli_per_circuit": float(est.max_toffoli_per_circuit),
        "num_circuits": int(est.num_circuits),
        "C_gs": float(C_gs),
        "metadata": {k: _to_toml_primitive(v) for k, v in est.metadata.items()},
    }
