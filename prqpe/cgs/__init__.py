"""Ground-state Trotter constant ``C_gs``.

The cost core needs a single number, the ground-state Trotter constant ``C_gs``,
that sets how fast the second-order product-formula bias shrinks with the step
size: ``|E0 - E0_tilde(delta)| = C_gs * delta^p`` with ``p`` the formula order
(``p = 2`` for the symmetric Trotter step).  This module estimates that constant
two ways:

* :func:`estimate_cgs` -- the direct fit on a small system.  For a grid of step
  sizes it builds the exact Trotter unitary, extracts the effective ground-state
  energy, and fits the bias power law in log-log space.
* :func:`cgs_rule` -- the extrapolation ``C_gs = A * lambda^b`` for systems too
  large to diagonalize, with ``(A, b)`` fit offline from :func:`estimate_cgs`
  values across small systems.

The fit helper :func:`fit_cgs` is exposed for callers that already hold the bias
errors.  Everything here depends only on ``numpy``/``scipy`` and the
:class:`~prqpe.representation.pauli.PauliHamiltonian` container; the direct fit
is small-systems-only because it diagonalizes dense operators.
"""

from __future__ import annotations

import numpy as np

from ..representation.pauli import PauliHamiltonian
from .trotter_matrix import (
    effective_ground_energy,
    ground_state,
    trot2_matrix,
)

#: Default geometric grid of step sizes for :func:`estimate_cgs`.
DEFAULT_DELTAS = np.geomspace(0.2, 0.02, 7)

__all__ = [
    "estimate_cgs",
    "fit_cgs",
    "cgs_rule",
    "DEFAULT_DELTAS",
]


def fit_cgs(errors, deltas) -> tuple[float, float, np.ndarray]:
    """Fit the bias power law ``error = C_gs * delta^p`` in log-log space.

    A straight-line least-squares fit of ``log(error)`` against ``log(delta)``
    gives the exponent ``p`` (the slope) and ``C_gs`` (the exponentiated
    intercept).

    Parameters
    ----------
    errors : array_like
        The Trotter bias ``|E0 - E0_tilde(delta)|`` at each step size.  All
        entries must be strictly positive.
    deltas : array_like
        The corresponding step sizes.

    Returns
    -------
    tuple
        ``(C_gs, p, residuals)``: the fitted constant, the fitted exponent, and
        the fit residuals ``log(error) - (p * log(delta) + log(C_gs))``.
    """
    deltas = np.asarray(deltas, dtype=float)
    errors = np.asarray(errors, dtype=float)
    if np.any(errors <= 0.0):
        raise ValueError(
            "non-positive Trotter bias encountered; widen the step-size grid"
        )
    log_d = np.log(deltas)
    log_e = np.log(errors)
    slope, intercept = np.polyfit(log_d, log_e, 1)
    residuals = log_e - (slope * log_d + intercept)
    return float(np.exp(intercept)), float(slope), residuals


def estimate_cgs(ham: PauliHamiltonian, deltas=None) -> float:
    """Estimate the ground-state Trotter constant ``C_gs`` of a Hamiltonian.

    The exact ground state is found by diagonalization.  For each step size the
    symmetric Trotter unitary ``S_2(delta)`` is built, its effective
    ground-state energy is extracted, and the bias ``|E0 - E0_tilde(delta)|`` is
    fit to ``C_gs * delta^p`` in log-log space.  The returned value is the
    fitted ``C_gs``.

    Parameters
    ----------
    ham : PauliHamiltonian
        The Hamiltonian.
    deltas : array_like, optional
        The step-size grid.  Defaults to :data:`DEFAULT_DELTAS`, a geometric
        grid small enough to sit in the asymptotic regime.

    Returns
    -------
    float
        The fitted ground-state Trotter constant ``C_gs``.
    """
    deltas = DEFAULT_DELTAS if deltas is None else np.asarray(deltas, dtype=float)
    e0, psi0 = ground_state(ham)
    errors = np.empty(len(deltas), dtype=float)
    for i, delta in enumerate(deltas):
        U = trot2_matrix(ham, delta)
        e0_tilde, _ = effective_ground_energy(U, delta, psi0)
        errors[i] = abs(e0 - e0_tilde)
    c_gs, _p, _residuals = fit_cgs(errors, deltas)
    return c_gs


def cgs_rule(lam: float, A: float, b: float) -> float:
    """Extrapolate ``C_gs = A * lambda**b`` from a fitted scaling rule.

    ``(A, b)`` is fit offline from :func:`estimate_cgs` values on small systems;
    ``lam`` is the Hamiltonian 1-norm ``lambda = sum_l |h_l|``.
    """
    return float(A) * float(lam) ** float(b)
