"""Hamiltonian-weight reduction for ground-state phase estimation.

The gate cost scales with the Hamiltonian 1-norm ``lambda = sum_l |h_l|``, which
for an electronic Hamiltonian in spatial-orbital chemists'-notation integrals
``(h1, g)`` has a closed form (see :func:`lambda_pauli`).  Two freedoms shrink it
without changing the physics on the target ``n``-electron sector:

* an **orbital rotation** ``U`` (a basis change, so it preserves the spectrum), and
* a **block-invariant shift** by the number-operator deficit ``(n_hat - n)``,
  which is identically zero on the ``n``-electron sector but cancels 1-norm on
  the other particle-number sectors.

The entry point :func:`reduce_hamiltonian_weight` jointly minimizes the
closed-form 1-norm over both freedoms with an annealed smoothing of ``|.|`` and
L-BFGS-B.  The shift follows Loaiza and Izmaylov, Quantum Sci. Technol. 8,
035019 (2023).  Depends only on ``numpy`` and ``scipy``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.linalg as sla
from scipy.optimize import minimize

__all__ = ["ReducedIntegrals", "reduce_hamiltonian_weight"]

#: Default smoothing-parameter schedule for the annealed objective.  Each stage
#: replaces ``|x|`` by ``sqrt(x^2 + mu^2) - mu`` and warm-starts from the
#: previous stage; the final ``mu = 0`` stage is the exact (nonsmooth) 1-norm.
DEFAULT_MU_SCHEDULE = (1e-1, 1e-2, 1e-3, 1e-5, 0.0)


# ----------------------------------------------------------------------------
# Closed-form Pauli 1-norm.
# ----------------------------------------------------------------------------

def effective_one_body(h1: np.ndarray, g: np.ndarray) -> np.ndarray:
    """The effective one-body matrix ``T`` folding the two-body mean field.

    With spatial-orbital chemists'-notation integrals ``g = (pq|rs)`` the
    one-body part of the normal-ordered Hamiltonian is

        T = h1 + sum_r g[:, :, r, r] - 0.5 * sum_r g[:, r, r, :].
    """
    return (h1
            + np.einsum("pqrr->pq", g)
            - 0.5 * np.einsum("prrq->pq", g))


def lambda_pauli(h1: np.ndarray, g: np.ndarray) -> dict:
    """Closed-form Jordan-Wigner / Bravyi-Kitaev 1-norm of ``(h1, g)``.

    Inputs are spatial-orbital, real, chemists'-notation integrals: ``h1`` is the
    symmetric bare one-body matrix ``(N, N)`` and ``g = (pq|rs)`` is the
    8-fold-symmetric two-electron tensor ``(N, N, N, N)``.  The 1-norm splits as

        lambda = lambda_T + lambda_V2 + lambda_V3,

    a one-body term ``lambda_T = sum |T|`` with ``T`` from
    :func:`effective_one_body`, an opposite-spin two-body term
    ``lambda_V2 = sum |g| / 4``, and a same-spin two-body term
    ``lambda_V3 = sum_{p!=r, q!=s} |D| / 8`` with the antisymmetrized
    ``D[p, q, r, s] = g[p, q, r, s] - g[p, s, r, q]``.  Returns a dict with the
    three parts and their ``total``.
    """
    N = h1.shape[0]
    T = effective_one_body(h1, g)
    D = g - g.transpose(0, 3, 2, 1)
    off = ~np.eye(N, dtype=bool)
    mask = off[:, None, :, None] & off[None, :, None, :]  # p != r AND q != s
    lam_T = float(np.abs(T).sum())
    lam_V2 = float(np.abs(g).sum() / 4.0)
    lam_V3 = float(np.abs(D)[mask].sum() / 8.0)
    return {
        "one_body": lam_T,
        "opposite_spin": lam_V2,
        "same_spin": lam_V3,
        "total": lam_T + lam_V2 + lam_V3,
    }


# ----------------------------------------------------------------------------
# Integral transforms: orbital rotation and block-invariant shift.
# ----------------------------------------------------------------------------

def rotate_integrals(h1: np.ndarray, g: np.ndarray,
                     U: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotate ``(h1, g)`` into the orbital basis defined by orthogonal ``U``.

    ``U`` mixes spatial orbitals; the change of basis is

        h1_rot = U.T @ h1 @ U,
        g_rot[a, b, c, d] = sum g[p, q, r, s] U[p, a] U[q, b] U[r, c] U[s, d].

    The four-index contraction is staged as two einsums for cost.  An orbital
    rotation is a unitary change of basis and preserves the full spectrum, so the
    core energy is unchanged (handled by the caller).
    """
    h1_rot = U.T @ h1 @ U
    g_rot = np.einsum("pqrs,pa,qb->abrs", g, U, U, optimize=True)
    g_rot = np.einsum("abrs,rc,sd->abcd", g_rot, U, U, optimize=True)
    return h1_rot, g_rot


def shifted_integrals(h1: np.ndarray, g: np.ndarray, s: np.ndarray,
                      nu: float, c: float,
                      n_elec: float) -> tuple[np.ndarray, np.ndarray]:
    """Block-invariant shift of ``(h1, g)`` by symmetry-operator parameters.

    The shift adds ``Delta H`` built from the number-operator deficit
    ``(n_hat - n)`` with a symmetric one-body operator ``s`` and scalars
    ``nu, c``.  Expanded into integral coefficients (``I`` the ``N x N`` identity,
    ``n = n_elec``):

        gp = g + einsum('pq,rs->pqrs', s, I) + einsum('pq,rs->pqrs', I, s)
               + nu * einsum('pq,rs->pqrs', I, I),
        hp = h1 + (1 - n) * s + (nu * (0.5 - n) + c) * I.

    Returns ``(hp, gp)`` only; the accompanying core-energy constant
    ``nu * n^2 / 2 - c * n`` is added by the caller.
    """
    N = h1.shape[0]
    eye = np.eye(N)
    gp = (g
          + np.einsum("pq,rs->pqrs", s, eye)
          + np.einsum("pq,rs->pqrs", eye, s)
          + nu * np.einsum("pq,rs->pqrs", eye, eye))
    hp = h1 + (1.0 - n_elec) * s + (nu * (0.5 - n_elec) + c) * eye
    return hp, gp


def shift_core_constant(nu: float, c: float, n_elec: float) -> float:
    """The core-energy constant added by the shift: ``nu * n^2 / 2 - c * n``."""
    return nu * n_elec ** 2 / 2.0 - c * n_elec


# ----------------------------------------------------------------------------
# Smoothed objective and analytic gradient.
# ----------------------------------------------------------------------------

def _sabs(x: np.ndarray, mu: float) -> np.ndarray:
    """Smoothed absolute value ``sqrt(x^2 + mu^2) - mu`` (equals ``|x|`` at mu=0)."""
    if mu == 0.0:
        return np.abs(x)
    return np.sqrt(x * x + mu * mu) - mu


def _sgrad(x: np.ndarray, mu: float) -> np.ndarray:
    """Derivative of :func:`_sabs`: ``x / sqrt(x^2 + mu^2)`` (sign at mu=0)."""
    if mu == 0.0:
        return np.sign(x)
    return x / np.sqrt(x * x + mu * mu)


def _lambda_smooth_value(hp: np.ndarray, gp: np.ndarray, mu: float,
                         mask: np.ndarray) -> float:
    """Smoothed 1-norm on shifted integrals ``(hp, gp)`` at smoothing ``mu``."""
    T = effective_one_body(hp, gp)
    D = gp - gp.transpose(0, 3, 2, 1)
    lam_T = _sabs(T, mu).sum()
    lam_V2 = 0.25 * _sabs(gp, mu).sum()
    lam_V3 = 0.125 * (_sabs(D, mu) * mask).sum()
    return float(lam_T + lam_V2 + lam_V3)


def _lambda_smooth_grad_integrals(hp: np.ndarray, gp: np.ndarray, mu: float,
                                  mask: np.ndarray
                                  ) -> tuple[np.ndarray, np.ndarray]:
    """Cotangents ``(dL/dhp, dL/dgp)`` of the smoothed 1-norm on ``(hp, gp)``."""
    N = hp.shape[0]
    eye = np.eye(N)
    T = effective_one_body(hp, gp)
    D = gp - gp.transpose(0, 3, 2, 1)

    dT = _sgrad(T, mu)                       # dL/dT
    dgp = 0.25 * _sgrad(gp, mu)             # opposite-spin term contribution

    # one-body term backprop: T = h1 + einsum('pqrr->pq') - 0.5 einsum('prrq->pq')
    dhp = dT.copy()
    dgp += np.einsum("pq,rs->pqrs", dT, eye)          # d g[p, q, r, r] += dT[p, q]
    dgp += -0.5 * np.einsum("ps,qr->pqrs", dT, eye)   # d g[p, r, r, q] += -0.5 dT[p, q]

    # same-spin term backprop: D = gp - gp.transpose(0, 3, 2, 1)
    dD = 0.125 * (_sgrad(D, mu) * mask)
    dgp += dD - dD.transpose(0, 3, 2, 1)
    return dhp, dgp


# ----------------------------------------------------------------------------
# Parameter packing.
# ----------------------------------------------------------------------------

@dataclass(frozen=True)
class _Layout:
    """Flat-vector layout for the optimization parameters ``theta``.

    ``theta = concat([X_lower (if orbital), s_upper, [nu], [c]])`` with
    ``X_lower`` the strict lower triangle of the antisymmetric rotation generator
    and ``s_upper`` the upper triangle (incl. diagonal) of the symmetric shift
    matrix.
    """

    N: int
    optimize_orbitals: bool

    @property
    def itril(self):
        return np.tril_indices(self.N, -1)

    @property
    def iu(self):
        return np.triu_indices(self.N)

    @property
    def nX(self):
        return self.N * (self.N - 1) // 2 if self.optimize_orbitals else 0

    @property
    def nS(self):
        return self.N * (self.N + 1) // 2

    @property
    def n_par(self):
        return self.nX + self.nS + 2

    def unpack(self, theta):
        """``theta -> (X, s, nu, c)`` with antisymmetric ``X`` and symmetric ``s``."""
        N = self.N
        theta = np.asarray(theta)
        X = np.zeros((N, N), dtype=theta.dtype)
        off = 0
        if self.optimize_orbitals:
            X[self.itril] = theta[:self.nX]
            X = X - X.T
            off = self.nX
        s = np.zeros((N, N), dtype=theta.dtype)
        s[self.iu] = theta[off:off + self.nS]
        s = s + s.T - np.diag(np.diag(s))
        nu = theta[off + self.nS]
        c = theta[off + self.nS + 1]
        return X, s, nu, c


# ----------------------------------------------------------------------------
# The full pipeline (rotate then shift) and its gradient.
# ----------------------------------------------------------------------------

def _pipeline(theta, layout, h1, g, n_elec):
    """Apply rotation then shift; return ``(hp, gp, U, X, s, nu, c)``."""
    N = layout.N
    X, s, nu, c = layout.unpack(theta)
    if layout.optimize_orbitals:
        U = sla.expm(X)  # X is already the antisymmetric generator (X - X^T)
        h_r, g_r = rotate_integrals(h1, g, U)
    else:
        U = np.eye(N)
        h_r, g_r = h1, g
    hp, gp = shifted_integrals(h_r, g_r, s, nu, c, n_elec)
    return hp, gp, U, X, s, nu, c


def _shift_block_grad(dhp, dgp, layout, n_elec):
    """Backprop integral cotangents through the shift to ``(d s_upper, d nu, d c)``."""
    N = layout.N
    eye = np.eye(N)
    # hp = h_r + (1 - n) s + (nu (0.5 - n) + c) I
    # gp = g_r + s (x) I + I (x) s + nu I (x) I
    ds = ((1.0 - n_elec) * dhp
          + np.einsum("pqrs,rs->pq", dgp, eye)
          + np.einsum("pqrs,pq->rs", dgp, eye))
    # fold the lower triangle onto the upper (diagonal once) to match the symmetric parametrization
    ds_sym = ds + ds.T - np.diag(np.diag(ds))
    d_s_upper = ds_sym[layout.iu]
    d_nu = float((0.5 - n_elec) * np.trace(dhp)
                 + np.einsum("pqrs,pq,rs->", dgp, eye, eye))
    d_c = float(np.trace(dhp))
    return d_s_upper, d_nu, d_c


def _make_funcs(layout, h1, g, n_elec, mask):
    """Build ``(value, grad)`` closures for the smoothed objective at a given mu.

    The shift-block gradient is analytic.  The orbital-block gradient (the
    expm pushforward) is obtained by a complex-step derivative on the rotation
    parameters only, exact to machine precision.
    """
    N = layout.N

    def value(theta, mu):
        hp, gp, *_ = _pipeline(theta, layout, h1, g, n_elec)
        return _lambda_smooth_value(hp, gp, mu, mask)

    def grad(theta, mu):
        theta = np.asarray(theta, dtype=float)
        hp, gp, U, X, s, nu, c = _pipeline(theta, layout, h1, g, n_elec)
        dhp, dgp = _lambda_smooth_grad_integrals(hp, gp, mu, mask)
        d_s_upper, d_nu, d_c = _shift_block_grad(dhp, dgp, layout, n_elec)

        if layout.optimize_orbitals:
            d_X = np.zeros(layout.nX)
            step = 1e-30
            base = theta.copy()
            for k in range(layout.nX):
                pert = base.astype(complex)
                pert[k] += 1j * step
                val = _value_complex(pert, layout, h1, g, n_elec, mu, mask)
                d_X[k] = val.imag / step
            return np.concatenate([d_X, d_s_upper, [d_nu], [d_c]])
        return np.concatenate([d_s_upper, [d_nu], [d_c]])

    return value, grad


def _value_complex(theta, layout, h1, g, n_elec, mu, mask):
    """Complex-step-friendly smoothed value (holomorphic sqrt, not np.abs).

    Used only for the rotation-block complex-step gradient, so mu is the
    positive smoothing of the current anneal stage and sqrt is holomorphic.
    """
    hp, gp, *_ = _pipeline(theta, layout, h1, g, n_elec)
    T = effective_one_body(hp, gp)
    D = gp - gp.transpose(0, 3, 2, 1)
    sab = lambda z: np.sqrt(z * z + mu * mu) - mu
    return sab(T).sum() + 0.25 * sab(gp).sum() + 0.125 * (sab(D) * mask).sum()


# ----------------------------------------------------------------------------
# Result container.
# ----------------------------------------------------------------------------

@dataclass(frozen=True)
class ReducedIntegrals:
    """Reduced integrals and provenance from :func:`reduce_hamiltonian_weight`."""

    h1: np.ndarray             # (N, N) reduced bare one-body integrals (rotated + shifted)
    g: np.ndarray              # (N, N, N, N) reduced chemists' (pq|rs), 8-fold symmetric
    e_core: float              # core energy incl. shift constant
    n_elec: int                # target electron number (the preserved sector)
    lambda_reduced: float      # final closed-form Pauli 1-norm at the optimum
    lambda_base: float         # 1-norm of the input integrals (no rotation/shift)
    U: np.ndarray = field(default=None)         # orbital rotation actually applied
    s: np.ndarray = field(default=None)         # symmetric shift matrix at optimum
    nu: float = 0.0                             # quadratic number shift at optimum
    c: float = 0.0                              # linear number shift at optimum


# ----------------------------------------------------------------------------
# Public entry point.
# ----------------------------------------------------------------------------

def reduce_hamiltonian_weight(
    integrals,
    *,
    optimize_orbitals: bool = True,
    nu_nonneg: bool = False,
    mu_schedule=DEFAULT_MU_SCHEDULE,
    maxiter: int = 3000,
    x0: np.ndarray = None,
    rng_seed: int = 0,
) -> ReducedIntegrals:
    """Minimize the closed-form Pauli 1-norm over orbital and shift freedoms.

    ``integrals`` exposes ``.h1`` ``(N, N)``, ``.g`` ``(N, N, N, N)``, ``.e_core``
    and ``.n_elec``.  The 1-norm of :func:`lambda_pauli` is minimized over an
    orbital rotation ``U = expm(X - X^T)`` and a block-invariant shift
    ``(s, nu, c)`` on the rotated-then-shifted integrals; the nonsmooth ``|.|`` is
    annealed through ``mu_schedule`` with L-BFGS-B (warm-started between stages),
    the final stage being the exact 1-norm.  The returned :class:`ReducedIntegrals`
    carries the reduced spatial-orbital integrals and the shifted core energy
    ``e_core + nu * n^2 / 2 - c * n``; mapping them to qubits is the caller's job.

    Parameters
    ----------
    integrals : object with ``.h1``, ``.g``, ``.e_core``, ``.n_elec``.
    optimize_orbitals : include the orbital rotation (else fit the shift only).
    nu_nonneg : constrain the quadratic number-shift scalar ``nu >= 0``.
    mu_schedule : decreasing smoothing parameters; the final value should be ``0``.
    maxiter : maximum L-BFGS-B iterations per annealing stage.
    x0 : explicit flat start vector of length ``n_par`` (overrides ``rng_seed``).
    rng_seed : seed for the small random start perturbation; distinct seeds reach
        distinct local optima of the nonconvex objective (``0`` = zero start).
    """
    h1 = np.asarray(integrals.h1, dtype=float)
    g = np.asarray(integrals.g, dtype=float)
    e_core = float(integrals.e_core)
    n_elec = int(integrals.n_elec)
    N = h1.shape[0]

    layout = _Layout(N=N, optimize_orbitals=optimize_orbitals)
    off = ~np.eye(N, dtype=bool)
    mask = off[:, None, :, None] & off[None, :, None, :]

    lambda_base = lambda_pauli(h1, g)["total"]
    value, grad = _make_funcs(layout, h1, g, n_elec, mask)

    def fun(theta, mu):
        return value(theta, mu), grad(theta, mu)

    if x0 is not None:
        x = np.asarray(x0, dtype=float).copy()
        if x.shape[0] != layout.n_par:
            raise ValueError("x0 has the wrong length for the parameter layout")
    elif rng_seed:
        # A small random perturbation of the zero start; the scale is tiny so the
        # first eval is still near the base 1-norm while distinct seeds reach
        # distinct local optima of the nonconvex objective.
        rng = np.random.default_rng(rng_seed)
        x = 1e-2 * rng.standard_normal(layout.n_par)
    else:
        x = np.zeros(layout.n_par)

    bounds = None
    if nu_nonneg:
        bounds = [(None, None)] * layout.n_par
        bounds[layout.nX + layout.nS] = (0.0, None)  # nu slot

    res = None
    for mu in mu_schedule:
        res = minimize(
            fun, x, args=(mu,), jac=True, method="L-BFGS-B", bounds=bounds,
            options=dict(maxiter=maxiter, maxfun=10 * maxiter,
                         ftol=1e-15, gtol=1e-12, maxcor=30),
        )
        x = res.x

    # The achieved 1-norm at the optimum.
    hp, gp, U, X, s, nu, c = _pipeline(x, layout, h1, g, n_elec)
    lambda_reduced = lambda_pauli(hp, gp)["total"]

    return ReducedIntegrals(
        h1=hp,
        g=gp,
        e_core=e_core + shift_core_constant(nu, c, n_elec),
        n_elec=n_elec,
        lambda_reduced=float(lambda_reduced),
        lambda_base=float(lambda_base),
        U=U,
        s=s,
        nu=float(nu),
        c=float(c),
    )
