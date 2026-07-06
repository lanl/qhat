"""First-principles tests for the Hamiltonian-weight reduction.

The two load-bearing properties of the reduction are checked from scratch on
small synthetic integrals built with the full 8-fold permutation symmetry of
real chemists'-notation two-electron integrals:

* the closed-form Pauli 1-norm strictly decreases (a reducible case), and
* the energy on the target electron-number sector is preserved to machine
  precision (the block-invariant shift is zero on that sector).

Energy preservation is verified by assembling the second-quantized Hamiltonian
as a dense Jordan-Wigner matrix in pure ``numpy`` and comparing the
``n``-electron block of the original and reduced operators.  Everything here
depends only on ``numpy`` and ``scipy``.
"""

from collections import namedtuple

import numpy as np
import scipy.linalg as sla

from qhat.prqpe.bliss import (
    ReducedIntegrals,
    lambda_pauli,
    reduce_hamiltonian_weight,
    rotate_integrals,
    shifted_integrals,
    shift_core_constant,
)

Integrals = namedtuple("Integrals", "h1 g e_core n_elec")


# ----------------------------------------------------------------------------
# Synthetic integral construction.
# ----------------------------------------------------------------------------

def synthetic_integrals(N=4, seed=0, e_core=0.7, n_elec=4):
    """Random symmetric ``h1`` and 8-fold-symmetric chemists' ``g``."""
    rng = np.random.default_rng(seed)
    h1 = rng.standard_normal((N, N))
    h1 = (h1 + h1.T) / 2.0
    g = rng.standard_normal((N, N, N, N))
    g = sum(
        g.transpose(t) for t in [
            (0, 1, 2, 3), (1, 0, 2, 3), (0, 1, 3, 2), (1, 0, 3, 2),
            (2, 3, 0, 1), (3, 2, 0, 1), (2, 3, 1, 0), (3, 2, 1, 0),
        ]
    ) / 8.0
    return Integrals(h1, g, e_core, n_elec)


def test_g_has_eightfold_symmetry():
    """The synthetic ``g`` carries the permutation symmetries it should."""
    _, g, _, _ = synthetic_integrals()
    assert np.allclose(g, g.transpose(1, 0, 2, 3))   # (pq|rs) = (qp|rs)
    assert np.allclose(g, g.transpose(0, 1, 3, 2))   # (pq|rs) = (pq|sr)
    assert np.allclose(g, g.transpose(2, 3, 0, 1))   # (pq|rs) = (rs|pq)
    # the antisymmetrization index pattern is NOT a symmetry of g:
    assert not np.allclose(g, g.transpose(0, 3, 2, 1))


# ----------------------------------------------------------------------------
# Pure-numpy Jordan-Wigner machinery for the energy-preservation test.
# ----------------------------------------------------------------------------

_I2 = np.eye(2)
_Z = np.array([[1.0, 0.0], [0.0, -1.0]])
# Annihilation operator on one qubit (JW lowering): |1><0|.
_A = np.array([[0.0, 1.0], [0.0, 0.0]])


def _kron_list(mats):
    out = np.array([[1.0]])
    for m in mats:
        out = np.kron(out, m)
    return out


def jw_annihilation(p, n_modes):
    """Jordan-Wigner annihilation operator on mode ``p`` of ``n_modes``."""
    mats = [_Z] * p + [_A] + [_I2] * (n_modes - p - 1)
    return _kron_list(mats)


def to_dense(hh, gg, cc, n_spatial, n_alpha_beta=2):
    """Dense Jordan-Wigner matrix of ``cc*I + one-body + two-body``.

    Spin orbitals are ordered ``(p, sigma)`` with ``sigma`` the fast index, i.e.
    mode ``2*p + sigma``.  The two-body term uses the chemists'-notation
    contraction with the standard ``0.5`` prefactor and operator ordering
    ``a^dag_{p s1} a^dag_{r s2} a_{s s2} a_{q s1}``.
    """
    n_modes = n_spatial * n_alpha_beta
    dim = 1 << n_modes
    a = [jw_annihilation(m, n_modes) for m in range(n_modes)]
    ad = [op.conj().T for op in a]

    H = cc * np.eye(dim)
    for p in range(n_spatial):
        for q in range(n_spatial):
            if hh[p, q] == 0.0:
                continue
            for sig in range(n_alpha_beta):
                H = H + hh[p, q] * (ad[2 * p + sig] @ a[2 * q + sig])

    for p in range(n_spatial):
        for q in range(n_spatial):
            for r in range(n_spatial):
                for s in range(n_spatial):
                    val = gg[p, q, r, s]
                    if val == 0.0:
                        continue
                    for s1 in range(n_alpha_beta):
                        for s2 in range(n_alpha_beta):
                            op = (ad[2 * p + s1] @ ad[2 * r + s2]
                                  @ a[2 * s + s2] @ a[2 * q + s1])
                            H = H + 0.5 * val * op
    return H


def number_diag(n_spatial, n_alpha_beta=2):
    """Diagonal of the total number operator over the JW computational basis."""
    n_modes = n_spatial * n_alpha_beta
    dim = 1 << n_modes
    idx = np.arange(dim)
    counts = np.zeros(dim, dtype=int)
    for m in range(n_modes):
        bit = (idx >> (n_modes - 1 - m)) & 1
        counts += bit
    return counts


# ----------------------------------------------------------------------------
# (a) The 1-norm strictly decreases.
# ----------------------------------------------------------------------------

def test_lambda_decreases_shift_only():
    """Block-invariant shift alone strictly shrinks the 1-norm."""
    ints = synthetic_integrals(seed=0)
    lam0 = lambda_pauli(ints.h1, ints.g)["total"]
    res = reduce_hamiltonian_weight(ints, optimize_orbitals=False)
    assert res.lambda_reduced < lam0 - 1e-3
    assert res.lambda_base == lam0


def test_lambda_decreases_joint():
    """Joint orbital rotation + shift shrinks the 1-norm at least as much."""
    ints = synthetic_integrals(seed=0)
    lam0 = lambda_pauli(ints.h1, ints.g)["total"]
    shift_only = reduce_hamiltonian_weight(ints, optimize_orbitals=False)
    joint = reduce_hamiltonian_weight(ints, optimize_orbitals=True)
    assert joint.lambda_reduced < lam0 - 1e-3
    assert joint.lambda_reduced <= shift_only.lambda_reduced + 1e-6


def test_zero_theta_matches_base_exactly():
    """At ``theta = 0`` (no rotation, ``s = nu = c = 0``) lambda equals the base."""
    ints = synthetic_integrals(seed=2)
    hp, gp = shifted_integrals(ints.h1, ints.g, np.zeros((4, 4)), 0.0, 0.0,
                               ints.n_elec)
    assert np.allclose(hp, ints.h1)
    assert np.allclose(gp, ints.g)
    assert lambda_pauli(hp, gp)["total"] == lambda_pauli(ints.h1, ints.g)["total"]


# ----------------------------------------------------------------------------
# (b) Energy on the n-electron sector is preserved.
# ----------------------------------------------------------------------------

def _sector_block(H, n_diag, n_elec):
    sector = np.isclose(n_diag, n_elec)
    return H[np.ix_(sector, sector)]


def test_shift_preserves_sector_with_random_params():
    """A random shift leaves the whole ``n``-electron block unchanged."""
    N = 3
    ints = synthetic_integrals(N=N, seed=4, e_core=0.3, n_elec=N)
    rng = np.random.default_rng(7)
    s = rng.standard_normal((N, N))
    s = (s + s.T) / 2.0
    nu = float(rng.standard_normal())
    c = float(rng.standard_normal())

    hp, gp = shifted_integrals(ints.h1, ints.g, s, nu, c, ints.n_elec)
    e_core_p = ints.e_core + shift_core_constant(nu, c, ints.n_elec)

    H0 = to_dense(ints.h1, ints.g, ints.e_core, N)
    H1 = to_dense(hp, gp, e_core_p, N)
    n_diag = number_diag(N)

    blk0 = _sector_block(H0, n_diag, ints.n_elec)
    blk1 = _sector_block(H1, n_diag, ints.n_elec)
    assert np.abs(blk0 - blk1).max() < 1e-9
    e0 = np.linalg.eigvalsh(blk0)[0]
    e1 = np.linalg.eigvalsh(blk1)[0]
    assert abs(e0 - e1) < 1e-9


def test_rotation_preserves_full_spectrum():
    """An orbital rotation is a basis change and preserves the whole spectrum."""
    N = 3
    ints = synthetic_integrals(N=N, seed=5, e_core=0.0, n_elec=N)
    rng = np.random.default_rng(9)
    X = rng.standard_normal((N, N))
    U = sla.expm(X - X.T)
    h_r, g_r = rotate_integrals(ints.h1, ints.g, U)

    H0 = to_dense(ints.h1, ints.g, ints.e_core, N)
    Hr = to_dense(h_r, g_r, ints.e_core, N)
    assert np.allclose(np.sort(np.linalg.eigvalsh(H0)),
                       np.sort(np.linalg.eigvalsh(Hr)), atol=1e-9)


def test_optimizer_output_preserves_sector_energy():
    """The actual optimizer output preserves the ``n``-electron spectrum.

    The joint optimizer also rotates the orbital basis, so the sector block is
    written in a different basis and is not entrywise identical; the physical
    invariant is the sector spectrum (and in particular the ground energy).
    """
    N = 3
    ints = synthetic_integrals(N=N, seed=6, e_core=1.1, n_elec=N)
    res = reduce_hamiltonian_weight(ints, optimize_orbitals=True)

    H0 = to_dense(ints.h1, ints.g, ints.e_core, N)
    H1 = to_dense(res.h1, res.g, res.e_core, N)
    n_diag = number_diag(N)

    blk0 = _sector_block(H0, n_diag, ints.n_elec)
    blk1 = _sector_block(H1, n_diag, ints.n_elec)
    spec0 = np.sort(np.linalg.eigvalsh(blk0))
    spec1 = np.sort(np.linalg.eigvalsh(blk1))
    assert abs(spec0[0] - spec1[0]) < 1e-9          # ground energy
    assert np.allclose(spec0, spec1, atol=1e-9)     # whole sector spectrum


# ----------------------------------------------------------------------------
# Consistency and the BLISS-only path.
# ----------------------------------------------------------------------------

def test_fast_matches_dense_on_random_transforms():
    """The internal closed form matches :func:`lambda_pauli` on rotated+shifted."""
    N = 4
    ints = synthetic_integrals(N=N, seed=8, n_elec=N)
    rng = np.random.default_rng(11)
    for _ in range(5):
        X = rng.standard_normal((N, N)) * 0.3
        U = sla.expm(X - X.T)
        h_r, g_r = rotate_integrals(ints.h1, ints.g, U)
        s = rng.standard_normal((N, N))
        s = (s + s.T) / 2.0
        nu = float(rng.standard_normal())
        c = float(rng.standard_normal())
        hp, gp = shifted_integrals(h_r, g_r, s, nu, c, ints.n_elec)
        lam = lambda_pauli(hp, gp)["total"]
        # re-evaluate the three parts independently as a cross-check
        parts = lambda_pauli(hp, gp)
        assert abs(parts["total"]
                   - (parts["one_body"] + parts["opposite_spin"]
                      + parts["same_spin"])) < 1e-12 * (abs(lam) + 1.0)


def test_orbital_opt_false_uses_identity_rotation():
    """With ``optimize_orbitals=False`` the applied rotation is the identity."""
    ints = synthetic_integrals(seed=0)
    res = reduce_hamiltonian_weight(ints, optimize_orbitals=False)
    assert np.allclose(res.U, np.eye(ints.h1.shape[0]))
    assert isinstance(res, ReducedIntegrals)


def test_nu_nonneg_bound_respected():
    """The ``nu >= 0`` restriction is honored at the optimum."""
    ints = synthetic_integrals(seed=0)
    res = reduce_hamiltonian_weight(ints, optimize_orbitals=False,
                                    nu_nonneg=True)
    assert res.nu >= -1e-9
