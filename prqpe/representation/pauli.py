"""Real-coefficient Pauli-string Hamiltonian container.

A Hamiltonian is stored as a weighted sum of Pauli strings,

    H = sum_l h_l P_l,   h_l real,   P_l a tensor product of {I, X, Y, Z},

with 1-norm ``lambda = sum_l |h_l|``.

Each Pauli string is held as a pair of integer bitmasks ``(x, z)`` with the phase
convention

    P(x, z) = i^popcount(x AND z) * X^x Z^z

(the X factor applied after the Z factor).  This is Hermitian and self-inverse, and
maps the letter ``Y`` on a qubit to ``x_q = z_q = 1`` since ``Y = i X Z``.  Bit
``n_qubits - 1 - q`` of a computational-basis index belongs to qubit ``q`` -- qubit 0
is the most significant bit -- so the dense and sparse matrices built here agree
entrywise with the standard sparse-operator ordering used by qubit-Hamiltonian
toolkits.

Action on a basis state ``|b>``::

    P(x, z) |b> = i^popcount(x & z) * (-1)^popcount(b & z) |b XOR x>,

hence ``(P v)[c] = i^popcount(x & z) * (-1)^popcount((c ^ x) & z) * v[c ^ x]``: a bit
permutation followed by a diagonal phase, vectorized over terms below.

Pauli rotations use ``P^2 = I``::

    exp(-i phi P) |v> = cos(phi) |v> - i sin(phi) P |v>.

This module depends only on ``numpy`` and ``scipy``.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def popcount(a):
    """Number of set bits, elementwise, for a nonnegative integer array."""
    a = np.asarray(a, dtype=np.uint64)
    if hasattr(np, "bitwise_count"):
        return np.bitwise_count(a).astype(np.int64)
    out = np.zeros(a.shape, dtype=np.int64)
    x = a.copy()
    while np.any(x):
        out += (x & np.uint64(1)).astype(np.int64)
        x >>= np.uint64(1)
    return out


class PauliHamiltonian:
    """A real-coefficient Pauli-string Hamiltonian ``H = sum_l h_l P_l``.

    Parameters
    ----------
    n_qubits : int
        Number of qubits the Hamiltonian acts on.
    x_masks, z_masks : integer arrays of shape ``(L,)``
        The X and Z bitmasks of each Pauli string, in the convention of the
        module docstring.
    coeffs : float array of shape ``(L,)``
        The real coefficients ``h_l``.
    """

    def __init__(self, n_qubits, x_masks, z_masks, coeffs):
        self.n_qubits = int(n_qubits)
        self.dim = 1 << self.n_qubits
        self.x_masks = np.asarray(x_masks, dtype=np.int64).ravel()
        self.z_masks = np.asarray(z_masks, dtype=np.int64).ravel()
        coeffs = np.asarray(coeffs)
        if np.iscomplexobj(coeffs):
            if np.max(np.abs(coeffs.imag), initial=0.0) > 1e-10:
                raise ValueError("Pauli coefficients must be real")
            coeffs = coeffs.real
        self.coeffs = np.asarray(coeffs, dtype=np.float64).ravel()
        if not (len(self.x_masks) == len(self.z_masks) == len(self.coeffs)):
            raise ValueError("mask/coefficient length mismatch")
        self._perms = None
        self._phases = None

    # ------------------------------------------------------------------ basic
    @property
    def L(self):
        """Number of Pauli terms."""
        return len(self.coeffs)

    @property
    def lambda_weight(self):
        """The 1-norm ``lambda = sum_l |h_l|``."""
        return float(np.abs(self.coeffs).sum())

    def subset(self, idx):
        """A new Hamiltonian holding the terms ``idx`` (order preserved)."""
        idx = np.asarray(idx)
        return PauliHamiltonian(
            self.n_qubits, self.x_masks[idx], self.z_masks[idx], self.coeffs[idx]
        )

    def sorted_by_weight(self):
        """A copy with terms ordered by decreasing ``|h_l|`` (stable sort)."""
        order = np.argsort(-np.abs(self.coeffs), kind="stable")
        return self.subset(order)

    # ------------------------------------------------------------------ tables
    def _tables(self):
        """Per-term permutation and diagonal phase.

        Returns ``(perms, phases)`` of shape ``(L, dim)`` such that
        ``(P_l v)[c] = phases[l, c] * v[perms[l, c]]``: ``perms`` is the
        column permutation ``c XOR x_l`` and ``phases`` carries the head phase
        ``i^popcount(x_l & z_l)`` times the basis sign ``(-1)^popcount(c & z_l)``.
        Cached after the first call.
        """
        if self._perms is None:
            c = np.arange(self.dim, dtype=np.int64)
            self._perms = c[None, :] ^ self.x_masks[:, None]
            head = 1j ** (popcount(self.x_masks & self.z_masks) % 4)
            signs = 1.0 - 2.0 * (popcount(self._perms & self.z_masks[:, None]) % 2)
            self._phases = head[:, None] * signs
        return self._perms, self._phases

    # -------------------------------------------------------- matrix operators
    def pauli_apply_matrix(self, l, M):
        """Left-multiply a dense ``(dim, dim)`` matrix: ``P_l @ M``."""
        perms, phases = self._tables()
        return phases[l][:, None] * M[perms[l], :]

    def rotation_apply_matrix(self, l, phi, M):
        """``exp(-i phi P_l) @ M`` for dense ``M`` (assembles product-formula unitaries)."""
        return np.cos(phi) * M - 1j * np.sin(phi) * self.pauli_apply_matrix(l, M)

    def term_matrix(self, l):
        """The sparse matrix of ``P_l``."""
        perms, phases = self._tables()
        data = phases[l].copy()
        rows = np.arange(self.dim)
        return sp.csr_matrix((data, (rows, perms[l])), shape=(self.dim, self.dim))

    def to_sparse(self):
        """``H = sum_l h_l P_l`` as a scipy CSR matrix (small ``n_qubits`` only)."""
        perms, phases = self._tables()
        data = (self.coeffs[:, None] * phases).ravel()
        rows = np.tile(np.arange(self.dim), self.L)
        cols = perms.ravel()
        return sp.coo_matrix(
            (data, (rows, cols)), shape=(self.dim, self.dim)
        ).tocsr()

    def to_dense(self):
        """``H`` as a dense complex array (small ``n_qubits`` only)."""
        return self.to_sparse().toarray()
