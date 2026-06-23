"""Dense second-order Trotter unitaries and exact ground-state extraction.

The ground-state Trotter constant is read off the bias of the symmetric
second-order product formula.  This module assembles that formula as an explicit
matrix and pulls the effective ground-state energy out of it by diagonalization:

* :func:`trot2_matrix` -- the dense symmetric Trotter unitary ``S_2(delta)``;
* :func:`effective_ground_energy` -- the energy and overlap of the Trotter
  eigenstate closest to the exact ground state;
* :func:`ground_state` -- the exact ground state ``(E0, |psi0>)`` of the
  Hamiltonian.

Every routine here diagonalizes a dense ``dim x dim`` operator, so it is meant
for small systems (``dim`` up to a few thousand).  It depends only on ``numpy``,
``scipy``, and the :class:`~prqpe.representation.pauli.PauliHamiltonian`
container.
"""

from __future__ import annotations

import numpy as np

from ..representation.pauli import PauliHamiltonian


def trot2_matrix(ham: PauliHamiltonian, delta: float) -> np.ndarray:
    """The symmetric second-order Trotter unitary ``S_2(delta)`` as a dense array.

    For ``H = sum_l h_l P_l`` the symmetric (Strang) splitting is

        S_2(delta) = (prod_l exp(-i delta H_l / 2))
                     (prod_l reversed exp(-i delta H_l / 2)),

    a half step through the terms in stored order followed by a half step back.
    The product is built by left-multiplying the identity with each rotation via
    :meth:`PauliHamiltonian.rotation_apply_matrix`.

    Parameters
    ----------
    ham : PauliHamiltonian
        The Hamiltonian whose Trotter step is assembled.
    delta : float
        The Trotter step size.

    Returns
    -------
    numpy.ndarray
        The dense ``(dim, dim)`` unitary.  Building it costs ``O(L * dim^2)``
        work.
    """
    M = np.eye(ham.dim, dtype=complex)
    half = 0.5 * float(delta)
    # Rightmost factor first: ascending half-sweep, then descending half-sweep.
    for l in range(ham.L):
        M = ham.rotation_apply_matrix(l, half * ham.coeffs[l], M)
    for l in range(ham.L - 1, -1, -1):
        M = ham.rotation_apply_matrix(l, half * ham.coeffs[l], M)
    return M


def effective_ground_energy(U: np.ndarray, delta: float, gs) -> tuple[float, float]:
    """Effective ground-state energy of a Trotter unitary ``U = S_2(delta)``.

    Writing ``U = exp(-i delta H_eff)``, the eigenphases of ``U`` are
    ``delta`` times the effective energies.  The eigenvector with the largest
    overlap with the exact ground state ``gs`` is taken to be the Trotter image
    of the true ground state, and its energy is read from the eigenphase as
    ``E0_tilde = -arg(eigenvalue) / delta``.

    Parameters
    ----------
    U : numpy.ndarray
        A dense ``(dim, dim)`` Trotter unitary.
    delta : float
        The Trotter step that produced ``U``.
    gs : array_like
        The exact ground-state vector, length ``dim``.

    Returns
    -------
    tuple of float
        ``(E0_tilde, overlap_sq)``: the effective ground-state energy and the
        squared overlap of the selected eigenvector with ``gs``.
    """
    evals, evecs = np.linalg.eig(U)
    overlaps = np.abs(np.asarray(gs, dtype=complex).conj() @ evecs) ** 2
    j = int(np.argmax(overlaps))
    e0_tilde = -np.angle(evals[j]) / float(delta)
    return float(e0_tilde), float(overlaps[j])


def ground_state(ham: PauliHamiltonian) -> tuple[float, np.ndarray]:
    """The exact ground state ``(E0, |psi0>)`` of a Hamiltonian.

    Parameters
    ----------
    ham : PauliHamiltonian
        The Hamiltonian to diagonalize.

    Returns
    -------
    tuple
        ``(E0, psi0)``: the ground-state energy and the corresponding
        complex eigenvector.
    """
    evals, evecs = np.linalg.eigh(ham.to_dense())
    return float(evals[0]), np.ascontiguousarray(evecs[:, 0]).astype(complex)
