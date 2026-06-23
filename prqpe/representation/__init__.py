"""Hamiltonian-to-cost-input adapter.

The cost core consumes only the descending-sorted weights ``|h_l|``, the system
qubit count, and the ground-state Trotter constant.  This module owns the step that
turns a qubit Hamiltonian into the operator container that produces those weights:

* :func:`build_representation` -- convert a qubit Hamiltonian into a
  :class:`~prqpe.representation.pauli.PauliHamiltonian`.
* :func:`truncate_small_terms` -- drop the smallest-weight terms within a fixed
  weight budget.
* :func:`sorted_weights` -- the descending ``|h_l|`` array the cost core reads.

Only :func:`build_representation` touches an external qubit-Hamiltonian package,
and it imports it lazily.
"""

from __future__ import annotations

import numpy as np

from .pauli import PauliHamiltonian

#: Default cumulative weight budget for :func:`truncate_small_terms`.
DEFAULT_WEIGHT_THRESHOLD = 1e-3

__all__ = [
    "PauliHamiltonian",
    "build_representation",
    "truncate_small_terms",
    "sorted_weights",
    "DEFAULT_WEIGHT_THRESHOLD",
]


def build_representation(hamiltonian, representation: str = "pauli"
                         ) -> PauliHamiltonian:
    """Build a :class:`PauliHamiltonian` from a qubit Hamiltonian.

    Parameters
    ----------
    hamiltonian : PauliHamiltonian or qubit operator
        Either an already-built :class:`PauliHamiltonian` (returned unchanged) or
        an openfermion ``QubitOperator``.
    representation : str
        Only ``"pauli"`` is implemented.

    Returns
    -------
    PauliHamiltonian
        The Pauli-string Hamiltonian, with any identity term dropped (a constant
        energy shift carries no resource cost and is reinstated classically).
    """
    if representation != "pauli":
        raise ValueError(
            f"unknown representation {representation!r}; only 'pauli' is implemented")

    if isinstance(hamiltonian, PauliHamiltonian):
        return hamiltonian
    return _pauli_from_qubit_operator(hamiltonian)


def _pauli_from_qubit_operator(qop, tol: float = 1e-12) -> PauliHamiltonian:
    """Convert an openfermion ``QubitOperator`` to a :class:`PauliHamiltonian`.

    The letter-to-bitmask map follows :mod:`prqpe.representation.pauli`.  The
    identity term is dropped and terms with negligible real part are skipped.
    """
    from openfermion import count_qubits

    n_qubits = count_qubits(qop)
    xs, zs, cs = [], [], []
    for term, coeff in qop.terms.items():
        coeff = complex(coeff)
        if abs(coeff.imag) > 1e-9:
            raise ValueError(f"non-real Pauli coefficient {coeff} for {term}")
        if not term:  # identity term -> classical constant, no cost
            continue
        if abs(coeff.real) < tol:
            continue
        x = z = 0
        for q, letter in term:
            bit = 1 << (n_qubits - 1 - q)
            if letter in ("X", "Y"):
                x |= bit
            if letter in ("Z", "Y"):
                z |= bit
        xs.append(x)
        zs.append(z)
        cs.append(coeff.real)
    return PauliHamiltonian(n_qubits, xs, zs, cs)


def truncate_small_terms(ham: PauliHamiltonian,
                         weight_threshold: float = DEFAULT_WEIGHT_THRESHOLD
                         ) -> PauliHamiltonian:
    """Drop the smallest-weight terms within a cumulative weight budget.

    Working upward from the smallest ``|h_l|``, terms are dropped while the
    *cumulative* dropped weight stays within ``weight_threshold``.  The first
    term that would push the cumulative dropped weight past the budget is kept,
    along with everything larger; the kept Hamiltonian therefore differs from the
    input by at most ``weight_threshold`` in 1-norm.

    Parameters
    ----------
    ham : PauliHamiltonian
        The Hamiltonian to truncate.
    weight_threshold : float
        Maximum total dropped weight (same units as the coefficients).

    Returns
    -------
    PauliHamiltonian
        A new Hamiltonian with the kept terms, ordered by decreasing ``|h_l|``.
    """
    ordered = ham.sorted_by_weight()
    weights = np.abs(ordered.coeffs)
    # Cumulative weight accumulated from the smallest term upward.
    cum_from_small = np.cumsum(weights[::-1])
    # Number of trailing (smallest) terms whose cumulative weight stays within
    # the budget: searchsorted from the right finds the largest such count.
    n_drop = int(np.searchsorted(cum_from_small, weight_threshold, side="right"))
    n_keep = ordered.L - n_drop
    return ordered.subset(np.arange(n_keep))


def sorted_weights(ham: PauliHamiltonian) -> np.ndarray:
    """The descending ``|h_l|`` array the cost core consumes.

    Equivalent to ``np.abs(ham.sorted_by_weight().coeffs)``.
    """
    return np.abs(ham.sorted_by_weight().coeffs)
