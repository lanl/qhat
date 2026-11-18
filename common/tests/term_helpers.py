"""Shared helpers for building Hamiltonian terms in tests."""

from qhat.common.terms.single_pauli_term import SinglePauliTerm


def single_pauli_terms(pauli_terms):
    """Convert a sequence of (pauli_string, coefficient) tuples to a tuple of SinglePauliTerm."""
    return tuple(SinglePauliTerm(pauli, coefficient) for pauli, coefficient in pauli_terms)
