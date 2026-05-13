"""
Shared test utilities for Pauli string operations.

This module provides common helper functions used across multiple test files.
"""

import numpy as np
from scipy.linalg import expm


def get_pauli_matrix(char):
    """Return the 2x2 Pauli matrix for a given character.

    Args:
        char: One of 'I', 'X', 'Y', 'Z'

    Returns:
        2x2 complex numpy array representing the Pauli matrix

    Raises:
        ValueError: If char is not a valid Pauli character
    """
    if char == 'I':
        return np.array([[1, 0], [0, 1]], dtype=complex)
    elif char == 'X':
        return np.array([[0, 1], [1, 0]], dtype=complex)
    elif char == 'Y':
        return np.array([[0, -1j], [1j, 0]], dtype=complex)
    elif char == 'Z':
        return np.array([[1, 0], [0, -1]], dtype=complex)
    else:
        raise ValueError(f"Invalid Pauli character: {char}")


def pauli_string_to_matrix(pauli_string):
    """Convert a Pauli string to its matrix representation via tensor product.

    Args:
        pauli_string: String of Pauli operators (e.g., "XYZ")

    Returns:
        2^n x 2^n complex numpy array, where n is the length of the string

    Example:
        >>> P = pauli_string_to_matrix("XY")
        >>> P.shape
        (4, 4)
    """
    result = get_pauli_matrix(pauli_string[0])
    for char in pauli_string[1:]:
        result = np.kron(result, get_pauli_matrix(char))
    return result


def analytical_evolution(pauli_string, coefficient, time, hbar=1.0):
    """Compute exp(-i * coefficient * PauliString * time / hbar) analytically.

    This computes the exact unitary evolution for a single Pauli string term
    using matrix exponentiation.

    Args:
        pauli_string: String of Pauli operators (e.g., "XYZ")
        coefficient: Scalar coefficient
        time: Evolution time
        hbar: Reduced Planck constant (default 1.0)

    Returns:
        2^n x 2^n complex numpy array representing the evolution operator

    Example:
        >>> U = analytical_evolution("X", 0.5, 1.0)
        >>> U.shape
        (2, 2)
    """
    P = pauli_string_to_matrix(pauli_string)
    return expm(-1j * coefficient * P * time / hbar)


def analytical_commuting_evolution(pauli_terms, time, hbar=1.0):
    """Compute the exact evolution for commuting terms analytically.

    Since terms commute: exp(-i*(A+B)*t) = exp(-i*A*t) * exp(-i*B*t)

    Args:
        pauli_terms: Sequence of (pauli_string, coefficient) tuples
        time: Evolution time
        hbar: Reduced Planck constant (default 1.0)

    Returns:
        2^n x 2^n complex numpy array representing the evolution operator

    Raises:
        ValueError: If pauli_terms is empty

    Example:
        >>> terms = [("XI", 0.5), ("IX", 0.3)]
        >>> U = analytical_commuting_evolution(terms, 1.0)
        >>> U.shape
        (4, 4)
    """
    if len(pauli_terms) == 0:
        raise ValueError("Must have at least one term")

    num_qubits = len(pauli_terms[0][0])
    dim = 2 ** num_qubits

    # Start with identity
    U = np.eye(dim, dtype=complex)

    # Multiply by each individual evolution
    for pauli_string, coefficient in pauli_terms:
        U_term = analytical_evolution(pauli_string, coefficient, time, hbar)
        U = U_term @ U

    return U
