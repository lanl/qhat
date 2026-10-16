"""Dense product-formula Trotter diagnostics.

These routines are intentionally small wrappers around QHAT's existing
Trotter coefficient expansion and Pauli matrix utilities.  They are useful for
numerical JW/BK comparisons on small systems where dense matrices are practical.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
from scipy.linalg import expm, norm

from .trotter_sequences import (
    expand_ramped_trotterization,
    get_trotterization_coefficients,
)
from .pauli_keys import pauli_key_to_matrix

TermList = Sequence[tuple[tuple[tuple[int, str], ...], complex]]


def exact_unitary(hamiltonian: np.ndarray, time: float) -> np.ndarray:
    """Compute exact time evolution ``exp(-i time H)``."""
    return expm(-1j * time * hamiltonian)


def dense_trotter_unitary(
    term_list: TermList,
    time: float,
    num_steps: int,
    n_qubits: int,
    method=1,
    combine_terms: bool = True,
) -> np.ndarray:
    """Build a dense Trotter unitary from an ordered sparse-Pauli term list.

    ``method`` is passed to
    :func:`qhat.common.trotter_flattened.get_trotterization_coefficients`, so
    first-, second-, fourth-, and other QHAT-supported ramped formulas share the
    same coefficient convention as the bloq implementation.
    """
    if num_steps <= 0:
        raise ValueError(f"num_steps must be positive, got {num_steps}.")
    if n_qubits < 0:
        raise ValueError(f"n_qubits must be nonnegative, got {n_qubits}.")

    dim = 2**n_qubits
    unitary = np.eye(dim, dtype=complex)

    if len(term_list) == 0:
        return unitary

    coefficients = get_trotterization_coefficients(method)
    expanded_sequence = expand_ramped_trotterization(
        num_terms=len(term_list),
        coefficients=coefficients,
        num_steps=num_steps,
        combine_terms=combine_terms,
    )
    dt = time / num_steps

    matrix_cache: dict[tuple[tuple[int, str], ...], np.ndarray] = {}

    for term_idx, scale in expanded_sequence:
        key, coeff = term_list[term_idx]
        if key not in matrix_cache:
            matrix_cache[key] = pauli_key_to_matrix(key, n_qubits)
        unitary = unitary @ expm(-1j * dt * scale * coeff * matrix_cache[key])

    return unitary


def trotter_error_matrix(
    hamiltonian: np.ndarray,
    term_list: TermList,
    time: float,
    num_steps: int,
    n_qubits: int,
    method=1,
    combine_terms: bool = True,
) -> np.ndarray:
    """Return ``U_trotter - U_exact`` for the supplied product formula."""
    return dense_trotter_unitary(
        term_list=term_list,
        time=time,
        num_steps=num_steps,
        n_qubits=n_qubits,
        method=method,
        combine_terms=combine_terms,
    ) - exact_unitary(hamiltonian, time)


def trotter_error_norm(
    hamiltonian: np.ndarray,
    term_list: TermList,
    time: float,
    num_steps: int,
    n_qubits: int,
    method=1,
    combine_terms: bool = True,
    ord=2,
) -> float:
    """Return ``||U_trotter - U_exact||``."""
    error = trotter_error_matrix(
        hamiltonian=hamiltonian,
        term_list=term_list,
        time=time,
        num_steps=num_steps,
        n_qubits=n_qubits,
        method=method,
        combine_terms=combine_terms,
    )
    return float(norm(error, ord=ord))


def first_order_trotter_unitary(term_list, t, r, n_qubits):
    """Compatibility wrapper for the uploaded first-order API."""
    return dense_trotter_unitary(term_list, t, r, n_qubits, method=1)


def second_order_trotter_unitary(term_list, t, r, n_qubits):
    """Compatibility wrapper for the uploaded second-order API."""
    return dense_trotter_unitary(term_list, t, r, n_qubits, method=2)


def first_order_trotter_error_matrix(H, term_list, t, r, n_qubits):
    """Compatibility wrapper for the uploaded first-order API."""
    return trotter_error_matrix(H, term_list, t, r, n_qubits, method=1)


def first_order_trotter_error_norm(H, term_list, t, r, n_qubits, ord=2):
    """Compatibility wrapper for the uploaded first-order API."""
    return trotter_error_norm(H, term_list, t, r, n_qubits, method=1, ord=ord)


def second_order_trotter_error_matrix(H, term_list, t, r, n_qubits):
    """Compatibility wrapper for the uploaded second-order API."""
    return trotter_error_matrix(H, term_list, t, r, n_qubits, method=2)


def second_order_trotter_error_norm(H, term_list, t, r, n_qubits, ord=2):
    """Compatibility wrapper for the uploaded second-order API."""
    return trotter_error_norm(H, term_list, t, r, n_qubits, method=2, ord=ord)
