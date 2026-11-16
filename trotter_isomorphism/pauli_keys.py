"""Utilities for sparse OpenFermion-style Pauli keys.

The dense matrix conversion intentionally reuses QHAT's existing dense Pauli
string helper instead of carrying a second copy of the Pauli matrices here.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np

from qhat.common.pauli_utils import pauli_string_to_matrix, validate_pauli_string

PauliKey = tuple[tuple[int, str], ...]


def clean_qubit_operator(op, atol: float = 1e-12):
    """Return an OpenFermion ``QubitOperator`` with tiny coefficients removed."""
    from openfermion import QubitOperator

    out = QubitOperator()

    for term, coeff in op.terms.items():
        if abs(coeff) > atol:
            out += QubitOperator(term, coeff)

    out.compress(abs_tol=atol)
    return out


def pauli_key(term: Iterable[tuple[int, str]]) -> PauliKey:
    """Return a sorted, hashable Pauli-key representation."""
    return tuple(sorted((int(q), str(op)) for q, op in term))


def pauli_key_to_dense_string(key: Sequence[tuple[int, str]], n_qubits: int) -> str:
    """Convert ``((0, "X"), (2, "Z"))`` into a dense string like ``"XIZ"``."""
    if n_qubits < 0:
        raise ValueError(f"n_qubits must be nonnegative, got {n_qubits}.")

    dense = ["I"] * n_qubits

    for q, op in key:
        if q < 0 or q >= n_qubits:
            raise ValueError(f"Qubit index {q} is outside 0..{n_qubits - 1}.")
        if op not in {"X", "Y", "Z"}:
            raise ValueError(f"Invalid sparse Pauli operator {op!r}; expected X, Y, or Z.")
        if dense[q] != "I":
            raise ValueError(f"Duplicate Pauli operator supplied for qubit {q}.")
        dense[q] = op

    return "".join(dense)


def dense_string_to_pauli_key(pauli_string: str) -> PauliKey:
    """Convert a dense Pauli string into an OpenFermion-style sparse key."""
    validate_pauli_string(pauli_string)
    return tuple((idx, op) for idx, op in enumerate(pauli_string) if op != "I")


def pauli_key_to_matrix(key: Sequence[tuple[int, str]], n_qubits: int) -> np.ndarray:
    """Convert an OpenFermion-style Pauli key into a dense matrix."""
    return pauli_string_to_matrix(pauli_key_to_dense_string(key, n_qubits))


def qubit_operator_to_matrix(qop, n_qubits: int) -> np.ndarray:
    """Convert an OpenFermion ``QubitOperator`` into a dense matrix."""
    dim = 2**n_qubits
    matrix = np.zeros((dim, dim), dtype=complex)

    for key, coeff in qop.terms.items():
        matrix += coeff * pauli_key_to_matrix(key, n_qubits)

    return matrix


def format_pauli_key(key: Sequence[tuple[int, str]]) -> str:
    """Format a sparse Pauli key as ``X0 Z2``; identity is ``I``."""
    if tuple(key) == ():
        return "I"

    return " ".join(f"{op}{q}" for q, op in key)


def anticommutes(p: Sequence[tuple[int, str]], q: Sequence[tuple[int, str]]) -> bool:
    """Return ``True`` iff two sparse Pauli keys anticommute."""
    pd = dict(p)
    qd = dict(q)
    count = 0

    for idx in set(pd) & set(qd):
        if pd[idx] != qd[idx]:
            count += 1

    return bool(count % 2)


def pauli_weight(key: Sequence[tuple[int, str]]) -> int:
    """Return the number of non-identity factors in a sparse Pauli key."""
    return len(key)
