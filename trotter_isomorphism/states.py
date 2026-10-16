"""Computational-basis state helpers."""

from __future__ import annotations

from typing import Sequence

import numpy as np


def bits_from_index(index: int, n_qubits: int) -> np.ndarray:
    """Convert a basis index into bits ``[q0, q1, ...]``.

    Qubit 0 is the leftmost Kronecker factor, matching
    ``qhat.common.pauli_utils.pauli_string_to_matrix``.
    """
    if index < 0:
        raise ValueError(f"index must be nonnegative, got {index}.")
    if n_qubits < 0:
        raise ValueError(f"n_qubits must be nonnegative, got {n_qubits}.")
    if index >= 2**n_qubits:
        raise ValueError(f"index {index} does not fit in {n_qubits} qubits.")

    return np.array(
        [(index >> (n_qubits - 1 - q)) & 1 for q in range(n_qubits)],
        dtype=np.uint8,
    )


def index_from_bits(bits: Sequence[int]) -> int:
    """Convert bits ``[q0, q1, ...]`` into a computational-basis index."""
    idx = 0
    n_qubits = len(bits)

    for q, bit in enumerate(bits):
        if bit not in (0, 1):
            raise ValueError(f"Bit values must be 0 or 1, got {bit!r}.")
        idx |= int(bit) << (n_qubits - 1 - q)

    return idx


def bitstring_state(bits: Sequence[int]) -> np.ndarray:
    """Create the computational-basis state ``|bits>``."""
    bits = tuple(int(bit) for bit in bits)
    dim = 2 ** len(bits)
    state = np.zeros(dim, dtype=complex)
    state[index_from_bits(bits)] = 1.0
    return state
