"""Operator-norm and state-overlap Trotter comparison helpers."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .product_formula import (
    TermList,
    dense_trotter_unitary,
    exact_unitary,
    trotter_error_norm,
)
from .states import bitstring_state


def parse_bitstring(bitstring: str) -> list[int]:
    """Parse a string like ``"1100"`` into ``[1, 1, 0, 0]``."""
    bitstring = bitstring.strip()

    if not all(ch in {"0", "1"} for ch in bitstring):
        raise ValueError("Bitstring must contain only 0 and 1.")

    return [int(ch) for ch in bitstring]


def normalize_state(state: np.ndarray, atol: float = 1e-14) -> np.ndarray:
    """Return a normalized copy of ``state``."""
    state = np.asarray(state, dtype=complex)
    state_norm = np.linalg.norm(state)

    if state_norm < atol:
        raise ValueError("Cannot normalize a zero or near-zero state.")

    return state / state_norm


def trotter_overlap_amplitude(
    hamiltonian: np.ndarray,
    term_list: TermList,
    psi0: np.ndarray,
    time: float,
    num_steps: int,
    n_qubits: int,
    method=1,
    combine_terms: bool = True,
) -> complex:
    """Compute ``<psi0| U_trotter^dagger U_exact |psi0>``."""
    psi0 = normalize_state(psi0)
    exact = exact_unitary(hamiltonian, time)
    trotter = dense_trotter_unitary(
        term_list=term_list,
        time=time,
        num_steps=num_steps,
        n_qubits=n_qubits,
        method=method,
        combine_terms=combine_terms,
    )
    return complex(psi0.conj().T @ (trotter.conj().T @ exact @ psi0))


def fidelity_from_overlap_amplitude(amplitude: complex) -> float:
    """Return ``|amplitude|^2``."""
    return float(abs(amplitude) ** 2)


def infidelity_from_overlap_amplitude(amplitude: complex) -> float:
    """Return ``1 - |amplitude|^2``."""
    return float(1.0 - fidelity_from_overlap_amplitude(amplitude))


def compare_state_dependent_trotter_error(
    psi0_jw: np.ndarray,
    C_jw_to_bk: np.ndarray,
    H_jw: np.ndarray,
    H_bk: np.ndarray,
    jw_terms: TermList,
    bk_terms: TermList,
    time: float,
    num_steps: int,
    n_qubits: int,
    method=1,
    combine_terms: bool = True,
) -> dict[str, object]:
    """Compare JW and BK state-dependent Trotter overlap diagnostics."""
    psi0_jw = normalize_state(psi0_jw)
    psi0_bk = normalize_state(C_jw_to_bk @ psi0_jw)

    amp_jw = trotter_overlap_amplitude(
        hamiltonian=H_jw,
        term_list=jw_terms,
        psi0=psi0_jw,
        time=time,
        num_steps=num_steps,
        n_qubits=n_qubits,
        method=method,
        combine_terms=combine_terms,
    )
    amp_bk = trotter_overlap_amplitude(
        hamiltonian=H_bk,
        term_list=bk_terms,
        psi0=psi0_bk,
        time=time,
        num_steps=num_steps,
        n_qubits=n_qubits,
        method=method,
        combine_terms=combine_terms,
    )

    fidelity_jw = fidelity_from_overlap_amplitude(amp_jw)
    fidelity_bk = fidelity_from_overlap_amplitude(amp_bk)
    infidelity_jw = 1.0 - fidelity_jw
    infidelity_bk = 1.0 - fidelity_bk

    return {
        "psi0_jw": psi0_jw,
        "psi0_bk": psi0_bk,
        "amp_jw": amp_jw,
        "amp_bk": amp_bk,
        "fidelity_jw": fidelity_jw,
        "fidelity_bk": fidelity_bk,
        "infidelity_jw": infidelity_jw,
        "infidelity_bk": infidelity_bk,
        "amplitude_difference": abs(amp_jw - amp_bk),
        "fidelity_difference": abs(fidelity_jw - fidelity_bk),
        "infidelity_difference": abs(infidelity_jw - infidelity_bk),
    }


def compare_jw_bk_operator_norm_errors(
    H_jw: np.ndarray,
    H_bk: np.ndarray,
    jw_terms: TermList,
    bk_terms: TermList,
    time: float,
    num_steps: int,
    n_qubits: int,
    method=1,
    combine_terms: bool = True,
    ord=2,
) -> dict[str, float]:
    """Compare JW and BK dense Trotter error norms."""
    err_jw = trotter_error_norm(
        hamiltonian=H_jw,
        term_list=jw_terms,
        time=time,
        num_steps=num_steps,
        n_qubits=n_qubits,
        method=method,
        combine_terms=combine_terms,
        ord=ord,
    )
    err_bk = trotter_error_norm(
        hamiltonian=H_bk,
        term_list=bk_terms,
        time=time,
        num_steps=num_steps,
        n_qubits=n_qubits,
        method=method,
        combine_terms=combine_terms,
        ord=ord,
    )

    return {
        "error_jw": err_jw,
        "error_bk": err_bk,
        "error_difference": abs(err_jw - err_bk),
    }


def bitstring_overlap_comparison(
    psi0: str | Sequence[int],
    C_jw_to_bk: np.ndarray,
    H_jw: np.ndarray,
    H_bk: np.ndarray,
    jw_terms: TermList,
    bk_terms: TermList,
    time: float,
    num_steps: int,
    n_qubits: int,
    method=1,
    combine_terms: bool = True,
) -> Mapping[str, object]:
    """Run a state-overlap comparison from a JW computational-basis bitstring."""
    bits = parse_bitstring(psi0) if isinstance(psi0, str) else list(psi0)
    if len(bits) != n_qubits:
        raise ValueError(f"Expected {n_qubits} bits, got {len(bits)}.")

    return compare_state_dependent_trotter_error(
        psi0_jw=bitstring_state(bits),
        C_jw_to_bk=C_jw_to_bk,
        H_jw=H_jw,
        H_bk=H_bk,
        jw_terms=jw_terms,
        bk_terms=bk_terms,
        time=time,
        num_steps=num_steps,
        n_qubits=n_qubits,
        method=method,
        combine_terms=combine_terms,
    )
