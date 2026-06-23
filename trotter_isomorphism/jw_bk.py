"""JW/BK Pauli-string isomorphism and Clifford construction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from openfermion import count_qubits
from openfermion.transforms import bravyi_kitaev, jordan_wigner

from qhat.trotter_isomorphism.majorana import (
    fermion_operator_to_majorana_operator_dict,
    majorana_monomial_as_fermion_operator,
)
from qhat.trotter_isomorphism.pauli_keys import (
    PauliKey,
    anticommutes,
    clean_qubit_operator,
    format_pauli_key,
    pauli_key,
    qubit_operator_to_matrix,
)
from qhat.trotter_isomorphism.states import bits_from_index, index_from_bits


@dataclass(frozen=True)
class JWBKAnalysis:
    """Reusable data for matched JW/BK comparison experiments."""

    fermion_op: Any
    n_qubits: int
    jw_to_bk: dict[PauliKey, PauliKey]
    bk_to_jw: dict[PauliKey, PauliKey]
    details: dict[str, Any]
    jw_direct: Any
    bk_direct: Any
    jw_order: list[PauliKey]
    jw_terms: list[tuple[PauliKey, complex]]
    bk_terms: list[tuple[PauliKey, complex]]
    H_jw: np.ndarray | None = None
    H_bk: np.ndarray | None = None

    def as_dict(self) -> dict[str, Any]:
        """Return a dictionary matching the original standalone scripts."""
        return {
            "fermion_op": self.fermion_op,
            "n_qubits": self.n_qubits,
            "jw_to_bk": self.jw_to_bk,
            "bk_to_jw": self.bk_to_jw,
            "details": self.details,
            "jw_direct": self.jw_direct,
            "bk_direct": self.bk_direct,
            "jw_order": self.jw_order,
            "jw_terms": self.jw_terms,
            "bk_terms": self.bk_terms,
            "H_jw": self.H_jw,
            "H_bk": self.H_bk,
        }


def count_qubits_in_fermion_operator(fermion_op) -> int:
    """Infer the number of fermionic modes/qubits from an OpenFermion operator."""
    return count_qubits(fermion_op)


def qubit_support_from_majorana_labels(
    majorana_terms: dict[tuple[int, ...], complex],
    transform,
    n_qubits: int,
    atol: float = 1e-12,
) -> dict[tuple[int, ...], tuple[PauliKey, complex]]:
    """Map each Majorana monomial label to one Pauli string under a transform."""
    out = {}

    for label, majorana_coeff in majorana_terms.items():
        fop = majorana_monomial_as_fermion_operator(label)

        if transform is jordan_wigner:
            qop = clean_qubit_operator(transform(fop), atol=atol)
        else:
            qop = clean_qubit_operator(transform(fop, n_qubits=n_qubits), atol=atol)

        if len(qop.terms) != 1:
            raise RuntimeError(
                f"Majorana monomial {label} mapped to {len(qop.terms)} Pauli terms, "
                f"expected exactly 1. Terms: {qop.terms}"
            )

        [(pterm, phase_coeff)] = list(qop.terms.items())
        out[label] = (pauli_key(pterm), majorana_coeff * phase_coeff)

    return out


def find_jw_bk_isomorphism(fermion_op, n_qubits: int | None = None, atol: float = 1e-10):
    """Find the collected-support isomorphism between JW and BK Pauli terms."""
    if n_qubits is None:
        n_qubits = count_qubits_in_fermion_operator(fermion_op)

    majorana_terms = fermion_operator_to_majorana_operator_dict(fermion_op, atol=atol)
    jw_by_label = qubit_support_from_majorana_labels(
        majorana_terms,
        jordan_wigner,
        n_qubits,
        atol=atol,
    )
    bk_by_label = qubit_support_from_majorana_labels(
        majorana_terms,
        bravyi_kitaev,
        n_qubits,
        atol=atol,
    )

    jw_to_bk_candidates: dict[PauliKey, set[PauliKey]] = {}
    bk_to_jw_candidates: dict[PauliKey, set[PauliKey]] = {}
    details_by_label = {}

    for label in majorana_terms:
        jw_p, jw_coeff = jw_by_label[label]
        bk_p, bk_coeff = bk_by_label[label]
        jw_to_bk_candidates.setdefault(jw_p, set()).add(bk_p)
        bk_to_jw_candidates.setdefault(bk_p, set()).add(jw_p)
        details_by_label[label] = {
            "majorana_coeff": majorana_terms[label],
            "jw_pauli": jw_p,
            "jw_coeff_from_label": jw_coeff,
            "bk_pauli": bk_p,
            "bk_coeff_from_label": bk_coeff,
        }

    jw_to_bk = {}
    for jw_p, bk_set in jw_to_bk_candidates.items():
        if len(bk_set) != 1:
            raise RuntimeError(f"JW Pauli {jw_p} maps to multiple BK Paulis: {bk_set}")
        jw_to_bk[jw_p] = next(iter(bk_set))

    bk_to_jw = {}
    for bk_p, jw_set in bk_to_jw_candidates.items():
        if len(jw_set) != 1:
            raise RuntimeError(f"BK Pauli {bk_p} maps to multiple JW Paulis: {jw_set}")
        bk_to_jw[bk_p] = next(iter(jw_set))

    jw_direct = clean_qubit_operator(jordan_wigner(fermion_op), atol=atol)
    bk_direct = clean_qubit_operator(bravyi_kitaev(fermion_op, n_qubits=n_qubits), atol=atol)

    jw_support_direct = {pauli_key(p) for p in jw_direct.terms}
    bk_support_direct = {pauli_key(p) for p in bk_direct.terms}
    mapped_jw_support = {jw_to_bk[p] for p in jw_to_bk}

    if jw_support_direct != set(jw_to_bk):
        missing = jw_support_direct - set(jw_to_bk)
        extra = set(jw_to_bk) - jw_support_direct
        raise RuntimeError(
            "JW support mismatch after Majorana labeling.\n"
            f"Missing from label map: {missing}\n"
            f"Extra in label map: {extra}"
        )

    if bk_support_direct != mapped_jw_support:
        missing = bk_support_direct - mapped_jw_support
        extra = mapped_jw_support - bk_support_direct
        raise RuntimeError(
            "BK support mismatch after applying JW->BK map.\n"
            f"Missing from mapped support: {missing}\n"
            f"Extra in mapped support: {extra}"
        )

    details = {
        "n_qubits": n_qubits,
        "majorana_terms": majorana_terms,
        "details_by_label": details_by_label,
        "jw_direct": jw_direct,
        "bk_direct": bk_direct,
    }
    return jw_to_bk, bk_to_jw, details


def verify_graph_isomorphism(jw_to_bk: dict[PauliKey, PauliKey]) -> bool:
    """Verify that a JW -> BK map preserves pairwise anticommutation."""
    jw_vertices = list(jw_to_bk.keys())

    for i, p in enumerate(jw_vertices):
        for q in jw_vertices[i + 1 :]:
            lhs = anticommutes(p, q)
            rhs = anticommutes(jw_to_bk[p], jw_to_bk[q])

            if lhs != rhs:
                raise RuntimeError(
                    "Commutation mismatch:\n"
                    f"  JW: {p}, {q}, anticommutes={lhs}\n"
                    f"  BK: {jw_to_bk[p]}, {jw_to_bk[q]}, anticommutes={rhs}"
                )

    return True


def build_matched_term_lists(jw_to_bk, jw_direct, bk_direct, jw_order):
    """Build matched JW and BK ordered term lists for Trotter comparisons."""
    jw_terms = [(p, jw_direct.terms[p]) for p in jw_order]
    bk_terms = [(jw_to_bk[p], bk_direct.terms[jw_to_bk[p]]) for p in jw_order]
    return jw_terms, bk_terms


def default_jw_order(jw_to_bk, remove_identity: bool = False) -> list[PauliKey]:
    """Return a stable JW order suitable for inducing the BK order."""
    order = sorted(jw_to_bk.keys(), key=format_pauli_key)
    if remove_identity:
        order = [p for p in order if p != ()]
    return order


def build_jw_bk_analysis(
    fermion_op,
    n_qubits: int | None = None,
    atol: float = 1e-10,
    remove_identity: bool = False,
    build_matrices: bool = True,
    jw_order: Sequence[PauliKey] | None = None,
) -> JWBKAnalysis:
    """Construct matched JW/BK data from a second-quantized Hamiltonian."""
    if n_qubits is None:
        n_qubits = count_qubits_in_fermion_operator(fermion_op)

    jw_to_bk, bk_to_jw, details = find_jw_bk_isomorphism(
        fermion_op,
        n_qubits=n_qubits,
        atol=atol,
    )
    verify_graph_isomorphism(jw_to_bk)

    jw_direct = details["jw_direct"]
    bk_direct = details["bk_direct"]

    if jw_order is None:
        jw_order = default_jw_order(jw_to_bk, remove_identity=remove_identity)
    else:
        jw_order = list(jw_order)
        if remove_identity:
            jw_order = [p for p in jw_order if p != ()]

    missing = [p for p in jw_order if p not in jw_to_bk]
    if missing:
        raise ValueError(f"jw_order contains Pauli keys not present in jw_to_bk: {missing}")

    jw_terms, bk_terms = build_matched_term_lists(jw_to_bk, jw_direct, bk_direct, jw_order)

    H_jw = qubit_operator_to_matrix(jw_direct, n_qubits) if build_matrices else None
    H_bk = qubit_operator_to_matrix(bk_direct, n_qubits) if build_matrices else None

    return JWBKAnalysis(
        fermion_op=fermion_op,
        n_qubits=n_qubits,
        jw_to_bk=jw_to_bk,
        bk_to_jw=bk_to_jw,
        details=details,
        jw_direct=jw_direct,
        bk_direct=bk_direct,
        jw_order=list(jw_order),
        jw_terms=jw_terms,
        bk_terms=bk_terms,
        H_jw=H_jw,
        H_bk=H_bk,
    )


def gf2_inverse(A) -> np.ndarray:
    """Invert a binary matrix over GF(2)."""
    A = np.array(A, dtype=np.uint8) % 2
    n = A.shape[0]

    if A.shape != (n, n):
        raise ValueError("A must be a square matrix.")

    aug = np.concatenate([A.copy(), np.eye(n, dtype=np.uint8)], axis=1)

    for col in range(n):
        pivot = None

        for row in range(col, n):
            if aug[row, col] == 1:
                pivot = row
                break

        if pivot is None:
            raise ValueError("Matrix is not invertible over GF(2).")

        if pivot != col:
            aug[[col, pivot]] = aug[[pivot, col]]

        for row in range(n):
            if row != col and aug[row, col] == 1:
                aug[row, :] ^= aug[col, :]

    return aug[:, n:] % 2


def z_image_matrix_from_jw_to_bk(jw_to_bk, n_qubits: int) -> np.ndarray:
    """Extract the binary Z-image matrix from a JW -> BK Pauli map."""
    A = np.zeros((n_qubits, n_qubits), dtype=np.uint8)

    for j in range(n_qubits):
        jw_zj = ((j, "Z"),)

        if jw_zj not in jw_to_bk:
            raise ValueError(
                f"Single-qubit Z{j} not found in jw_to_bk map. "
                "This method needs Z_j terms in the JW support."
            )

        bk_image = jw_to_bk[jw_zj]

        for q, op in bk_image:
            if op != "Z":
                raise ValueError(f"Image of Z{j} is not Z-only: {format_pauli_key(bk_image)}")
            A[q, j] = 1

    return A


def build_clifford_from_z_map(jw_to_bk, n_qubits: int):
    """Build the JW -> BK computational-basis permutation unitary."""
    A = z_image_matrix_from_jw_to_bk(jw_to_bk=jw_to_bk, n_qubits=n_qubits)
    B = gf2_inverse(A).T % 2
    dim = 2**n_qubits
    C = np.zeros((dim, dim), dtype=complex)

    for jw_index in range(dim):
        n_bits = bits_from_index(jw_index, n_qubits)
        b_bits = (B @ n_bits) % 2
        bk_index = index_from_bits(b_bits)
        C[bk_index, jw_index] = 1.0

    return C, B, A
