"""Correctness tests for the Pauli representation layer.

These verify the operator container against hand-computed dense matrices in the
stated ``(x, z)``-bitmask convention, the weight bookkeeping (``lambda``, sorting),
and the cumulative-budget truncation -- including its off-by-one
boundary.  The openfermion conversion path is exercised only when openfermion is
installed; everything else is ``numpy``-only.
"""

import numpy as np
import pytest

from qhat.prqpe.representation import (
    build_representation,
    sorted_weights,
    truncate_small_terms,
)
from qhat.prqpe.representation.pauli import PauliHamiltonian, popcount

# Single-qubit Pauli matrices, qubit 0 = the only (most significant) bit.
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


# ----------------------------------------------------------------------------
# Single-qubit dense matrices against the (x, z) convention.
# ----------------------------------------------------------------------------

def test_single_qubit_x():
    # x = 1, z = 0 -> P = i^0 X^1 Z^0 = X.
    ham = PauliHamiltonian(1, [1], [0], [1.0])
    np.testing.assert_allclose(ham.to_dense(), X)


def test_single_qubit_z():
    # x = 0, z = 1 -> P = i^0 X^0 Z^1 = Z.
    ham = PauliHamiltonian(1, [0], [1], [1.0])
    np.testing.assert_allclose(ham.to_dense(), Z)


def test_single_qubit_y():
    # x = z = 1 -> P = i^popcount(1) X Z = i X Z = Y.
    ham = PauliHamiltonian(1, [1], [1], [1.0])
    np.testing.assert_allclose(ham.to_dense(), Y)


def test_single_qubit_identity():
    # x = z = 0 -> P = I.
    ham = PauliHamiltonian(1, [0], [0], [1.0])
    np.testing.assert_allclose(ham.to_dense(), I2)


def test_coefficient_scales_matrix():
    ham = PauliHamiltonian(1, [0], [1], [-2.5])
    np.testing.assert_allclose(ham.to_dense(), -2.5 * Z)


# ----------------------------------------------------------------------------
# Two-qubit term, qubit 0 = MSB (left tensor factor).
# ----------------------------------------------------------------------------

def test_two_qubit_xz():
    # X on qubit 0 (bit 1) tensor Z on qubit 1 (bit 0): x = 0b10 = 2, z = 0b01 = 1.
    # popcount(2 & 1) = 0, so the head phase is i^0 = 1 and P = X (x) Z (x) Z.
    ham = PauliHamiltonian(2, [2], [1], [1.0])
    np.testing.assert_allclose(ham.to_dense(), np.kron(X, Z))


def test_two_qubit_yy():
    # Y on both qubits: x = 0b11 = 3, z = 0b11 = 3, head phase i^popcount(3) = i^2 = -1.
    ham = PauliHamiltonian(2, [3], [3], [1.0])
    np.testing.assert_allclose(ham.to_dense(), np.kron(Y, Y))


def test_two_qubit_sum_is_hermitian():
    ham = PauliHamiltonian(2, [2, 3], [1, 3], [0.5, -0.75])
    dense = ham.to_dense()
    np.testing.assert_allclose(dense, 0.5 * np.kron(X, Z) - 0.75 * np.kron(Y, Y))
    np.testing.assert_allclose(dense, dense.conj().T)


def test_to_sparse_matches_to_dense():
    ham = PauliHamiltonian(2, [2, 3, 0], [1, 3, 2], [0.5, -0.75, 1.25])
    np.testing.assert_allclose(ham.to_sparse().toarray(), ham.to_dense())


# ----------------------------------------------------------------------------
# popcount and the phase head.
# ----------------------------------------------------------------------------

def test_popcount_matches_python():
    vals = np.array([0, 1, 2, 3, 7, 255, 1023], dtype=np.int64)
    expected = np.array([bin(int(v)).count("1") for v in vals])
    np.testing.assert_array_equal(popcount(vals), expected)


# ----------------------------------------------------------------------------
# pauli_apply_matrix and rotation_apply_matrix against explicit matrices.
# ----------------------------------------------------------------------------

def test_pauli_apply_matrix_left_multiplies():
    rng = np.random.default_rng(1)
    ham = PauliHamiltonian(2, [2, 3], [1, 3], [1.0, 1.0])
    M = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))
    for l in range(ham.L):
        P = ham.term_matrix(l).toarray()
        np.testing.assert_allclose(ham.pauli_apply_matrix(l, M), P @ M)


def test_rotation_apply_matrix_matches_expm():
    from scipy.linalg import expm

    rng = np.random.default_rng(2)
    ham = PauliHamiltonian(2, [2, 3], [1, 3], [1.0, -1.0])
    M = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))
    phi = 0.37
    for l in range(ham.L):
        P = ham.term_matrix(l).toarray()
        expected = expm(-1j * phi * P) @ M
        np.testing.assert_allclose(ham.rotation_apply_matrix(l, phi, M), expected,
                                   atol=1e-12)


# ----------------------------------------------------------------------------
# Weight bookkeeping: lambda, sorting.
# ----------------------------------------------------------------------------

def test_lambda_is_sum_abs_coeffs():
    coeffs = [0.5, -2.0, 1.5, -0.25]
    ham = PauliHamiltonian(2, [0, 1, 2, 3], [0, 1, 2, 3], coeffs)
    assert ham.lambda_weight == pytest.approx(sum(abs(c) for c in coeffs))


def test_sorted_by_weight_is_descending():
    coeffs = np.array([0.5, -2.0, 1.5, -0.25, 3.0])
    ham = PauliHamiltonian(3, [0, 1, 2, 3, 4], [0, 1, 2, 3, 4], coeffs)
    ordered = ham.sorted_by_weight()
    w = np.abs(ordered.coeffs)
    assert np.all(np.diff(w) <= 0)
    np.testing.assert_array_equal(w, np.sort(np.abs(coeffs))[::-1])


def test_sorted_weights_helper_contract():
    coeffs = np.array([0.5, -2.0, 1.5, -0.25, 3.0])
    ham = PauliHamiltonian(3, [0, 1, 2, 3, 4], [0, 1, 2, 3, 4], coeffs)
    np.testing.assert_array_equal(
        sorted_weights(ham), np.abs(ham.sorted_by_weight().coeffs))
    assert np.all(np.diff(sorted_weights(ham)) <= 0)


def test_subset_preserves_order():
    ham = PauliHamiltonian(2, [0, 1, 2, 3], [3, 2, 1, 0], [1.0, 2.0, 3.0, 4.0])
    sub = ham.subset([3, 0])
    np.testing.assert_array_equal(sub.coeffs, [4.0, 1.0])
    np.testing.assert_array_equal(sub.x_masks, [3, 0])
    np.testing.assert_array_equal(sub.z_masks, [0, 3])


# ----------------------------------------------------------------------------
# truncate_small_terms: cumulative-budget rule and its boundary.
# ----------------------------------------------------------------------------

def test_truncate_drops_smallest_within_budget():
    # Tail weights 1e-4 + 2e-4 + 3e-4 = 6e-4 <= 1e-3; adding the next (5e-4)
    # would exceed the budget, so exactly three terms are dropped.
    coeffs = np.array([1.0, 0.5, 5e-4, 3e-4, 2e-4, 1e-4])
    ham = PauliHamiltonian(3, np.arange(6), np.arange(6), coeffs)
    out = truncate_small_terms(ham, weight_threshold=1e-3)
    np.testing.assert_array_equal(np.abs(out.coeffs), [1.0, 0.5, 5e-4])
    assert out.L == 3


def test_truncate_boundary_exactly_at_budget_drops_term():
    # Cumulative dropped weight hits the budget exactly: side="right" keeps it
    # within budget, so the term at the boundary is dropped.
    coeffs = np.array([1.0, 6e-4, 4e-4])  # tail 4e-4, then 4e-4 + 6e-4 = 1e-3
    ham = PauliHamiltonian(2, [0, 1, 2], [0, 1, 2], coeffs)
    out = truncate_small_terms(ham, weight_threshold=1e-3)
    # 4e-4 dropped (cum 4e-4 <= 1e-3); 6e-4 dropped (cum 1e-3 <= 1e-3); 1.0 kept.
    np.testing.assert_array_equal(np.abs(out.coeffs), [1.0])
    assert out.L == 1


def test_truncate_off_by_one_just_over_budget_keeps_term():
    # The same shape, but the second-smallest is nudged so the cumulative drop
    # would be 1.0001e-3 > 1e-3: that term must be kept (off-by-one boundary).
    coeffs = np.array([1.0, 6.01e-4, 4e-4])  # 4e-4 + 6.01e-4 = 1.001e-3 > 1e-3
    ham = PauliHamiltonian(2, [0, 1, 2], [0, 1, 2], coeffs)
    out = truncate_small_terms(ham, weight_threshold=1e-3)
    np.testing.assert_array_equal(np.abs(out.coeffs), [1.0, 6.01e-4])
    assert out.L == 2


def test_truncate_drops_nothing_when_budget_below_smallest():
    coeffs = np.array([1.0, 0.5, 0.01])
    ham = PauliHamiltonian(2, [0, 1, 2], [0, 1, 2], coeffs)
    out = truncate_small_terms(ham, weight_threshold=1e-3)
    assert out.L == ham.L
    np.testing.assert_array_equal(np.abs(out.coeffs), [1.0, 0.5, 0.01])


def test_truncate_total_dropped_weight_within_budget():
    rng = np.random.default_rng(7)
    coeffs = rng.standard_normal(50) * 1e-3
    ham = PauliHamiltonian(6, np.arange(50), np.arange(50), coeffs)
    budget = 5e-3
    out = truncate_small_terms(ham, weight_threshold=budget)
    dropped = ham.lambda_weight - out.lambda_weight
    assert dropped <= budget + 1e-15
    assert np.all(np.diff(np.abs(out.coeffs)) <= 0)


def test_truncate_output_is_sorted_descending():
    coeffs = np.array([0.5, 1.0, 2e-4, 0.3])
    ham = PauliHamiltonian(2, [0, 1, 2, 3], [0, 1, 2, 3], coeffs)
    out = truncate_small_terms(ham, weight_threshold=1e-3)
    assert np.all(np.diff(np.abs(out.coeffs)) <= 0)


# ----------------------------------------------------------------------------
# build_representation: pass-through and the unimplemented arm.
# ----------------------------------------------------------------------------

def test_build_representation_passthrough():
    ham = PauliHamiltonian(2, [2, 3], [1, 3], [0.5, -0.75])
    assert build_representation(ham) is ham
    assert build_representation(ham, representation="pauli") is ham


def test_build_representation_unknown_raises():
    ham = PauliHamiltonian(1, [1], [0], [1.0])
    with pytest.raises(ValueError):
        build_representation(ham, representation="qubitization")


# ----------------------------------------------------------------------------
# openfermion conversion path (skipped if openfermion is absent).
# ----------------------------------------------------------------------------

def test_build_representation_from_qubit_operator():
    openfermion = pytest.importorskip("openfermion")
    from openfermion import QubitOperator, get_sparse_operator

    qop = (0.5 * QubitOperator("X0 Z1")
           - 0.75 * QubitOperator("Y0 Y1")
           + 1.3 * QubitOperator(""))  # identity term -> dropped
    ham = build_representation(qop)
    assert ham.n_qubits == 2
    assert ham.L == 2  # identity dropped

    # Matches openfermion's own sparse operator on the non-identity part.
    nonconst = 0.5 * QubitOperator("X0 Z1") - 0.75 * QubitOperator("Y0 Y1")
    ref = get_sparse_operator(nonconst, n_qubits=2).toarray()
    np.testing.assert_allclose(ham.to_dense(), ref, atol=1e-12)


def test_qubit_operator_lambda_matches_coeffs():
    pytest.importorskip("openfermion")
    from openfermion import QubitOperator

    qop = 0.3 * QubitOperator("X0") - 1.2 * QubitOperator("Z0 Z1")
    ham = build_representation(qop)
    assert ham.lambda_weight == pytest.approx(0.3 + 1.2)
