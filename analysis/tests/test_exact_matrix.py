"""
Integration tests for exact Hamiltonian matrix computation.

Tests the full pipeline from Hamiltonian to exact matrix.
"""

import numpy as np
import pytest
import tempfile
import os

from qhat.analysis.hamiltonian import Hamiltonian, LinearCombinationOfPauliStrings
from qhat.analysis.matrix_operations import PauliStringOperator


# =================================================================================================
# Test Hamiltonian.to_matrix() Method
# =================================================================================================

def test_hamiltonian_to_matrix_small_system():
    """Test converting a small Hamiltonian to dense matrix."""
    # Create a simple 2-qubit Hamiltonian
    pauli_data = {
        'ZZ': 0.5,
        'XX': -0.5,
        'II': 1.0
    }
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    # Convert to matrix (should be dense for 2 qubits)
    H_matrix = hamiltonian.to_matrix()

    # Verify it's a dense numpy array
    assert isinstance(H_matrix, np.ndarray)
    assert H_matrix.shape == (4, 4)

    # Verify it's Hermitian
    assert np.allclose(H_matrix, H_matrix.conj().T)


def test_hamiltonian_to_matrix_large_system():
    """Test converting a large Hamiltonian to matrix-free operator."""
    # Create a 20-qubit Hamiltonian (dimension 2^20 = 1,048,576)
    pauli_data = {'X' * 20: 1.0}
    lcps = LinearCombinationOfPauliStrings(num_qubits=20, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    # Convert to matrix (should be matrix-free for 20 qubits)
    H_op = hamiltonian.to_matrix()

    # Verify it's a matrix-free operator
    assert isinstance(H_op, PauliStringOperator)
    assert H_op.shape == (2**20, 2**20)


def test_hamiltonian_to_matrix_force_dense():
    """Test forcing dense representation."""
    pauli_data = {'ZZ': 0.5}
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    H_matrix = hamiltonian.to_matrix(force_dense=True)
    assert isinstance(H_matrix, np.ndarray)


def test_hamiltonian_to_matrix_force_sparse():
    """Test forcing sparse representation."""
    pauli_data = {'ZZ': 0.5}
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    H_op = hamiltonian.to_matrix(force_sparse=True)
    assert isinstance(H_op, PauliStringOperator)


def test_hamiltonian_to_matrix_known_eigenvalues():
    """Test that exact matrix has correct eigenvalues for known Hamiltonian."""
    # Simple 2-qubit Hamiltonian: H = Z⊗Z
    # Eigenvalues: [1, -1, -1, 1] (for basis |00⟩, |01⟩, |10⟩, |11⟩)
    pauli_data = {'ZZ': 1.0}
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    H_matrix = hamiltonian.to_matrix()

    eigenvalues = np.linalg.eigvalsh(H_matrix)
    expected = np.array([-1, -1, 1, 1])  # Sorted

    assert np.allclose(sorted(eigenvalues), sorted(expected))


def test_hamiltonian_to_matrix_identity():
    """Test Hamiltonian that's just identity."""
    pauli_data = {'II': 2.5}
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    H_matrix = hamiltonian.to_matrix()

    expected = 2.5 * np.eye(4, dtype=complex)
    assert np.allclose(H_matrix, expected)


def test_hamiltonian_to_matrix_pauli_x():
    """Test Hamiltonian with Pauli X."""
    # H = X⊗I on 2 qubits
    pauli_data = {'XI': 1.0}
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    H_matrix = hamiltonian.to_matrix()

    # Verify it's Hermitian
    assert np.allclose(H_matrix, H_matrix.conj().T)

    # X⊗I eigenvalues should be ±1 (twice each, since I tensor'd with it)
    eigenvalues = np.linalg.eigvalsh(H_matrix)
    assert np.allclose(sorted(eigenvalues), [-1, -1, 1, 1])


def test_hamiltonian_to_matrix_multiple_terms():
    """Test Hamiltonian with multiple Pauli terms."""
    pauli_data = {
        'XX': 1.0,
        'YY': 1.0,
        'ZZ': -1.0,
        'II': 0.5
    }
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    H_matrix = hamiltonian.to_matrix()

    # Verify Hermiticity
    assert np.allclose(H_matrix, H_matrix.conj().T)

    # Verify eigenvalues are real
    eigenvalues = np.linalg.eigvalsh(H_matrix)
    assert all(np.isreal(eigenvalues))


# =================================================================================================
# Test Matrix Properties
# =================================================================================================

def test_exact_matrix_hermiticity():
    """Test that exact matrix is always Hermitian."""
    test_cases = [
        {'ZZ': 1.0},
        {'XX': 0.5, 'YY': 0.3, 'ZZ': -0.2},
        {'XII': 1.0, 'IXI': 1.0, 'IIX': 1.0},
    ]

    for pauli_data in test_cases:
        num_qubits = len(list(pauli_data.keys())[0])
        lcps = LinearCombinationOfPauliStrings(num_qubits=num_qubits, dense=pauli_data)
        hamiltonian = Hamiltonian(lcps)

        # Force dense to make Hermiticity check easy
        H_matrix = hamiltonian.to_matrix(force_dense=True)

        assert np.allclose(H_matrix, H_matrix.conj().T), \
            f"Matrix not Hermitian for {pauli_data}"


def test_exact_matrix_real_eigenvalues():
    """Test that exact matrix has real eigenvalues (Hermitian property)."""
    pauli_data = {'XX': 1.0, 'YY': 0.5, 'ZZ': -0.3}
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    H_matrix = hamiltonian.to_matrix()
    eigenvalues = np.linalg.eigvalsh(H_matrix)

    # All eigenvalues should be real (imaginary part ~0)
    assert np.allclose(eigenvalues.imag, 0)


def test_exact_matrix_energy_bounds():
    """Test that eigenvalues are within expected energy bounds."""
    pauli_data = {'ZI': 1.0, 'IZ': 1.0, 'ZZ': -2.0}
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    H_matrix = hamiltonian.to_matrix()
    eigenvalues = np.linalg.eigvalsh(H_matrix)

    # For this Hamiltonian, eigenvalues should be in [-4, 4]
    # (rough bound: sum of absolute values)
    assert all(eigenvalues >= -4.1)  # Small tolerance
    assert all(eigenvalues <= 4.1)


# =================================================================================================
# Test Matrix-Free Operator Behavior
# =================================================================================================

def test_matrix_free_operator_matvec():
    """Test matrix-free operator gives correct results."""
    pauli_data = {'XX': 1.0, 'ZZ': -1.0}
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    # Get both dense and sparse versions
    H_dense = hamiltonian.to_matrix(force_dense=True)
    H_sparse = hamiltonian.to_matrix(force_sparse=True)

    # Test on a few different states
    test_states = [
        np.array([1, 0, 0, 0], dtype=complex),
        np.array([0, 1, 0, 0], dtype=complex),
        np.array([1, 1, 1, 1], dtype=complex) / 2,
    ]

    for state in test_states:
        result_dense = H_dense @ state
        result_sparse = H_sparse @ state
        assert np.allclose(result_dense, result_sparse), \
            "Dense and sparse results don't match"


# =================================================================================================
# Integration with Existing Code
# =================================================================================================

def test_hamiltonian_from_pauli_strings():
    """Test creating Hamiltonian from Pauli strings and converting to matrix."""
    # This mimics the typical workflow
    pauli_strings_sparse = {
        ((0, 'Z'), (1, 'Z')): 0.5,  # ZZ with sparse format
        ((0, 'X'), (1, 'X')): -0.5,  # XX
        (): 1.0  # Identity
    }
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, sparse=pauli_strings_sparse)
    hamiltonian = Hamiltonian(lcps)

    H_matrix = hamiltonian.to_matrix()

    assert H_matrix.shape == (4, 4)
    assert np.allclose(H_matrix, H_matrix.conj().T)


def test_three_qubit_system():
    """Test with a 3-qubit system."""
    pauli_data = {
        'ZZZ': 1.0,
        'XXX': 0.5,
        'III': -0.5
    }
    lcps = LinearCombinationOfPauliStrings(num_qubits=3, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    H_matrix = hamiltonian.to_matrix()

    assert H_matrix.shape == (8, 8)
    assert np.allclose(H_matrix, H_matrix.conj().T)

    # Check that identity term is present
    eigenvalues = np.linalg.eigvalsh(H_matrix)
    # With -0.5 identity, all eigenvalues shifted by -0.5
    assert np.min(eigenvalues) <= -0.5


# =================================================================================================
# Edge Cases
# =================================================================================================

def test_single_qubit_hamiltonian():
    """Test with a single qubit."""
    pauli_data = {'X': 1.0}
    lcps = LinearCombinationOfPauliStrings(num_qubits=1, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    H_matrix = hamiltonian.to_matrix()

    assert H_matrix.shape == (2, 2)
    # Pauli X eigenvalues are ±1
    eigenvalues = np.linalg.eigvalsh(H_matrix)
    assert np.allclose(sorted(eigenvalues), [-1, 1])


def test_empty_hamiltonian():
    """Test with only identity term."""
    pauli_data = {'II': 1.0}
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    hamiltonian = Hamiltonian(lcps)

    H_matrix = hamiltonian.to_matrix()

    expected = np.eye(4, dtype=complex)
    assert np.allclose(H_matrix, expected)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
