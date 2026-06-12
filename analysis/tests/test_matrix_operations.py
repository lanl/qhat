"""
Unit tests for matrix operations module.

Tests Pauli string to matrix conversion for both dense and sparse representations.
"""

import numpy as np
import pytest

from qhat.analysis.matrix_operations import (
    pauli_char_to_matrix,
    pauli_string_to_matrix,
    pauli_dict_to_matrix,
    PauliStringOperator,
    get_operator_type,
    create_hamiltonian_operator,
    PAULI_I, PAULI_X, PAULI_Y, PAULI_Z
)


# =================================================================================================
# Test Pauli Character to Matrix
# =================================================================================================

def test_pauli_char_to_matrix_identity():
    """Test conversion of 'I' to identity matrix."""
    result = pauli_char_to_matrix('I')
    expected = PAULI_I
    assert np.allclose(result, expected)


def test_pauli_char_to_matrix_x():
    """Test conversion of 'X' to Pauli-X matrix."""
    result = pauli_char_to_matrix('X')
    expected = PAULI_X
    assert np.allclose(result, expected)


def test_pauli_char_to_matrix_y():
    """Test conversion of 'Y' to Pauli-Y matrix."""
    result = pauli_char_to_matrix('Y')
    expected = PAULI_Y
    assert np.allclose(result, expected)


def test_pauli_char_to_matrix_z():
    """Test conversion of 'Z' to Pauli-Z matrix."""
    result = pauli_char_to_matrix('Z')
    expected = PAULI_Z
    assert np.allclose(result, expected)


def test_pauli_char_to_matrix_case_insensitive():
    """Test that lowercase Pauli characters work."""
    assert np.allclose(pauli_char_to_matrix('x'), PAULI_X)
    assert np.allclose(pauli_char_to_matrix('y'), PAULI_Y)
    assert np.allclose(pauli_char_to_matrix('z'), PAULI_Z)


def test_pauli_char_to_matrix_invalid():
    """Test that invalid characters raise ValueError."""
    with pytest.raises(ValueError, match="Invalid Pauli character"):
        pauli_char_to_matrix('A')


# =================================================================================================
# Test Pauli String to Matrix
# =================================================================================================

def test_pauli_string_to_matrix_single():
    """Test single-qubit Pauli strings."""
    assert np.allclose(pauli_string_to_matrix('I'), PAULI_I)
    assert np.allclose(pauli_string_to_matrix('X'), PAULI_X)
    assert np.allclose(pauli_string_to_matrix('Y'), PAULI_Y)
    assert np.allclose(pauli_string_to_matrix('Z'), PAULI_Z)


def test_pauli_string_to_matrix_two_qubit():
    """Test two-qubit Pauli strings."""
    # X ⊗ X
    result = pauli_string_to_matrix('XX')
    expected = np.kron(PAULI_X, PAULI_X)
    assert np.allclose(result, expected)

    # Z ⊗ Z
    result = pauli_string_to_matrix('ZZ')
    expected = np.kron(PAULI_Z, PAULI_Z)
    assert np.allclose(result, expected)

    # X ⊗ Y
    result = pauli_string_to_matrix('XY')
    expected = np.kron(PAULI_X, PAULI_Y)
    assert np.allclose(result, expected)


def test_pauli_string_to_matrix_shape():
    """Test that output shape is 2^n × 2^n."""
    assert pauli_string_to_matrix('X').shape == (2, 2)
    assert pauli_string_to_matrix('XX').shape == (4, 4)
    assert pauli_string_to_matrix('XXX').shape == (8, 8)
    assert pauli_string_to_matrix('XXXX').shape == (16, 16)


def test_pauli_string_to_matrix_empty():
    """Test empty string returns 1×1 identity."""
    result = pauli_string_to_matrix('')
    assert result.shape == (1, 1)
    assert result[0, 0] == 1.0


def test_pauli_string_to_matrix_hermitian():
    """Test that all Pauli matrices are Hermitian."""
    for pauli_string in ['X', 'Y', 'Z', 'XX', 'YZ', 'ZXY']:
        mat = pauli_string_to_matrix(pauli_string)
        assert np.allclose(mat, mat.conj().T), f"{pauli_string} not Hermitian"


# =================================================================================================
# Test Pauli Dict to Matrix
# =================================================================================================

def test_pauli_dict_to_matrix_single_term():
    """Test Hamiltonian with single term."""
    pauli_dict = {'ZZ': 0.5}
    result = pauli_dict_to_matrix(pauli_dict, 2)

    expected = 0.5 * pauli_string_to_matrix('ZZ')
    assert np.allclose(result, expected)


def test_pauli_dict_to_matrix_multiple_terms():
    """Test Hamiltonian with multiple terms."""
    pauli_dict = {
        'XX': 0.5,
        'YY': 0.5,
        'ZZ': -1.0,
        'II': 2.0
    }
    result = pauli_dict_to_matrix(pauli_dict, 2)

    expected = (0.5 * pauli_string_to_matrix('XX') +
                0.5 * pauli_string_to_matrix('YY') +
                -1.0 * pauli_string_to_matrix('ZZ') +
                2.0 * pauli_string_to_matrix('II'))
    assert np.allclose(result, expected)


def test_pauli_dict_to_matrix_hermitian():
    """Test that resulting Hamiltonian is Hermitian."""
    pauli_dict = {'XX': 0.5, 'YY': 0.3, 'ZZ': -0.2, 'II': 1.0}
    H = pauli_dict_to_matrix(pauli_dict, 2)

    assert np.allclose(H, H.conj().T)


def test_pauli_dict_to_matrix_shape():
    """Test output shape is correct."""
    assert pauli_dict_to_matrix({'II': 1.0}, 2).shape == (4, 4)
    assert pauli_dict_to_matrix({'III': 1.0}, 3).shape == (8, 8)


def test_pauli_dict_to_matrix_length_mismatch():
    """Test error when Pauli string length doesn't match num_qubits."""
    pauli_dict = {'XXX': 1.0}  # 3 qubits
    with pytest.raises(ValueError, match="has length 3, but num_qubits=2"):
        pauli_dict_to_matrix(pauli_dict, 2)


# =================================================================================================
# Test PauliStringOperator (Matrix-Free)
# =================================================================================================

def test_pauli_operator_init():
    """Test PauliStringOperator initialization."""
    pauli_dict = {'ZZ': 0.5}
    op = PauliStringOperator(pauli_dict, 2)

    assert op.shape == (4, 4)
    assert op.dtype == np.complex128
    assert op.num_qubits == 2
    assert op.dimension == 4


def test_pauli_operator_matvec_identity():
    """Test matvec with identity operator."""
    pauli_dict = {'II': 1.0}
    op = PauliStringOperator(pauli_dict, 2)

    v = np.array([1, 2, 3, 4], dtype=complex)
    result = op @ v

    assert np.allclose(result, v)


def test_pauli_operator_matvec_pauli_z():
    """Test matvec with Pauli Z on single qubit."""
    pauli_dict = {'ZI': 1.0}  # Z on first qubit
    op = PauliStringOperator(pauli_dict, 2)

    # |00⟩
    v = np.array([1, 0, 0, 0], dtype=complex)
    result = op @ v
    expected = np.array([1, 0, 0, 0], dtype=complex)
    assert np.allclose(result, expected)

    # |10⟩ (first qubit is |1⟩)
    v = np.array([0, 1, 0, 0], dtype=complex)
    result = op @ v
    expected = np.array([0, -1, 0, 0], dtype=complex)
    assert np.allclose(result, expected)


def test_pauli_operator_matvec_pauli_x():
    """Test matvec with Pauli X."""
    pauli_dict = {'XI': 1.0}  # X on first qubit
    op = PauliStringOperator(pauli_dict, 2)

    # |00⟩ → |10⟩
    v = np.array([1, 0, 0, 0], dtype=complex)
    result = op @ v
    expected = np.array([0, 1, 0, 0], dtype=complex)
    assert np.allclose(result, expected)


def test_pauli_operator_matvec_pauli_y():
    """Test matvec with Pauli Y."""
    pauli_dict = {'YI': 1.0}  # Y on first qubit
    op = PauliStringOperator(pauli_dict, 2)

    # |00⟩ → -i|10⟩
    v = np.array([1, 0, 0, 0], dtype=complex)
    result = op @ v
    expected = np.array([0, -1j, 0, 0], dtype=complex)
    assert np.allclose(result, expected)


def test_pauli_operator_matches_dense():
    """Test that operator gives same result as dense matrix."""
    pauli_dict = {
        'XX': 0.5,
        'YY': 0.3,
        'ZZ': -0.2,
        'II': 1.0
    }

    # Dense matrix
    H_dense = pauli_dict_to_matrix(pauli_dict, 2)

    # Matrix-free operator
    H_op = PauliStringOperator(pauli_dict, 2)

    # Test on several different states
    test_states = [
        np.array([1, 0, 0, 0], dtype=complex),
        np.array([0, 1, 0, 0], dtype=complex),
        np.array([0, 0, 1, 0], dtype=complex),
        np.array([0, 0, 0, 1], dtype=complex),
        np.array([1, 1, 1, 1], dtype=complex) / 2,
    ]

    for v in test_states:
        result_dense = H_dense @ v
        result_sparse = H_op @ v
        assert np.allclose(result_dense, result_sparse, atol=1e-10)


def test_pauli_operator_multiple_terms():
    """Test operator with multiple Pauli terms."""
    pauli_dict = {'XX': 1.0, 'ZZ': -1.0}
    op = PauliStringOperator(pauli_dict, 2)

    v = np.array([1, 0, 0, 0], dtype=complex)
    result = op @ v

    # XX|00⟩ = |11⟩, ZZ|00⟩ = |00⟩
    # (XX - ZZ)|00⟩ = |11⟩ - |00⟩
    expected = np.array([-1, 0, 0, 1], dtype=complex)
    assert np.allclose(result, expected)


def test_pauli_operator_dimension_mismatch():
    """Test error when vector dimension doesn't match."""
    pauli_dict = {'ZZ': 0.5}
    op = PauliStringOperator(pauli_dict, 2)

    v = np.array([1, 0, 0], dtype=complex)  # Wrong dimension
    with pytest.raises(ValueError, match="dimension mismatch"):
        op @ v


def test_pauli_operator_rmatvec():
    """Test adjoint operation for Hermitian operators."""
    pauli_dict = {'XX': 0.5, 'YY': 0.3}
    op = PauliStringOperator(pauli_dict, 2)

    v = np.array([1, 0, 0, 0], dtype=complex)

    # For Hermitian operators, H† = H
    result_forward = op @ v
    result_adjoint = op._rmatvec(v)

    assert np.allclose(result_forward, result_adjoint)


# =================================================================================================
# Test Operator Type Selection
# =================================================================================================

def test_get_operator_type_small():
    """Test that small systems use dense representation with 1 GB threshold."""
    # 10 qubits: 2^10 = 1024, memory = 1024^2 * 16 bytes = 16 MB << 1 GB
    assert get_operator_type(10, memory_threshold_gb=1.0) == "dense"
    # 13 qubits: 2^13 = 8192, memory = 8192^2 * 16 bytes = 1 GB
    assert get_operator_type(13, memory_threshold_gb=1.0) == "dense"


def test_get_operator_type_large():
    """Test that large systems use sparse representation with 1 GB threshold."""
    # 14 qubits: 2^14 = 16384, memory = 16384^2 * 16 bytes = 4 GB > 1 GB
    assert get_operator_type(14, memory_threshold_gb=1.0) == "sparse"
    # 20 qubits: much larger
    assert get_operator_type(20, memory_threshold_gb=1.0) == "sparse"


def test_get_operator_type_custom_threshold():
    """Test with different memory thresholds."""
    # With 0.5 GB threshold: 12 qubits = 256 MB, should use dense
    assert get_operator_type(12, memory_threshold_gb=0.5) == "dense"
    # With 0.5 GB threshold: 13 qubits = 1 GB, should use sparse
    assert get_operator_type(13, memory_threshold_gb=0.5) == "sparse"
    # With 0.1 GB threshold: 11 qubits = 64 MB, should use dense
    assert get_operator_type(11, memory_threshold_gb=0.1) == "dense"
    # With 0.1 GB threshold: 12 qubits = 256 MB, should use sparse
    assert get_operator_type(12, memory_threshold_gb=0.1) == "sparse"


def test_get_operator_type_force_dense():
    """Test forcing dense representation."""
    assert get_operator_type(20, memory_threshold_gb=1.0, force_dense=True) == "dense"


def test_get_operator_type_force_sparse():
    """Test forcing sparse representation."""
    assert get_operator_type(10, memory_threshold_gb=1.0, force_sparse=True) == "sparse"


def test_get_operator_type_conflicting_forces():
    """Test error when forcing both dense and sparse."""
    with pytest.raises(ValueError, match="Cannot force both"):
        get_operator_type(10, memory_threshold_gb=1.0, force_dense=True, force_sparse=True)


def test_get_operator_type_very_large():
    """Test log-space calculation for very large qubit counts."""
    # 50 qubits: 2^50 dimension, way too large for dense
    assert get_operator_type(50, memory_threshold_gb=1.0) == "sparse"
    # 40 qubits: 2^40 dimension, would require ~16 PB
    assert get_operator_type(40, memory_threshold_gb=1.0) == "sparse"
    # Even with huge threshold (1 PB), 50 qubits should be sparse
    assert get_operator_type(50, memory_threshold_gb=1e6) == "sparse"


# =================================================================================================
# Test Create Hamiltonian Operator
# =================================================================================================

def test_create_hamiltonian_operator_dense():
    """Test creating dense operator for small system."""
    pauli_dict = {'XX': 1.0, 'ZZ': -1.0}
    op = create_hamiltonian_operator(pauli_dict, 2, memory_threshold_gb=1.0)

    assert isinstance(op, np.ndarray)
    assert op.shape == (4, 4)


def test_create_hamiltonian_operator_sparse():
    """Test creating sparse operator for large system."""
    pauli_dict = {'X' * 20: 1.0}
    op = create_hamiltonian_operator(pauli_dict, 20, memory_threshold_gb=1.0)

    assert isinstance(op, PauliStringOperator)
    assert op.shape == (2**20, 2**20)


def test_create_hamiltonian_operator_forced_dense():
    """Test forcing dense for a system that would normally be sparse."""
    # Use 12 qubits with 0.1 GB threshold
    # 12 qubits: 4096^2 * 16 bytes = 268 MB > 0.1 GB threshold (so normally sparse)
    # But 268 MB < 1 GB (safe to allocate)
    pauli_dict = {'X' * 12: 1.0}
    op = create_hamiltonian_operator(pauli_dict, 12, memory_threshold_gb=0.1, force_dense=True)

    assert isinstance(op, np.ndarray)
    assert op.shape == (2**12, 2**12)


def test_create_hamiltonian_operator_forced_sparse():
    """Test forcing sparse for small system."""
    pauli_dict = {'XX': 1.0}
    op = create_hamiltonian_operator(pauli_dict, 2, memory_threshold_gb=1.0, force_sparse=True)

    assert isinstance(op, PauliStringOperator)


# =================================================================================================
# Integration Tests
# =================================================================================================

def test_known_hamiltonian_hydrogen():
    """Test with a known Hamiltonian: hydrogen-like."""
    # Simple 2-qubit Hamiltonian
    pauli_dict = {
        'ZI': -0.5,
        'IZ': -0.5,
        'ZZ': 0.25,
        'XX': 0.5,
        'II': 1.0
    }

    H = pauli_dict_to_matrix(pauli_dict, 2)

    # Check Hermiticity
    assert np.allclose(H, H.conj().T)

    # Check eigenvalues are real (Hermitian property)
    eigenvalues = np.linalg.eigvalsh(H)
    assert all(np.isreal(eigenvalues))


def test_large_coefficient():
    """Test with large coefficients."""
    pauli_dict = {'XX': 1e10, 'ZZ': -1e10}
    H = pauli_dict_to_matrix(pauli_dict, 2)

    assert not np.isnan(H).any()
    assert not np.isinf(H).any()


def test_complex_coefficients_raises_warning():
    """Test that complex coefficients work but may not be physical."""
    # Note: For physical Hamiltonians, coefficients should be real
    # But mathematically we can have complex coefficients
    pauli_dict = {'XX': 1.0 + 0.5j}
    H = pauli_dict_to_matrix(pauli_dict, 2)

    # Should not be Hermitian
    assert not np.allclose(H, H.conj().T)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
