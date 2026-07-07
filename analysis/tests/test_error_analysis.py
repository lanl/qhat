"""
Tests for error analysis functionality.

Tests cover:
- Eigenvalue error computation
- Matrix norm errors (Frobenius and spectral)
- State-dependent errors
- Integration with eigendecomposition
- File I/O for error results
"""

import numpy as np
import pytest
import tempfile
import os
from pathlib import Path

from qhat.analysis.analysis import error_analysis, validate_and_autocomplete_analysis_config
from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.file_io import save_state


# =================================================================================================
# Unit Tests: Eigenvalue Errors
# =================================================================================================

def test_eigenvalue_error_zero_when_identical():
    """Test that eigenvalue error is zero when matrices are identical."""
    from qhat.analysis.file_io import save_eigendecomposition, load_eigendecomposition

    # Create identical eigendecompositions
    eigenvalues = np.array([1.0, 2.0, 3.0])
    eigenvectors = np.eye(3, dtype=complex)

    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Save both decompositions
            save_eigendecomposition(
                'exact_eigendecomposition.npz', eigenvalues, eigenvectors,
                matrix_type='exact', num_eigenvalues=3, which_eigenvalues='smallest'
            )
            save_eigendecomposition(
                'approximate_eigendecomposition.npz', eigenvalues, eigenvectors,
                matrix_type='approximate', num_eigenvalues=3, which_eigenvalues='smallest'
            )

            # Load eigendecompositions
            exact_eigendecomp = load_eigendecomposition('exact_eigendecomposition.npz')
            approx_eigendecomp = load_eigendecomposition('approximate_eigendecomposition.npz')

            # Compute errors
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=None,
                unitary_matrix=None,
                exact_eigendecomp=exact_eigendecomp,
                approx_eigendecomp=approx_eigendecomp
            )

            # Verify zero errors
            assert 'eigenvalue_errors' in results
            np.testing.assert_array_almost_equal(
                results['eigenvalue_errors']['absolute_errors'],
                [0.0, 0.0, 0.0]
            )

        finally:
            os.chdir(original_dir)


def test_eigenvalue_error_nonzero_when_different():
    """Test that eigenvalue error is nonzero when matrices differ."""
    from qhat.analysis.file_io import save_eigendecomposition, load_eigendecomposition

    # Create different eigendecompositions
    exact_eigenvalues = np.array([1.0, 2.0, 3.0])
    approx_eigenvalues = np.array([1.1, 1.9, 3.2])
    eigenvectors = np.eye(3, dtype=complex)

    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Save both decompositions
            save_eigendecomposition(
                'exact_eigendecomposition.npz', exact_eigenvalues, eigenvectors,
                matrix_type='exact', num_eigenvalues=3, which_eigenvalues='smallest'
            )
            save_eigendecomposition(
                'approximate_eigendecomposition.npz', approx_eigenvalues, eigenvectors,
                matrix_type='approximate', num_eigenvalues=3, which_eigenvalues='smallest'
            )

            # Load eigendecompositions
            exact_eigendecomp = load_eigendecomposition('exact_eigendecomposition.npz')
            approx_eigendecomp = load_eigendecomposition('approximate_eigendecomposition.npz')

            # Compute errors
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=None,
                unitary_matrix=None,
                exact_eigendecomp=exact_eigendecomp,
                approx_eigendecomp=approx_eigendecomp
            )

            # Verify nonzero errors
            assert 'eigenvalue_errors' in results
            expected_errors = exact_eigenvalues - approx_eigenvalues
            np.testing.assert_array_almost_equal(
                results['eigenvalue_errors']['absolute_errors'],
                expected_errors.tolist()
            )

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Unit Tests: Matrix Norm Errors (Dense)
# =================================================================================================

def test_frobenius_norm_zero_when_identical():
    """Test Frobenius norm error is zero for identical matrices."""
    # 2x2 identity
    matrix = np.eye(2, dtype=complex)

    config = AnalysisConfiguration()
    config.error_matrix_norms = 'frobenius'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=matrix,
                unitary_matrix=matrix
            )

            assert 'matrix_frobenius_error' in results
            assert results['matrix_frobenius_error'] < 1e-10

        finally:
            os.chdir(original_dir)


def test_frobenius_norm_nonzero_when_different():
    """Test Frobenius norm error is nonzero for different matrices."""
    exact = np.eye(2, dtype=complex)
    approx = np.array([[1.0, 0.1], [0.0, 1.0]], dtype=complex)

    config = AnalysisConfiguration()
    config.error_matrix_norms = 'frobenius'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=exact,
                unitary_matrix=approx
            )

            assert 'matrix_frobenius_error' in results
            # Manual calculation: ||[[0, 0.1], [0, 0]]||_F = sqrt(0.01) = 0.1
            expected = np.linalg.norm(exact - approx, 'fro')
            np.testing.assert_almost_equal(
                results['matrix_frobenius_error'],
                expected
            )

        finally:
            os.chdir(original_dir)


def test_spectral_norm_zero_when_identical():
    """Test spectral norm error is zero for identical matrices."""
    matrix = np.eye(2, dtype=complex)

    config = AnalysisConfiguration()
    config.error_matrix_norms = 'spectral'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=matrix,
                unitary_matrix=matrix
            )

            assert 'matrix_spectral_error' in results
            assert results['matrix_spectral_error'] < 1e-10

        finally:
            os.chdir(original_dir)


def test_spectral_norm_nonzero_when_different():
    """Test spectral norm error is nonzero for different matrices."""
    exact = np.eye(2, dtype=complex)
    approx = np.array([[1.0, 0.3], [0.0, 1.0]], dtype=complex)

    config = AnalysisConfiguration()
    config.error_matrix_norms = 'spectral'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=exact,
                unitary_matrix=approx
            )

            assert 'matrix_spectral_error' in results
            # Manual calculation: ||[[0, 0.3], [0, 0]]||_2 = 0.3
            expected = np.linalg.norm(exact - approx, 2)
            np.testing.assert_almost_equal(
                results['matrix_spectral_error'],
                expected,
                decimal=5
            )

        finally:
            os.chdir(original_dir)


def test_both_matrix_norms():
    """Test computing both Frobenius and spectral norms."""
    exact = np.eye(2, dtype=complex)
    approx = np.array([[1.0, 0.2], [0.0, 1.0]], dtype=complex)

    config = AnalysisConfiguration()
    config.error_matrix_norms = ['frobenius', 'spectral']

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=exact,
                unitary_matrix=approx
            )

            # Both should be present
            assert 'matrix_frobenius_error' in results
            assert 'matrix_spectral_error' in results

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Unit Tests: State-Dependent Errors
# =================================================================================================

def test_state_error_zero_when_identical():
    """Test state error is zero for identical operators."""
    matrix = np.eye(2, dtype=complex)
    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.error_state_inputs = 'test_state.npy'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Save state
            save_state('test_state.npy', state)

            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=matrix,
                unitary_matrix=matrix
            )

            assert 'state_errors' in results
            assert len(results['state_errors']) == 1
            assert results['state_errors'][0]['absolute_error'] < 1e-10

        finally:
            os.chdir(original_dir)


def test_state_error_nonzero_when_different():
    """Test state error is nonzero for different operators."""
    exact = np.eye(2, dtype=complex)
    approx = np.array([[1.0, 0.1], [0.0, 1.0]], dtype=complex)
    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.error_state_inputs = 'test_state.npy'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Save state
            save_state('test_state.npy', state)

            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=exact,
                unitary_matrix=approx
            )

            assert 'state_errors' in results
            assert len(results['state_errors']) == 1

            # Manual calculation
            exact_final = exact @ state  # [1.0, 0.0]
            approx_final = approx @ state  # [1.0, 0.0]
            expected_error = np.linalg.norm(exact_final - approx_final)

            np.testing.assert_almost_equal(
                results['state_errors'][0]['absolute_error'],
                expected_error
            )

        finally:
            os.chdir(original_dir)


def test_state_error_multiple_states():
    """Test state errors for multiple input states."""
    matrix = np.eye(2, dtype=complex)
    state1 = np.array([1.0, 0.0], dtype=complex)
    state2 = np.array([0.0, 1.0], dtype=complex)

    config = AnalysisConfiguration()
    config.error_state_inputs = ['state1.npy', 'state2.npy']

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Save states
            save_state('state1.npy', state1)
            save_state('state2.npy', state2)

            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=matrix,
                unitary_matrix=matrix
            )

            assert 'state_errors' in results
            assert len(results['state_errors']) == 2
            assert results['state_errors'][0]['input_file'] == 'state1.npy'
            assert results['state_errors'][1]['input_file'] == 'state2.npy'

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Integration Tests: All Error Types Together
# =================================================================================================

def test_all_error_types_together():
    """Test computing all three error types in one analysis."""
    from qhat.analysis.file_io import save_eigendecomposition, load_eigendecomposition

    # Setup matrices
    exact = np.diag([1.0, 2.0])
    approx = np.diag([1.1, 1.9])
    state = np.array([1.0, 0.0], dtype=complex)
    eigenvalues_exact = np.array([1.0, 2.0])
    eigenvalues_approx = np.array([1.1, 1.9])
    eigenvectors = np.eye(2, dtype=complex)

    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    config.num_eigenvalues = 2  # Required for eigenvalue errors
    config.eigendecomposition_matrices = 'both'  # Will be auto-set by validation, but explicit here
    config.error_matrix_norms = ['frobenius', 'spectral']
    config.error_state_inputs = 'test_state.npy'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Save eigendecompositions
            save_eigendecomposition(
                'exact_eigendecomposition.npz', eigenvalues_exact, eigenvectors,
                matrix_type='exact', num_eigenvalues=2, which_eigenvalues='smallest'
            )
            save_eigendecomposition(
                'approximate_eigendecomposition.npz', eigenvalues_approx, eigenvectors,
                matrix_type='approximate', num_eigenvalues=2, which_eigenvalues='smallest'
            )

            # Load eigendecompositions
            exact_eigendecomp = load_eigendecomposition('exact_eigendecomposition.npz')
            approx_eigendecomp = load_eigendecomposition('approximate_eigendecomposition.npz')

            # Save state
            save_state('test_state.npy', state)

            # Compute all errors
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=exact,
                unitary_matrix=approx,
                exact_eigendecomp=exact_eigendecomp,
                approx_eigendecomp=approx_eigendecomp
            )

            # Verify all three types present
            assert 'eigenvalue_errors' in results
            assert 'matrix_frobenius_error' in results
            assert 'matrix_spectral_error' in results
            assert 'state_errors' in results

        finally:
            os.chdir(original_dir)


def test_error_output_file_created():
    """Test that error_analysis.npz file is created."""
    from qhat.analysis.file_io import save_eigendecomposition, load_eigendecomposition

    eigenvalues = np.array([1.0, 2.0])
    eigenvectors = np.eye(2, dtype=complex)

    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Save eigendecompositions
            save_eigendecomposition(
                'exact_eigendecomposition.npz', eigenvalues, eigenvectors,
                matrix_type='exact', num_eigenvalues=2, which_eigenvalues='smallest'
            )
            save_eigendecomposition(
                'approximate_eigendecomposition.npz', eigenvalues, eigenvectors,
                matrix_type='approximate', num_eigenvalues=2, which_eigenvalues='smallest'
            )

            # Load eigendecompositions
            exact_eigendecomp = load_eigendecomposition('exact_eigendecomposition.npz')
            approx_eigendecomp = load_eigendecomposition('approximate_eigendecomposition.npz')

            # Compute errors
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=None,
                unitary_matrix=None,
                exact_eigendecomp=exact_eigendecomp,
                approx_eigendecomp=approx_eigendecomp
            )

            # Verify file created
            assert 'output_file' in results
            assert os.path.exists('error_analysis.npz')

            # Verify can load
            data = np.load('error_analysis.npz')
            assert 'eigenvalue_absolute_errors' in data
            assert 'eigenvalue_relative_errors' in data

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Error Handling Tests
# =================================================================================================

def test_invalid_matrix_norm_type():
    """Test error when requesting invalid matrix norm type."""
    matrix = np.eye(2, dtype=complex)

    config = AnalysisConfiguration()
    config.error_matrix_norms = 'invalid'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with pytest.raises(ValueError, match="Unknown matrix norm type"):
        with tempfile.TemporaryDirectory() as tmpdir:
            original_dir = os.getcwd()
            os.chdir(tmpdir)
            try:
                error_analysis(
                    config,
                    hamiltonian=None,
                    algorithm=None,
                    exact_matrix=matrix,
                    unitary_matrix=matrix
                )
            finally:
                os.chdir(original_dir)


def test_missing_state_file():
    """Test error when state file doesn't exist."""
    matrix = np.eye(2, dtype=complex)

    config = AnalysisConfiguration()
    config.error_state_inputs = 'nonexistent.npy'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with pytest.raises(Exception):  # Will be FileNotFoundError or similar
        with tempfile.TemporaryDirectory() as tmpdir:
            original_dir = os.getcwd()
            os.chdir(tmpdir)
            try:
                error_analysis(
                    config,
                    hamiltonian=None,
                    algorithm=None,
                    exact_matrix=matrix,
                    unitary_matrix=matrix
                )
            finally:
                os.chdir(original_dir)


# =================================================================================================
# Matrix-Free Operator Tests
# =================================================================================================

def test_state_error_with_matrix_free_operators():
    """Test state errors work with matrix-free operators."""
    from qhat.analysis.matrix_operations import PauliStringOperator

    # Create matrix-free operators
    pauli_dict = {'II': 1.0}  # Identity operator
    exact_op = PauliStringOperator(pauli_dict, num_qubits=2)
    approx_op = PauliStringOperator(pauli_dict, num_qubits=2)

    state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.error_state_inputs = 'test_state.npy'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Save state
            save_state('test_state.npy', state)

            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=exact_op,
                unitary_matrix=approx_op
            )

            assert 'state_errors' in results
            # Error should be zero for identical operators
            assert results['state_errors'][0]['absolute_error'] < 1e-10

        finally:
            os.chdir(original_dir)


def test_frobenius_norm_with_matrix_free_small():
    """Test Frobenius norm computation with matrix-free operators (small system)."""
    from qhat.analysis.matrix_operations import PauliStringOperator

    # 2-qubit system (dimension 4) - still small enough for both paths
    pauli_dict_exact = {'II': 1.0}
    pauli_dict_approx = {'II': 1.0}

    exact_op = PauliStringOperator(pauli_dict_exact, num_qubits=2)
    approx_op = PauliStringOperator(pauli_dict_approx, num_qubits=2)

    config = AnalysisConfiguration()
    config.error_matrix_norms = 'frobenius'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=exact_op,
                unitary_matrix=approx_op
            )

            assert 'matrix_frobenius_error' in results
            # Identical operators should have zero error
            assert results['matrix_frobenius_error'] < 1e-10

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Relative Error Tests
# =================================================================================================

def test_eigenvalue_relative_error():
    """Test that relative eigenvalue errors are computed correctly."""
    from qhat.analysis.file_io import save_eigendecomposition, load_eigendecomposition

    exact_eigenvalues = np.array([10.0, 20.0])
    approx_eigenvalues = np.array([11.0, 18.0])
    eigenvectors = np.eye(2, dtype=complex)

    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            save_eigendecomposition(
                'exact_eigendecomposition.npz', exact_eigenvalues, eigenvectors,
                matrix_type='exact', num_eigenvalues=2, which_eigenvalues='smallest'
            )
            save_eigendecomposition(
                'approximate_eigendecomposition.npz', approx_eigenvalues, eigenvectors,
                matrix_type='approximate', num_eigenvalues=2, which_eigenvalues='smallest'
            )

            # Load eigendecompositions
            exact_eigendecomp = load_eigendecomposition('exact_eigendecomposition.npz')
            approx_eigendecomp = load_eigendecomposition('approximate_eigendecomposition.npz')

            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=None,
                unitary_matrix=None,
                exact_eigendecomp=exact_eigendecomp,
                approx_eigendecomp=approx_eigendecomp
            )

            # Check relative errors
            # (10 - 11) / |10| = -0.1, (20 - 18) / |20| = 0.1
            expected_relative = np.array([-0.1, 0.1])
            np.testing.assert_array_almost_equal(
                results['eigenvalue_errors']['relative_errors'],
                expected_relative.tolist()
            )

        finally:
            os.chdir(original_dir)


def test_state_relative_error():
    """Test that relative state errors are computed correctly."""
    exact = np.diag([2.0, 3.0])
    approx = np.diag([2.1, 2.9])
    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.error_state_inputs = 'test_state.npy'
    validate_and_autocomplete_analysis_config(config)  # Normalize config values

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            save_state('test_state.npy', state)

            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_matrix=exact,
                unitary_matrix=approx
            )

            # exact @ state = [2.0, 0.0], norm = 2.0
            # approx @ state = [2.1, 0.0], norm = 2.1
            # absolute error = 0.1
            # relative error = 0.1 / 2.0 = 0.05
            assert 'state_errors' in results
            assert results['state_errors'][0]['relative_error'] == pytest.approx(0.05)

        finally:
            os.chdir(original_dir)
