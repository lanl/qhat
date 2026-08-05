"""
Tests for error analysis functionality.

Tests cover:
- Eigenvalue error computation
- Matrix norm errors (Frobenius and spectral)
- State-dependent errors
- Integration with OperatorRepresentation
- File I/O for error results
"""

import numpy as np
import pytest
import tempfile
import os
import scipy.linalg

from qhat.analysis.error_analysis import error_analysis
from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.file_io import save_state
from qhat.analysis.operators import OperatorRepresentation


# =================================================================================================
# Helper Functions
# =================================================================================================

def create_test_hamiltonian():
    """Create a simple test Hamiltonian (Pauli Z)."""
    return np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


def create_identical_operator_representations(t=1.0):
    """Create exact_op and approx_op representing the same operator.

    Returns:
        tuple: (exact_op, approx_op) both OperatorRepresentation objects
    """
    H = create_test_hamiltonian()

    exact_op = OperatorRepresentation(
        data=H,
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='dense_matrix',
        tevol_hbar=t
    )

    approx_op = OperatorRepresentation(
        data=H,
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='dense_matrix',
        tevol_hbar=t
    )

    return exact_op, approx_op


def create_different_operator_representations(t=1.0):
    """Create exact_op and approx_op that are different.

    Returns:
        tuple: (exact_op, approx_op) both OperatorRepresentation objects
    """
    H_exact = create_test_hamiltonian()
    H_approx = np.array([[1.1, 0.0], [0.0, -0.9]], dtype=complex)

    exact_op = OperatorRepresentation(
        data=H_exact,
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='dense_matrix',
        tevol_hbar=t
    )

    approx_op = OperatorRepresentation(
        data=H_approx,
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='dense_matrix',
        tevol_hbar=t
    )

    return exact_op, approx_op


# =================================================================================================
# Unit Tests: Eigenvalue Errors
# =================================================================================================

def test_eigenenergy_error_zero_when_identical():
    """Test that eigenenergy error is zero when matrices are identical."""
    # Create identical eigendecompositions
    eigenvalues = np.array([1.0, 2.0, 3.0])
    eigenvectors = np.eye(3, dtype=complex)

    eigendata = {'eigenvalues': eigenvalues, 'eigenvectors': eigenvectors}

    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Create OperatorRepresentation objects from eigendecompositions
            exact_op = OperatorRepresentation(
                data=eigendata,
                operator_type='hamiltonian',
                energy_shifted=False,
                representation='eigendecomposition',
                tevol_hbar=1.0
            )

            approx_op = OperatorRepresentation(
                data=eigendata,
                operator_type='hamiltonian',
                energy_shifted=False,
                representation='eigendecomposition',
                tevol_hbar=1.0
            )

            # Compute errors
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=1.0
            )

            # Verify zero errors
            assert 'eigenenergy_errors' in results
            np.testing.assert_array_almost_equal(
                results['eigenenergy_errors']['absolute_errors'],
                [0.0, 0.0, 0.0]
            )

        finally:
            os.chdir(original_dir)


def test_eigenenergy_error_nonzero_when_different():
    """Test that eigenenergy error is nonzero when matrices differ."""
    # Create different eigendecompositions
    exact_eigenvalues = np.array([1.0, 2.0, 3.0])
    approx_eigenvalues = np.array([1.1, 1.9, 3.2])
    eigenvectors = np.eye(3, dtype=complex)

    exact_eigendata = {'eigenvalues': exact_eigenvalues, 'eigenvectors': eigenvectors}
    approx_eigendata = {'eigenvalues': approx_eigenvalues, 'eigenvectors': eigenvectors}

    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Create OperatorRepresentation objects from eigendecompositions
            exact_op = OperatorRepresentation(
                data=exact_eigendata,
                operator_type='hamiltonian',
                energy_shifted=False,
                representation='eigendecomposition',
                tevol_hbar=1.0
            )

            approx_op = OperatorRepresentation(
                data=approx_eigendata,
                operator_type='hamiltonian',
                energy_shifted=False,
                representation='eigendecomposition',
                tevol_hbar=1.0
            )

            # Compute errors
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=1.0
            )

            # Verify nonzero errors
            assert 'eigenenergy_errors' in results
            # New method computes errors via U eigenvalue ratios, which should give
            # approximately the same magnitude but may differ in sign/ordering due to
            # eigenvector matching. Check that errors are nonzero and reasonable magnitude.
            errors = np.array(results['eigenenergy_errors']['absolute_errors'])
            expected_errors = exact_eigenvalues - approx_eigenvalues

            # Check magnitudes are correct (allowing for reordering)
            np.testing.assert_array_almost_equal(
                np.sort(np.abs(errors)),
                np.sort(np.abs(expected_errors))
            )

        finally:
            os.chdir(original_dir)


def test_eigenvalue_relative_error():
    """Test that relative eigenvalue errors are computed correctly."""
    exact_eigenvalues = np.array([10.0, 20.0])
    approx_eigenvalues = np.array([11.0, 18.0])
    eigenvectors = np.eye(2, dtype=complex)

    exact_eigendata = {'eigenvalues': exact_eigenvalues, 'eigenvectors': eigenvectors}
    approx_eigendata = {'eigenvalues': approx_eigenvalues, 'eigenvectors': eigenvectors}

    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Create OperatorRepresentation objects from eigendecompositions
            exact_op = OperatorRepresentation(
                data=exact_eigendata,
                operator_type='hamiltonian',
                energy_shifted=False,
                representation='eigendecomposition',
                tevol_hbar=1.0
            )

            approx_op = OperatorRepresentation(
                data=approx_eigendata,
                operator_type='hamiltonian',
                energy_shifted=False,
                representation='eigendecomposition',
                tevol_hbar=1.0
            )

            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=1.0
            )

            # Check relative errors: should be approximately 10% for both eigenvalues
            # New method uses U eigenvalue ratios which may reorder results
            relative_errors = np.array(results['eigenenergy_errors']['relative_errors'])
            expected_relative = np.array([-0.1, 0.1])

            # Check magnitudes are correct (allowing for reordering)
            np.testing.assert_array_almost_equal(
                np.sort(np.abs(relative_errors)),
                np.sort(np.abs(expected_relative))
            )

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Unit Tests: Matrix Norm Errors (Dense)
# =================================================================================================

def test_frobenius_norm_zero_when_identical():
    """Test Frobenius norm error is zero for identical operators."""
    t = 1.0
    exact_op, approx_op = create_identical_operator_representations(t)

    config = AnalysisConfiguration()
    config.error_matrix_norms = 'frobenius'

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=t,
                energy_shift=0.0
            )

            assert 'matrix_frobenius_error' in results
            assert results['matrix_frobenius_error'] < 1e-10

        finally:
            os.chdir(original_dir)


def test_frobenius_norm_nonzero_when_different():
    """Test Frobenius norm error is nonzero for different operators."""
    t = 1.0
    exact_op, approx_op = create_different_operator_representations(t)

    config = AnalysisConfiguration()
    config.error_matrix_norms = 'frobenius'

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=t,
                energy_shift=0.0
            )

            assert 'matrix_frobenius_error' in results
            # Error should be nonzero since operators are different
            assert results['matrix_frobenius_error'] > 1e-10

        finally:
            os.chdir(original_dir)


def test_spectral_norm_zero_when_identical():
    """Test spectral norm error is zero for identical operators."""
    t = 1.0
    exact_op, approx_op = create_identical_operator_representations(t)

    config = AnalysisConfiguration()
    config.error_matrix_norms = 'spectral'

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=t,
                energy_shift=0.0
            )

            assert 'matrix_spectral_error' in results
            assert results['matrix_spectral_error'] < 1e-10

        finally:
            os.chdir(original_dir)


def test_spectral_norm_nonzero_when_different():
    """Test spectral norm error is nonzero for different operators."""
    t = 1.0
    exact_op, approx_op = create_different_operator_representations(t)

    config = AnalysisConfiguration()
    config.error_matrix_norms = 'spectral'

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=t,
                energy_shift=0.0
            )

            assert 'matrix_spectral_error' in results
            # Error should be nonzero since operators are different
            assert results['matrix_spectral_error'] > 1e-10

        finally:
            os.chdir(original_dir)


def test_both_matrix_norms():
    """Test computing both Frobenius and spectral norms."""
    t = 1.0
    exact_op, approx_op = create_different_operator_representations(t)

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
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=t,
                energy_shift=0.0
            )

            # Both should be present
            assert 'matrix_frobenius_error' in results
            assert 'matrix_spectral_error' in results

        finally:
            os.chdir(original_dir)


def test_invalid_matrix_norm_type():
    """Test error when requesting invalid matrix norm type."""
    t = 1.0
    exact_op, approx_op = create_identical_operator_representations(t)

    config = AnalysisConfiguration()
    config.error_matrix_norms = 'invalid'

    with pytest.raises(ValueError, match="Unknown matrix norm type"):
        with tempfile.TemporaryDirectory() as tmpdir:
            original_dir = os.getcwd()
            os.chdir(tmpdir)
            try:
                error_analysis(
                    config,
                    hamiltonian=None,
                    algorithm=None,
                    exact_op=exact_op,
                    approx_op=approx_op,
                    timestep=t,
                    energy_shift=0.0
                )
            finally:
                os.chdir(original_dir)


# =================================================================================================
# Unit Tests: State-Dependent Errors
# =================================================================================================

def test_state_error_zero_when_identical():
    """Test state error is zero for identical operators."""
    t = 1.0
    exact_op, approx_op = create_identical_operator_representations(t)
    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.error_state_inputs = 'test_state.npy'

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
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=t,
                energy_shift=0.0
            )

            assert 'state_errors' in results
            assert len(results['state_errors']) == 1
            assert results['state_errors'][0]['absolute_error'] < 1e-10

        finally:
            os.chdir(original_dir)


def test_state_error_nonzero_when_different():
    """Test state error is nonzero for different operators."""
    t = 1.0
    exact_op, approx_op = create_different_operator_representations(t)
    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.error_state_inputs = 'test_state.npy'

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
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=t,
                energy_shift=0.0
            )

            assert 'state_errors' in results
            assert len(results['state_errors']) == 1

            # Error should be nonzero since operators are different
            assert results['state_errors'][0]['absolute_error'] > 1e-10

        finally:
            os.chdir(original_dir)


def test_state_error_multiple_states():
    """Test state errors for multiple input states."""
    t = 1.0
    exact_op, approx_op = create_identical_operator_representations(t)
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
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=t,
                energy_shift=0.0
            )

            assert 'state_errors' in results
            assert len(results['state_errors']) == 2
            assert results['state_errors'][0]['input_file'] == 'state1.npy'
            assert results['state_errors'][1]['input_file'] == 'state2.npy'

        finally:
            os.chdir(original_dir)


def test_state_relative_error():
    """Test that relative state errors are computed correctly."""
    # Create two diagonal Hamiltonians that give different time evolutions
    t = 1.0
    H_exact = np.diag([2.0, 3.0])
    H_approx = np.diag([2.1, 2.9])
    state = np.array([1.0, 0.0], dtype=complex)

    exact_op = OperatorRepresentation(
        data=H_exact,
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='dense_matrix',
        tevol_hbar=t
    )

    approx_op = OperatorRepresentation(
        data=H_approx,
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='dense_matrix',
        tevol_hbar=t
    )

    config = AnalysisConfiguration()
    config.error_state_inputs = 'test_state.npy'

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            save_state('test_state.npy', state)

            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=t,
                energy_shift=0.0
            )

            # Check that error is computed (specific value depends on the time evolution)
            assert 'state_errors' in results
            assert results['state_errors'][0]['relative_error'] > 0

        finally:
            os.chdir(original_dir)


def test_missing_state_file():
    """Test error when state file doesn't exist."""
    matrix = np.eye(2, dtype=complex)

    exact_op = OperatorRepresentation(
        data=matrix,
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='dense_matrix',
        tevol_hbar=1.0
    )

    approx_op = OperatorRepresentation(
        data=matrix,
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='dense_matrix',
        tevol_hbar=1.0
    )

    config = AnalysisConfiguration()
    config.error_state_inputs = 'nonexistent.npy'

    with pytest.raises(Exception):  # Will be FileNotFoundError or similar
        with tempfile.TemporaryDirectory() as tmpdir:
            original_dir = os.getcwd()
            os.chdir(tmpdir)
            try:
                error_analysis(
                    config,
                    hamiltonian=None,
                    algorithm=None,
                    exact_op=exact_op,
                    approx_op=approx_op,
                    timestep=1.0
                )
            finally:
                os.chdir(original_dir)


# =================================================================================================
# Integration Tests: All Error Types Together
# =================================================================================================

def test_all_error_types_together():
    """Test computing all three error types in one analysis."""
    # Setup matrices
    t = 1.0
    H_exact = np.diag([1.0, 2.0])
    H_approx = np.diag([1.1, 1.9])
    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    config.error_matrix_norms = ['frobenius', 'spectral']
    config.error_state_inputs = 'test_state.npy'

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Create OperatorRepresentation objects
            exact_op = OperatorRepresentation(
                data=H_exact,
                operator_type='hamiltonian',
                energy_shifted=False,
                representation='dense_matrix',
                tevol_hbar=t
            )

            approx_op = OperatorRepresentation(
                data=H_approx,
                operator_type='hamiltonian',
                energy_shifted=False,
                representation='dense_matrix',
                tevol_hbar=t
            )

            # Save state
            save_state('test_state.npy', state)

            # Compute all errors
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=t,
                energy_shift=0.0
            )

            # Verify all three types present
            assert 'eigenenergy_errors' in results
            assert 'matrix_frobenius_error' in results
            assert 'matrix_spectral_error' in results
            assert 'state_errors' in results

        finally:
            os.chdir(original_dir)


def test_error_output_file_created():
    """Test that error_analysis.npz file is created."""
    eigenvalues = np.array([1.0, 2.0])
    eigenvectors = np.eye(2, dtype=complex)
    eigendata = {'eigenvalues': eigenvalues, 'eigenvectors': eigenvectors}

    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Create OperatorRepresentation objects from eigendecompositions
            exact_op = OperatorRepresentation(
                data=eigendata,
                operator_type='hamiltonian',
                energy_shifted=False,
                representation='eigendecomposition',
                tevol_hbar=1.0
            )

            approx_op = OperatorRepresentation(
                data=eigendata,
                operator_type='hamiltonian',
                energy_shifted=False,
                representation='eigendecomposition',
                tevol_hbar=1.0
            )

            # Compute errors
            results = error_analysis(
                config,
                hamiltonian=None,
                algorithm=None,
                exact_op=exact_op,
                approx_op=approx_op,
                timestep=1.0
            )

            # Verify file created
            assert 'output_file' in results
            assert os.path.exists('error_analysis.npz')

            # Verify can load
            data = np.load('error_analysis.npz')
            assert 'eigenenergy_absolute_errors' in data
            assert 'eigenenergy_relative_errors' in data

        finally:
            os.chdir(original_dir)
