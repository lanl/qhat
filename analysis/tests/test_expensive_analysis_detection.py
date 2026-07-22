"""
Tests for the shared functions that detect when expensive analyses are required.

These functions are used by both validation and analyze_algorithm to ensure
consistent detection of what computations are needed.
"""

import pytest
from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.analysis import (
    requires_exact_eigendecomposition,
    requires_approximate_eigendecomposition,
    requires_exact_matrix,
    requires_approximate_matrix
)


# =================================================================================================
# Tests for requires_exact_eigendecomposition
# =================================================================================================

def test_exact_eigendecomposition_not_required_by_default():
    """Test that exact eigendecomposition is not required by default."""
    config = AnalysisConfiguration()
    assert not requires_exact_eigendecomposition(config)


def test_exact_eigendecomposition_required_when_matrices_exact():
    """Test that exact eigendecomposition is required when eigendecomposition_matrices='exact'."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 5
    config.save_eigendecomposition_to_file(filename='exact.npz', operator='exact', form='hamiltonian', shift='unshifted')
    assert requires_exact_eigendecomposition(config)


def test_exact_eigendecomposition_required_when_matrices_both():
    """Test that exact eigendecomposition is required when eigendecomposition_matrices='both'."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 5
    config.save_eigendecomposition_to_file(filename='exact.npz', operator='exact', form='hamiltonian', shift='unshifted')
    config.save_eigendecomposition_to_file(filename='approx.npz', operator='approximate', form='time_evolution', shift='unshifted')
    assert requires_exact_eigendecomposition(config)


def test_exact_eigendecomposition_not_required_when_matrices_approximate():
    """Test that exact eigendecomposition is not required when eigendecomposition_matrices='approximate'."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 5
    config.save_eigendecomposition_to_file(filename='approx.npz', operator='approximate', form='time_evolution', shift='unshifted')
    assert not requires_exact_eigendecomposition(config)


def test_exact_eigendecomposition_required_for_eigenvalue_errors():
    """Test that exact eigendecomposition is required when eigenvalue errors are enabled."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    # Note: In practice, validation would require num_eigenvalues to be set,
    # but the function itself should return True based on enable_eigenvalue_errors alone
    assert requires_exact_eigendecomposition(config)


def test_exact_eigendecomposition_with_all_eigenvalues():
    """Test that exact eigendecomposition works with num_eigenvalues='all'."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 'all'
    config.save_eigendecomposition_to_file(filename='exact.npz', operator='exact', form='hamiltonian', shift='unshifted')
    assert requires_exact_eigendecomposition(config)


# =================================================================================================
# Tests for requires_approximate_eigendecomposition
# =================================================================================================

def test_approximate_eigendecomposition_not_required_by_default():
    """Test that approximate eigendecomposition is not required by default."""
    config = AnalysisConfiguration()
    assert not requires_approximate_eigendecomposition(config)


def test_approximate_eigendecomposition_required_when_matrices_approximate():
    """Test that approximate eigendecomposition is required when eigendecomposition_matrices='approximate'."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 5
    config.save_eigendecomposition_to_file(filename='approx.npz', operator='approximate', form='time_evolution', shift='unshifted')
    assert requires_approximate_eigendecomposition(config)


def test_approximate_eigendecomposition_required_when_matrices_both():
    """Test that approximate eigendecomposition is required when eigendecomposition_matrices='both'."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 5
    config.save_eigendecomposition_to_file(filename='exact.npz', operator='exact', form='hamiltonian', shift='unshifted')
    config.save_eigendecomposition_to_file(filename='approx.npz', operator='approximate', form='time_evolution', shift='unshifted')
    assert requires_approximate_eigendecomposition(config)


def test_approximate_eigendecomposition_not_required_when_matrices_exact():
    """Test that approximate eigendecomposition is not required when eigendecomposition_matrices='exact'."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 5
    config.save_eigendecomposition_to_file(filename='exact.npz', operator='exact', form='hamiltonian', shift='unshifted')
    assert not requires_approximate_eigendecomposition(config)


def test_approximate_eigendecomposition_required_for_eigenvalue_errors():
    """Test that approximate eigendecomposition is required when eigenvalue errors are enabled."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    assert requires_approximate_eigendecomposition(config)


# =================================================================================================
# Tests for requires_exact_matrix
# =================================================================================================

def test_exact_matrix_not_required_by_default():
    """Test that exact matrix is not required by default."""
    config = AnalysisConfiguration()
    assert not requires_exact_matrix(config)


def test_exact_matrix_required_for_output_file():
    """Test that exact matrix is required when exact_matrix_output_file is set."""
    config = AnalysisConfiguration()
    config.save_matrix_to_file(filename='exact.npz', operator='exact', form='hamiltonian', shift='unshifted')
    assert requires_exact_matrix(config)


def test_exact_matrix_required_for_exact_eigendecomposition():
    """Test that exact matrix is required when exact eigendecomposition is needed."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 5
    config.save_eigendecomposition_to_file(filename='exact.npz', operator='exact', form='hamiltonian', shift='unshifted')
    assert requires_exact_matrix(config)


def test_exact_matrix_required_for_both_eigendecomposition():
    """Test that exact matrix is required when both eigendecompositions are needed."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 5
    config.save_eigendecomposition_to_file(filename='exact.npz', operator='exact', form='hamiltonian', shift='unshifted')
    config.save_eigendecomposition_to_file(filename='approx.npz', operator='approximate', form='time_evolution', shift='unshifted')
    assert requires_exact_matrix(config)


def test_exact_matrix_required_for_matrix_norm_errors():
    """Test that exact matrix is required for matrix norm error analysis."""
    config = AnalysisConfiguration()
    config.error_matrix_norms = 'frobenius'
    assert requires_exact_matrix(config)


def test_exact_matrix_required_for_state_errors():
    """Test that exact matrix is required for state-dependent error analysis."""
    config = AnalysisConfiguration()
    config.error_state_inputs = 'state.npy'
    assert requires_exact_matrix(config)


def test_exact_matrix_required_for_eigenvalue_errors():
    """Test that exact matrix is required for eigenvalue error analysis (via eigendecomposition)."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    assert requires_exact_matrix(config)


# =================================================================================================
# Tests for requires_approximate_matrix
# =================================================================================================

def test_approximate_matrix_not_required_by_default():
    """Test that approximate matrix is not required by default."""
    config = AnalysisConfiguration()
    assert not requires_approximate_matrix(config)


def test_approximate_matrix_required_for_output_file():
    """Test that approximate matrix is required when matrix_output_file is set."""
    config = AnalysisConfiguration()
    config.algorithm_matrix_output_file = 'matrix.npz'
    assert requires_approximate_matrix(config)


def test_approximate_matrix_required_for_numerical_simulation():
    """Test that approximate matrix is required for numerical simulation."""
    config = AnalysisConfiguration()
    config.numerical_simulation_inputs = 'state.npy'
    assert requires_approximate_matrix(config)


def test_approximate_matrix_required_for_approximate_eigendecomposition():
    """Test that approximate matrix is required when approximate eigendecomposition is needed."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 5
    config.save_eigendecomposition_to_file(filename='approx.npz', operator='approximate', form='time_evolution', shift='unshifted')
    assert requires_approximate_matrix(config)


def test_approximate_matrix_required_for_both_eigendecomposition():
    """Test that approximate matrix is required when both eigendecompositions are needed."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 5
    config.save_eigendecomposition_to_file(filename='exact.npz', operator='exact', form='hamiltonian', shift='unshifted')
    config.save_eigendecomposition_to_file(filename='approx.npz', operator='approximate', form='time_evolution', shift='unshifted')
    assert requires_approximate_matrix(config)


def test_approximate_matrix_required_for_matrix_norm_errors():
    """Test that approximate matrix is required for matrix norm error analysis."""
    config = AnalysisConfiguration()
    config.error_matrix_norms = 'frobenius'
    assert requires_approximate_matrix(config)


def test_approximate_matrix_required_for_state_errors():
    """Test that approximate matrix is required for state-dependent error analysis."""
    config = AnalysisConfiguration()
    config.error_state_inputs = 'state.npy'
    assert requires_approximate_matrix(config)


def test_approximate_matrix_required_for_eigenvalue_errors():
    """Test that approximate matrix is required for eigenvalue error analysis (via eigendecomposition)."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    assert requires_approximate_matrix(config)


# =================================================================================================
# Integration tests: verify consistency
# =================================================================================================

def test_eigenvalue_errors_requires_both_matrices():
    """Test that eigenvalue errors correctly require both matrices."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True

    assert requires_exact_matrix(config), "Eigenvalue errors should require exact matrix"
    assert requires_approximate_matrix(config), "Eigenvalue errors should require approximate matrix"


def test_matrix_norm_errors_requires_both_matrices():
    """Test that matrix norm errors correctly require both matrices."""
    config = AnalysisConfiguration()
    config.error_matrix_norms = ['frobenius', 'spectral']

    assert requires_exact_matrix(config), "Matrix norm errors should require exact matrix"
    assert requires_approximate_matrix(config), "Matrix norm errors should require approximate matrix"


def test_state_errors_requires_both_matrices():
    """Test that state errors correctly require both matrices."""
    config = AnalysisConfiguration()
    config.error_state_inputs = ['state1.npy', 'state2.npy']

    assert requires_exact_matrix(config), "State errors should require exact matrix"
    assert requires_approximate_matrix(config), "State errors should require approximate matrix"


def test_eigendecomposition_dependency_on_matrices():
    """Test that eigendecomposition correctly implies matrix computation."""
    config = AnalysisConfiguration()
    config.num_eigenvalues = 5
    config.save_eigendecomposition_to_file(filename='exact.npz', operator='exact', form='hamiltonian', shift='unshifted')
    config.save_eigendecomposition_to_file(filename='approx.npz', operator='approximate', form='time_evolution', shift='unshifted')

    # Both eigendecompositions should be required
    assert requires_exact_eigendecomposition(config)
    assert requires_approximate_eigendecomposition(config)

    # And both matrices should be required (because eigendecomposition needs them)
    assert requires_exact_matrix(config)
    assert requires_approximate_matrix(config)
