"""
Tests for eigendecomposition analysis functionality.

Tests cover:
- Eigendecomposition of small known matrices
- Full vs partial decomposition modes
- Orthonormality of eigenvectors
- Eigenvalue ordering
- Different matrix types (exact, approximate)
- File I/O for eigendecompositions
"""

import numpy as np
import pytest
import tempfile
import os
from pathlib import Path

from qhat.analysis.file_io import save_eigendecomposition, load_eigendecomposition
from qhat.analysis.matrix_operations import PauliStringOperator


# =================================================================================================
# Unit Tests: Eigendecomposition I/O
# =================================================================================================

def test_save_load_eigendecomposition():
    """Test saving and loading eigendecomposition results."""
    # Create test data
    eigenvalues = np.array([1.0, 2.0, 3.0])
    eigenvectors = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)

    with tempfile.NamedTemporaryFile(suffix='.npz', delete=False) as f:
        output_file = f.name

    try:
        # Save
        save_eigendecomposition(
            output_file, eigenvalues, eigenvectors,
            matrix_type='exact',
            num_eigenvalues=3,
            which_eigenvalues='smallest'
        )

        # Load
        loaded = load_eigendecomposition(output_file)

        # Verify
        np.testing.assert_array_almost_equal(loaded['eigenvalues'], eigenvalues)
        np.testing.assert_array_almost_equal(loaded['eigenvectors'], eigenvectors)
        assert loaded['matrix_type'] == 'exact'
        assert loaded['which_eigenvalues'] == 'smallest'

    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def test_save_eigendecomposition_invalid_extension():
    """Test that save_eigendecomposition rejects non-.npz extensions."""
    eigenvalues = np.array([1.0])
    eigenvectors = np.array([[1.0]])

    with pytest.raises(ValueError, match="must end with .npz"):
        save_eigendecomposition(
            "test.txt", eigenvalues, eigenvectors,
            matrix_type='exact', num_eigenvalues=1, which_eigenvalues='smallest'
        )


def test_load_eigendecomposition_invalid_extension():
    """Test that load_eigendecomposition rejects non-.npz extensions."""
    with pytest.raises(ValueError, match="must end with .npz"):
        load_eigendecomposition("test.txt")


# =================================================================================================
# Unit Tests: Eigendecomposition of Known Matrices
# =================================================================================================

def test_eigendecomposition_identity():
    """Test eigendecomposition of identity matrix."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    # Create 2x2 identity matrix
    identity = np.eye(2, dtype=complex)

    config = AnalysisConfiguration()
    config.num_eigenvalues = "all"
    config.eigendecomposition_matrices = "approximate"

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(config,
                exact_matrix=None,
                unitary_matrix=identity
            )

            # Verify results
            assert 'approximate_eigendecomposition' in results
            assert results['approximate_eigendecomposition']['num_eigenvalues'] == 2

            # Load and check eigenvalues
            loaded = load_eigendecomposition(
                results['approximate_eigendecomposition']['file']
            )
            eigenvalues = loaded['eigenvalues']

            # Identity has all eigenvalues = 1
            np.testing.assert_array_almost_equal(eigenvalues, [1.0, 1.0])

        finally:
            os.chdir(original_dir)


def test_eigendecomposition_pauli_z():
    """Test eigendecomposition of Pauli Z matrix."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    # Pauli Z = [[1, 0], [0, -1]]
    pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)

    config = AnalysisConfiguration()
    config.num_eigenvalues = "all"
    config.eigendecomposition_matrices = "exact"

    # Create a mock hamiltonian object with num_qubits method
    class MockHamiltonian:
        def num_qubits(self):
            return 1

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(config,
                exact_matrix=pauli_z,
                unitary_matrix=None
            )

            # Load and check eigenvalues
            loaded = load_eigendecomposition(
                results['exact_eigendecomposition']['file']
            )
            eigenvalues = loaded['eigenvalues']

            # Pauli Z has eigenvalues -1, 1 (sorted: -1, 1)
            np.testing.assert_array_almost_equal(sorted(eigenvalues), [-1.0, 1.0])

        finally:
            os.chdir(original_dir)


def test_eigendecomposition_orthonormality():
    """Test that eigenvectors are orthonormal."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    # Create a simple Hermitian matrix
    matrix = np.array([[2, 1], [1, 2]], dtype=complex)

    config = AnalysisConfiguration()
    config.num_eigenvalues = "all"
    config.eigendecomposition_matrices = "approximate"

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(config,
                exact_matrix=None,
                unitary_matrix=matrix
            )

            # Load eigenvectors
            loaded = load_eigendecomposition(
                results['approximate_eigendecomposition']['file']
            )
            eigenvectors = loaded['eigenvectors']

            # Check orthonormality: V† V should be identity
            gram = eigenvectors.conj().T @ eigenvectors
            np.testing.assert_array_almost_equal(gram, np.eye(2))

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Unit Tests: Full vs Partial Decomposition
# =================================================================================================

def test_full_eigendecomposition():
    """Test full eigendecomposition mode."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    # 3x3 matrix
    matrix = np.diag([1.0, 2.0, 3.0])

    config = AnalysisConfiguration()
    config.num_eigenvalues = "all"
    config.eigendecomposition_matrices = "approximate"

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(config,
                exact_matrix=None,
                unitary_matrix=matrix
            )

            # Should compute all 3 eigenvalues
            assert results['approximate_eigendecomposition']['num_eigenvalues'] == 3

            loaded = load_eigendecomposition(
                results['approximate_eigendecomposition']['file']
            )
            np.testing.assert_array_almost_equal(loaded['eigenvalues'], [1.0, 2.0, 3.0])

        finally:
            os.chdir(original_dir)


def test_partial_eigendecomposition_smallest():
    """Test partial eigendecomposition with smallest eigenvalues."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    # 5x5 diagonal matrix with known eigenvalues
    matrix = np.diag([5.0, 2.0, 8.0, 1.0, 3.0])

    config = AnalysisConfiguration()
    config.num_eigenvalues = 3
    config.which_eigenvalues = 'smallest'
    config.eigendecomposition_matrices = "approximate"

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(config,
                exact_matrix=None,
                unitary_matrix=matrix
            )

            # Should compute 3 smallest eigenvalues: 1, 2, 3
            loaded = load_eigendecomposition(
                results['approximate_eigendecomposition']['file']
            )
            eigenvalues = sorted(loaded['eigenvalues'])
            np.testing.assert_array_almost_equal(eigenvalues, [1.0, 2.0, 3.0])

        finally:
            os.chdir(original_dir)


def test_partial_eigendecomposition_largest():
    """Test partial eigendecomposition with largest eigenvalues."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    # 5x5 diagonal matrix
    matrix = np.diag([5.0, 2.0, 8.0, 1.0, 3.0])

    config = AnalysisConfiguration()
    config.num_eigenvalues = 3
    config.which_eigenvalues = 'largest'
    config.eigendecomposition_matrices = "approximate"

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(config,
                exact_matrix=None,
                unitary_matrix=matrix
            )

            # Should compute 3 largest eigenvalues: 5, 8, 3
            loaded = load_eigendecomposition(
                results['approximate_eigendecomposition']['file']
            )
            eigenvalues = sorted(loaded['eigenvalues'], reverse=True)
            np.testing.assert_array_almost_equal(eigenvalues, [8.0, 5.0, 3.0])

        finally:
            os.chdir(original_dir)


def test_partial_eigendecomposition_both():
    """Test partial eigendecomposition with both smallest and largest."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    # 7x7 diagonal matrix
    matrix = np.diag([5.0, 2.0, 8.0, 1.0, 3.0, 9.0, 0.5])

    config = AnalysisConfiguration()
    config.num_eigenvalues = 2
    config.which_eigenvalues = 'both'
    config.eigendecomposition_matrices = "approximate"

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(config,
                exact_matrix=None,
                unitary_matrix=matrix
            )

            # Should compute 2 smallest + 2 largest = 4 total
            # Smallest: 0.5, 1.0; Largest: 8.0, 9.0
            assert results['approximate_eigendecomposition']['num_eigenvalues'] == 4

            loaded = load_eigendecomposition(
                results['approximate_eigendecomposition']['file']
            )
            eigenvalues = sorted(loaded['eigenvalues'])
            # Should have 0.5, 1.0, 8.0, 9.0
            expected = sorted([0.5, 1.0, 8.0, 9.0])
            np.testing.assert_array_almost_equal(eigenvalues, expected)

        finally:
            os.chdir(original_dir)


def test_case_insensitive_all():
    """Test that 'all', 'All', 'ALL' all work for full decomposition."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    matrix = np.eye(2, dtype=complex)

    for value in ["all", "All", "ALL"]:
        config = AnalysisConfiguration()
        config.num_eigenvalues = value
        config.eigendecomposition_matrices = "approximate"

        with tempfile.TemporaryDirectory() as tmpdir:
            original_dir = os.getcwd()
            os.chdir(tmpdir)
            try:
                results = eigendecomposition_analysis(config,
                    exact_matrix=None,
                    unitary_matrix=matrix
                )

                # Should compute all eigenvalues
                assert results['approximate_eigendecomposition']['num_eigenvalues'] == 2

            finally:
                os.chdir(original_dir)


# =================================================================================================
# Integration Tests: Matrix Type Selection
# =================================================================================================

def test_eigendecomposition_exact_only():
    """Test eigendecomposition with only exact matrix."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    exact_matrix = np.diag([1.0, 2.0])

    class MockHamiltonian:
        def num_qubits(self):
            return 1

    config = AnalysisConfiguration()
    config.num_eigenvalues = "all"
    config.eigendecomposition_matrices = "exact"

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(config,
                exact_matrix=exact_matrix,
                unitary_matrix=None
            )

            # Should only have exact results
            assert 'exact_eigendecomposition' in results
            assert 'approximate_eigendecomposition' not in results

        finally:
            os.chdir(original_dir)


def test_eigendecomposition_approximate_only():
    """Test eigendecomposition with only approximate matrix."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    approx_matrix = np.diag([1.0, 2.0])

    config = AnalysisConfiguration()
    config.num_eigenvalues = "all"
    config.eigendecomposition_matrices = "approximate"

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(config,
                exact_matrix=None,
                unitary_matrix=approx_matrix
            )

            # Should only have approximate results
            assert 'approximate_eigendecomposition' in results
            assert 'exact_eigendecomposition' not in results

        finally:
            os.chdir(original_dir)


def test_eigendecomposition_both_matrices():
    """Test eigendecomposition with both exact and approximate matrices."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    exact_matrix = np.diag([1.0, 2.0])
    approx_matrix = np.diag([1.1, 1.9])

    class MockHamiltonian:
        def num_qubits(self):
            return 1

    config = AnalysisConfiguration()
    config.num_eigenvalues = "all"
    config.eigendecomposition_matrices = "both"

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(config,
                exact_matrix=exact_matrix,
                unitary_matrix=approx_matrix
            )

            # Should have both results
            assert 'exact_eigendecomposition' in results
            assert 'approximate_eigendecomposition' in results

            # Verify files exist
            assert os.path.exists('exact_eigendecomposition.npz')
            assert os.path.exists('approximate_eigendecomposition.npz')

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Error Handling Tests
# =================================================================================================

def test_full_decomposition_too_large():
    """Test that full decomposition rejects matrix-free operators."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    # Create a mock matrix-free operator
    class MockOperator:
        def __init__(self):
            self.shape = (65536, 65536)
        def matvec(self, v):
            return v

    config = AnalysisConfiguration()
    config.num_eigenvalues = "all"
    config.eigendecomposition_matrices = "exact"

    with pytest.raises(ValueError, match="Full eigendecomposition not supported for matrix-free operators"):
        eigendecomposition_analysis(config,
            exact_matrix=MockOperator(),
            unitary_matrix=None
        )


def test_invalid_num_eigenvalues():
    """Test error handling for invalid num_eigenvalues."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    matrix = np.eye(4, dtype=complex)

    config = AnalysisConfiguration()
    config.num_eigenvalues = 0
    config.eigendecomposition_matrices = "approximate"

    with pytest.raises(ValueError, match="must be positive"):
        eigendecomposition_analysis(config,
            exact_matrix=None,
            unitary_matrix=matrix
        )


def test_num_eigenvalues_exceeds_dimension():
    """Test error when requesting more eigenvalues than dimension."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    matrix = np.eye(4, dtype=complex)

    config = AnalysisConfiguration()
    config.num_eigenvalues = 10  # More than dimension
    config.eigendecomposition_matrices = "approximate"

    with pytest.raises(ValueError, match="must be less than or equal to dimension"):
        eigendecomposition_analysis(config,
            exact_matrix=None,
            unitary_matrix=matrix
        )


def test_invalid_which_eigenvalues():
    """Test error handling for invalid which_eigenvalues."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    matrix = np.eye(4, dtype=complex)

    config = AnalysisConfiguration()
    config.num_eigenvalues = 2
    config.which_eigenvalues = 'invalid'
    config.eigendecomposition_matrices = "approximate"

    with pytest.raises(ValueError, match="must be 'smallest', 'largest', or 'both'"):
        eigendecomposition_analysis(config,
            exact_matrix=None,
            unitary_matrix=matrix
        )


# =================================================================================================
# Integration Test: Verify Full vs Partial Give Same Results
# =================================================================================================

def test_full_vs_partial_consistency():
    """Verify full and partial decomposition give same results for overlap."""
    from qhat.analysis.analysis import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    # Small matrix where we can compute both
    matrix = np.diag([1.0, 2.0, 3.0, 4.0])

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Full decomposition
            config_full = AnalysisConfiguration()
            config_full.num_eigenvalues = "all"
            config_full.eigendecomposition_matrices = "approximate"

            results_full = eigendecomposition_analysis(
                config_full,
                exact_matrix=None,
                unitary_matrix=matrix
            )

            loaded_full = load_eigendecomposition(
                results_full['approximate_eigendecomposition']['file']
            )

            # Partial decomposition (2 smallest)
            os.remove('approximate_eigendecomposition.npz')

            config_partial = AnalysisConfiguration()
            config_partial.num_eigenvalues = 2
            config_partial.which_eigenvalues = 'smallest'
            config_partial.eigendecomposition_matrices = "approximate"

            results_partial = eigendecomposition_analysis(
                config_partial,
                exact_matrix=None,
                unitary_matrix=matrix
            )

            loaded_partial = load_eigendecomposition(
                results_partial['approximate_eigendecomposition']['file']
            )

            # Check: partial should match first 2 of full
            np.testing.assert_array_almost_equal(
                sorted(loaded_partial['eigenvalues']),
                sorted(loaded_full['eigenvalues'])[:2]
            )

        finally:
            os.chdir(original_dir)
