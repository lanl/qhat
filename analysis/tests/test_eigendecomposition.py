"""
Tests for eigendecomposition analysis functionality (updated for new API).

Tests cover:
- Eigendecomposition I/O (save/load)
- Full spectrum eigendecomposition (only mode supported)
- Eigenenergy sorting
- Unitary eigenvalue to eigenenergy conversion
- Phase convention [0, 2π)
"""

import numpy as np
import pytest
import tempfile
import os

from qhat.analysis.file_io import save_eigendecomposition, load_eigendecomposition
from qhat.analysis.operators import convert_unitary_eigenvalues_to_eigenenergies


# =================================================================================================
# Unit Tests: Eigendecomposition I/O
# =================================================================================================

def test_save_load_eigendecomposition_exact():
    """Test saving and loading exact eigendecomposition results."""
    # Create test data (sorted eigenenergies)
    eigenenergies = np.array([1.0, 2.0, 3.0])
    eigenvectors = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)

    with tempfile.NamedTemporaryFile(suffix='.npz', delete=False) as f:
        output_file = f.name

    try:
        # Save
        save_eigendecomposition(
            output_file, eigenenergies, eigenvectors,
            matrix_type='exact'
        )

        # Load
        loaded = load_eigendecomposition(output_file)

        # Verify
        np.testing.assert_array_almost_equal(loaded['eigenenergies'], eigenenergies)
        np.testing.assert_array_almost_equal(loaded['eigenvectors'], eigenvectors)
        assert loaded['matrix_type'] == 'exact'
        assert loaded['num_eigenstates'] == 3

    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def test_save_load_eigendecomposition_approximate():
    """Test saving and loading approximate eigendecomposition with timestep."""
    eigenenergies = np.array([0.5, 1.0, 1.5])
    eigenvectors = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    timestep = np.pi
    unitary_eigenvalues = np.exp(-1j * eigenenergies * timestep)
    eigenphases = eigenenergies * timestep

    with tempfile.NamedTemporaryFile(suffix='.npz', delete=False) as f:
        output_file = f.name

    try:
        # Save with optional fields
        save_eigendecomposition(
            output_file, eigenenergies, eigenvectors,
            matrix_type='approximate',
            timestep=timestep,
            unitary_eigenvalues=unitary_eigenvalues,
            eigenphases=eigenphases
        )

        # Load
        loaded = load_eigendecomposition(output_file)

        # Verify
        np.testing.assert_array_almost_equal(loaded['eigenenergies'], eigenenergies)
        np.testing.assert_array_almost_equal(loaded['eigenvectors'], eigenvectors)
        assert loaded['matrix_type'] == 'approximate'
        assert loaded['timestep'] == timestep
        assert 'unitary_eigenvalues' in loaded
        assert 'eigenphases' in loaded

    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def test_save_eigendecomposition_invalid_extension():
    """Test that save_eigendecomposition rejects non-.npz extensions."""
    eigenenergies = np.array([1.0])
    eigenvectors = np.array([[1.0]])

    with pytest.raises(ValueError, match="must end with .npz"):
        save_eigendecomposition(
            "test.txt", eigenenergies, eigenvectors,
            matrix_type='exact'
        )


def test_load_eigendecomposition_invalid_extension():
    """Test that load_eigendecomposition rejects non-.npz extensions."""
    with pytest.raises(ValueError, match="must end with .npz"):
        load_eigendecomposition("test.txt")


# =================================================================================================
# Unit Tests: Phase-to-Energy Conversion
# =================================================================================================

def test_convert_unitary_eigenvalues_simple():
    """Test conversion with simple known values."""
    # For U = exp(-iHt), if E = 1.0 and t = π, then unitary eigenvalue = exp(-iπ) = -1
    timestep = np.pi
    eigenenergies_true = np.array([1.0])
    unitary_eigenvalues = np.exp(-1j * eigenenergies_true * timestep)

    eigenenergies, eigenphases = convert_unitary_eigenvalues_to_eigenenergies(
        unitary_eigenvalues, timestep
    )

    # Should recover the original energies
    np.testing.assert_array_almost_equal(eigenenergies, eigenenergies_true)
    # Phase should be π
    np.testing.assert_array_almost_equal(eigenphases, [np.pi])


def test_convert_unitary_eigenvalues_phase_convention():
    """Test that eigenphases are in [0, 2π) convention."""
    timestep = 1.0
    # Create a range of energies that will produce phases across [0, 2π)
    eigenenergies_true = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    unitary_eigenvalues = np.exp(-1j * eigenenergies_true * timestep)

    eigenenergies, eigenphases = convert_unitary_eigenvalues_to_eigenenergies(
        unitary_eigenvalues, timestep
    )

    # Check phase convention: all phases should be in [0, 2π)
    assert np.all(eigenphases >= 0.0)
    assert np.all(eigenphases < 2 * np.pi)

    # Should recover original energies
    np.testing.assert_array_almost_equal(eigenenergies, eigenenergies_true)


def test_convert_unitary_eigenvalues_multiple():
    """Test conversion with multiple eigenvalues."""
    timestep = np.pi
    eigenenergies_true = np.array([0.5, 1.0, 1.5, 2.0])
    unitary_eigenvalues = np.exp(-1j * eigenenergies_true * timestep)

    eigenenergies, eigenphases = convert_unitary_eigenvalues_to_eigenenergies(
        unitary_eigenvalues, timestep
    )

    # Should recover the original energies
    np.testing.assert_array_almost_equal(eigenenergies, eigenenergies_true)


# =================================================================================================
# Integration Tests: Full Eigendecomposition
# =================================================================================================

def test_eigendecomposition_sorting():
    """Test that eigendecomposition results are sorted by eigenenergy."""
    from qhat.analysis.matrix_eigendecomposition import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    # Create diagonal matrix with unsorted eigenvalues
    unsorted_energies = np.array([3.0, 1.0, 4.0, 2.0])
    matrix = np.diag(unsorted_energies)

    config = AnalysisConfiguration()
    # Use flexible API to request exact eigendecomposition
    config.save_operator_to_file(
        filename='exact_eig.npz',
        source='exact',
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='eigendecomposition'
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(
                config,
                exact_matrix=matrix,
                unitary_matrix=None
            )

            # Check that eigenenergies are sorted
            eigenenergies = results['exact_eigendecomposition']['eigenenergies']
            expected_sorted = np.array([1.0, 2.0, 3.0, 4.0])
            np.testing.assert_array_almost_equal(eigenenergies, expected_sorted)

        finally:
            os.chdir(original_dir)


def test_eigendecomposition_pauli_z():
    """Test eigendecomposition of Pauli Z matrix."""
    from qhat.analysis.matrix_eigendecomposition import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    # Pauli Z = [[1, 0], [0, -1]], eigenvalues: -1, +1
    pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)

    config = AnalysisConfiguration()
    config.save_operator_to_file(
        filename='exact_eig.npz',
        source='exact',
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='eigendecomposition'
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(
                config,
                exact_matrix=pauli_z,
                unitary_matrix=None
            )

            eigenenergies = results['exact_eigendecomposition']['eigenenergies']

            # Should be sorted: -1, +1
            np.testing.assert_array_almost_equal(eigenenergies, [-1.0, 1.0])

        finally:
            os.chdir(original_dir)


def test_eigendecomposition_orthonormality():
    """Test that eigenvectors are orthonormal."""
    from qhat.analysis.matrix_eigendecomposition import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration
    from qhat.analysis.file_io import load_eigendecomposition

    # Create a simple Hermitian matrix
    matrix = np.array([[2, 1], [1, 2]], dtype=complex)

    config = AnalysisConfiguration()
    config.save_operator_to_file(
        filename='exact_eig.npz',
        source='exact',
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='eigendecomposition'
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            results = eigendecomposition_analysis(
                config,
                exact_matrix=matrix,
                unitary_matrix=None
            )

            # Load eigenvectors from the saved file
            output_file = results['exact_eigendecomposition']['file']
            data = load_eigendecomposition(output_file)
            eigenvectors = data['eigenvectors']

            # Check orthonormality: V† V should be identity
            gram = eigenvectors.conj().T @ eigenvectors
            np.testing.assert_array_almost_equal(gram, np.eye(2))

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Error Handling Tests
# =================================================================================================

def test_eigendecomposition_matrix_free_error():
    """Test that matrix-free operators raise appropriate error."""
    from qhat.analysis.matrix_eigendecomposition import _eigendecompose_full
    from qhat.analysis.matrix_operations import PauliStringOperator

    # Create a mock matrix-free operator
    class MockOperator:
        def __init__(self):
            self.shape = (65536, 65536)
        def matvec(self, v):
            return v

    with pytest.raises(ValueError, match="not supported for matrix-free"):
        _eigendecompose_full(MockOperator(), 'exact')


def test_eigendecomposition_no_timestep_error():
    """Test that approximate eigendecomposition without timestep raises error."""
    from qhat.analysis.matrix_eigendecomposition import eigendecomposition_analysis
    from qhat.analysis.config_types import AnalysisConfiguration

    matrix = np.eye(2, dtype=complex)

    config = AnalysisConfiguration()
    # Use flexible API to request approximate eigendecomposition
    config.save_operator_to_file(
        filename='approx_eig.npz',
        source='approximate',
        operator_type='time_evolution',
        energy_shifted=False,
        representation='eigendecomposition'
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            with pytest.raises(ValueError, match="timestep is required"):
                eigendecomposition_analysis(
                    config,
                    timestep=None,  # Missing timestep
                    exact_matrix=None,
                    unitary_matrix=matrix
                )
        finally:
            os.chdir(original_dir)
