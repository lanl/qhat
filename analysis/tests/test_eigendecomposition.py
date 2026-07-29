"""
Tests for eigendecomposition I/O and operator conversions.

Tests cover:
- Eigendecomposition I/O (save/load)
- Operator conversions via OperatorRepresentation
"""

import numpy as np
import pytest
import tempfile
import os

from qhat.analysis.file_io import save_eigendecomposition, load_eigendecomposition
from qhat.analysis.operators import OperatorRepresentation


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
# Unit Tests: Operator Conversion
# =================================================================================================

def test_convert_unitary_eigenvalues_simple():
    """Test conversion with simple known values using OperatorRepresentation."""
    # For U = exp(-iHt), if E = 1.0 and t = π, then unitary eigenvalue = exp(-iπ) = -1
    tevol_hbar = np.pi
    eigenenergies_true = np.array([1.0])
    unitary_eigenvalues = np.exp(-1j * eigenenergies_true * tevol_hbar)

    # Create OperatorRepresentation from time-evolution eigenvalues
    U_matrix = np.diag(unitary_eigenvalues)
    op = OperatorRepresentation(
        data=U_matrix,
        operator_type='time_evolution',
        energy_shifted=False,
        representation='dense_matrix',
        tevol_hbar=tevol_hbar
    )

    # Get Hamiltonian eigendecomposition
    H_eigen = op.get(
        operator_type='hamiltonian',
        energy_shifted=False,
        representation='eigendecomposition'
    )

    # Should recover the original energies
    np.testing.assert_array_almost_equal(H_eigen['eigenvalues'], eigenenergies_true)
