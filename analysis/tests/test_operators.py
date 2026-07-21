"""
Tests for the OperatorRepresentation class.

These tests verify:
1. Basic conversions between operator types (H ↔ U)
2. Energy shift application and removal
3. Conversions between representations (dense matrix ↔ eigendecomposition)
4. Caching behavior
5. Round-trip conversions
6. Edge cases and error handling
"""

import numpy as np
import pytest
import scipy.linalg
from qhat.analysis.operators import OperatorRepresentation


class TestBasicCreation:
    """Test basic creation and initialization."""

    def test_create_from_dense_hamiltonian(self):
        """Test creating from a dense Hamiltonian matrix."""
        H = np.array([[1.0, 0.5], [0.5, 2.0]])
        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            energy_shifted=False,
            representation='dense_matrix',
            timestep=1.0
        )
        assert op is not None
        assert op.timestep == 1.0
        assert op.energy_shift == 0.0

    def test_create_from_eigendecomposition(self):
        """Test creating from eigendecomposition."""
        eigenvalues = np.array([1.0, 2.0])
        eigenvectors = np.eye(2)
        eigendata = {'eigenvalues': eigenvalues, 'eigenvectors': eigenvectors}

        op = OperatorRepresentation(
            data=eigendata,
            operator_type='hamiltonian',
            energy_shifted=False,
            representation='eigendecomposition',
            timestep=1.0
        )
        assert op is not None

    def test_invalid_dense_matrix_shape(self):
        """Test that non-square matrices are rejected."""
        H = np.array([[1.0, 0.5, 0.0], [0.5, 2.0, 0.0]])  # Not square
        with pytest.raises(ValueError, match="square"):
            OperatorRepresentation(
                data=H,
                operator_type='hamiltonian',
                representation='dense_matrix',
                timestep=1.0
            )

    def test_invalid_eigendecomposition_format(self):
        """Test that eigendecomposition requires dict with correct keys."""
        eigendata = {'eigenvalues': np.array([1.0, 2.0])}  # Missing eigenvectors
        with pytest.raises(ValueError, match="eigenvalues.*eigenvectors"):
            OperatorRepresentation(
                data=eigendata,
                operator_type='hamiltonian',
                representation='eigendecomposition',
                timestep=1.0
            )


class TestHamiltonianToTimeEvolution:
    """Test conversion from Hamiltonian to time-evolution operator."""

    def test_simple_hamiltonian_to_unitary(self):
        """Test H → U conversion for a simple case."""
        # Simple diagonal Hamiltonian
        H = np.diag([1.0, 2.0])
        timestep = 0.5

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            energy_shifted=False,
            representation='dense_matrix',
            timestep=timestep
        )

        # Get as time-evolution operator
        U = op.get(operator_type='time_evolution', representation='dense_matrix')

        # Expected: U = exp(-i*H*t) = diag(exp(-i*1*0.5), exp(-i*2*0.5))
        expected = np.diag([np.exp(-1j * 1.0 * 0.5), np.exp(-1j * 2.0 * 0.5)])

        np.testing.assert_allclose(U, expected, rtol=1e-14)

    def test_hermitian_hamiltonian_to_unitary(self):
        """Test that Hamiltonian → Unitary produces a unitary matrix."""
        H = np.array([[1.0, 0.5], [0.5, 2.0]])
        timestep = 1.0

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            representation='dense_matrix',
            timestep=timestep
        )

        U = op.get(operator_type='time_evolution', representation='dense_matrix')

        # Check unitarity: U†U = I
        identity = np.eye(2)
        unitarity_error = np.linalg.norm(U.conj().T @ U - identity, 'fro')
        assert unitarity_error < 1e-14

    def test_conversion_requires_timestep(self):
        """Test that H → U conversion fails without timestep."""
        H = np.array([[1.0, 0.5], [0.5, 2.0]])

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            representation='dense_matrix',
            timestep=None  # No timestep provided
        )

        with pytest.raises(ValueError, match="timestep"):
            op.get(operator_type='time_evolution')


class TestTimeEvolutionToHamiltonian:
    """Test conversion from time-evolution operator to Hamiltonian."""

    def test_simple_unitary_to_hamiltonian(self):
        """Test U → H conversion for a simple case."""
        # Create U from known H
        H_original = np.diag([1.0, 2.0])
        timestep = 0.5
        U = scipy.linalg.expm(-1j * H_original * timestep)

        op = OperatorRepresentation(
            data=U,
            operator_type='time_evolution',
            energy_shifted=False,
            representation='dense_matrix',
            timestep=timestep
        )

        # Convert back to H
        H_recovered = op.get(operator_type='hamiltonian', representation='dense_matrix')

        # Should recover original H (up to numerical precision)
        np.testing.assert_allclose(H_recovered, H_original, rtol=1e-12)

    def test_unitary_to_hamiltonian_is_hermitian(self):
        """Test that U → H produces a Hermitian matrix (for unitary U)."""
        # Create a unitary operator
        H = np.array([[1.0, 0.5], [0.5, 2.0]])
        timestep = 1.0
        U = scipy.linalg.expm(-1j * H * timestep)

        op = OperatorRepresentation(
            data=U,
            operator_type='time_evolution',
            representation='dense_matrix',
            timestep=timestep
        )

        H_recovered = op.get(operator_type='hamiltonian', representation='dense_matrix')

        # Check Hermiticity: H = H†
        hermiticity_error = np.linalg.norm(H_recovered - H_recovered.conj().T, 'fro')
        assert hermiticity_error < 1e-12


class TestEnergyShift:
    """Test energy shift application and removal."""

    def test_apply_shift_to_hamiltonian(self):
        """Test applying energy shift to Hamiltonian eigenvalues."""
        H = np.diag([1.0, 2.0])
        energy_shift = 0.5
        timestep = 1.0

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            energy_shifted=False,
            representation='dense_matrix',
            timestep=timestep,
            energy_shift=energy_shift
        )

        # Get shifted Hamiltonian
        H_shifted = op.get(
            operator_type='hamiltonian',
            energy_shifted=True,
            representation='dense_matrix'
        )

        # Expected: H_shifted = H + E*I (NEW CONVENTION, matches driver.py)
        expected = np.diag([1.0 + 0.5, 2.0 + 0.5])
        np.testing.assert_allclose(H_shifted, expected, rtol=1e-14)

    def test_remove_shift_from_hamiltonian(self):
        """Test removing energy shift from Hamiltonian eigenvalues."""
        H_shifted = np.diag([1.5, 2.5])  # This is H + 0.5*I where H = diag([1, 2]) (NEW CONVENTION)
        energy_shift = 0.5
        timestep = 1.0

        op = OperatorRepresentation(
            data=H_shifted,
            operator_type='hamiltonian',
            energy_shifted=True,
            representation='dense_matrix',
            timestep=timestep,
            energy_shift=energy_shift
        )

        # Get unshifted Hamiltonian
        H_unshifted = op.get(
            operator_type='hamiltonian',
            energy_shifted=False,
            representation='dense_matrix'
        )

        # Expected: H = H_shifted - E*I (NEW CONVENTION)
        expected = np.diag([1.0, 2.0])
        np.testing.assert_allclose(H_unshifted, expected, rtol=1e-14)

    def test_apply_shift_to_time_evolution(self):
        """Test applying energy shift to time-evolution operator."""
        # U_unshifted from H = diag([1, 2])
        H = np.diag([1.0, 2.0])
        timestep = 0.5
        energy_shift = 0.3
        U_unshifted = scipy.linalg.expm(-1j * H * timestep)

        op = OperatorRepresentation(
            data=U_unshifted,
            operator_type='time_evolution',
            energy_shifted=False,
            representation='dense_matrix',
            timestep=timestep,
            energy_shift=energy_shift
        )

        # Get shifted time-evolution
        U_shifted = op.get(
            operator_type='time_evolution',
            energy_shifted=True,
            representation='dense_matrix'
        )

        # Expected: U_shifted = exp(-i*E*t)*U_unshifted (NEW CONVENTION, matches driver.py)
        # H_shifted = H + E*I, so U_shifted = exp(-i*(H+E*I)*t) = exp(-i*E*t)*exp(-i*H*t)
        phase_factor = np.exp(-1j * energy_shift * timestep)
        expected = phase_factor * U_unshifted
        np.testing.assert_allclose(U_shifted, expected, rtol=1e-14)


class TestRepresentationConversion:
    """Test conversion between dense matrix and eigendecomposition."""

    def test_dense_to_eigendecomposition(self):
        """Test extracting eigendecomposition from dense matrix."""
        H = np.array([[1.0, 0.5], [0.5, 2.0]])

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            representation='dense_matrix',
            timestep=1.0
        )

        # Get eigendecomposition
        eigendata = op.get(representation='eigendecomposition')

        assert 'eigenvalues' in eigendata
        assert 'eigenvectors' in eigendata

        # Verify it's a valid eigendecomposition
        eigenvalues = eigendata['eigenvalues']
        eigenvectors = eigendata['eigenvectors']

        # Reconstruct H
        H_reconstructed = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.conj().T
        np.testing.assert_allclose(H_reconstructed, H, rtol=1e-12)

    def test_eigendecomposition_to_dense(self):
        """Test reconstructing dense matrix from eigendecomposition."""
        # Create eigendecomposition
        eigenvalues = np.array([0.5, 2.5])
        eigenvectors = np.array([[0.8507, -0.5257], [0.5257, 0.8507]])  # Orthonormal
        eigendata = {'eigenvalues': eigenvalues, 'eigenvectors': eigenvectors}

        op = OperatorRepresentation(
            data=eigendata,
            operator_type='hamiltonian',
            representation='eigendecomposition',
            timestep=1.0
        )

        # Get dense matrix
        H = op.get(representation='dense_matrix')

        # Expected: H = V @ diag(λ) @ V†
        expected = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.conj().T
        np.testing.assert_allclose(H, expected, rtol=1e-14)

    def test_round_trip_dense_eigendecomp_dense(self):
        """Test dense → eigendecomp → dense gives back original."""
        H = np.array([[1.0, 0.5], [0.5, 2.0]])

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            representation='dense_matrix',
            timestep=1.0
        )

        # Convert to eigendecomp and back
        eigendata = op.get(representation='eigendecomposition')
        H_recovered = op.get(representation='dense_matrix')

        np.testing.assert_allclose(H_recovered, H, rtol=1e-12)


class TestComplexConversions:
    """Test complex multi-step conversions."""

    def test_hamiltonian_unshifted_to_time_evolution_shifted(self):
        """Test H (unshifted) → U (shifted) conversion."""
        H = np.diag([1.0, 2.0])
        timestep = 0.5
        energy_shift = 0.3

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            energy_shifted=False,
            representation='dense_matrix',
            timestep=timestep,
            energy_shift=energy_shift
        )

        # Get shifted time-evolution operator
        U_shifted = op.get(
            operator_type='time_evolution',
            energy_shifted=True,
            representation='dense_matrix'
        )

        # Compute expected manually (NEW CONVENTION):
        # 1. H → U: U = exp(-i*H*t)
        # 2. Apply shift: U_shifted = exp(-i*E*t)*U (H_shifted = H + E*I)
        U = scipy.linalg.expm(-1j * H * timestep)
        phase_factor = np.exp(-1j * energy_shift * timestep)  # FIXED: was +1j
        expected = phase_factor * U

        np.testing.assert_allclose(U_shifted, expected, rtol=1e-14)

    def test_round_trip_hamiltonian_to_unitary_and_back(self):
        """Test H → U → H gives back original."""
        H_original = np.array([[1.0, 0.5], [0.5, 2.0]])
        timestep = 1.0

        op = OperatorRepresentation(
            data=H_original,
            operator_type='hamiltonian',
            representation='dense_matrix',
            timestep=timestep
        )

        # Convert H → U → H
        U = op.get(operator_type='time_evolution', representation='dense_matrix')

        # Create new operator from U
        op_U = OperatorRepresentation(
            data=U,
            operator_type='time_evolution',
            representation='dense_matrix',
            timestep=timestep
        )

        H_recovered = op_U.get(operator_type='hamiltonian', representation='dense_matrix')

        # Should match original (up to numerical precision)
        np.testing.assert_allclose(H_recovered, H_original, rtol=1e-10)

    def test_all_four_forms_consistency(self):
        """Test that all 4 forms (H/U × shifted/unshifted) are consistent."""
        H = np.diag([1.0, 2.0])
        timestep = 0.5
        energy_shift = 0.3

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            energy_shifted=False,
            representation='dense_matrix',
            timestep=timestep,
            energy_shift=energy_shift
        )

        # Get all four forms
        H_unshifted = op.get(operator_type='hamiltonian', energy_shifted=False)
        H_shifted = op.get(operator_type='hamiltonian', energy_shifted=True)
        U_unshifted = op.get(operator_type='time_evolution', energy_shifted=False)
        U_shifted = op.get(operator_type='time_evolution', energy_shifted=True)

        # Check relationships (NEW CONVENTION):
        # 1. H_shifted = H_unshifted + E*I (FIXED: was -)
        expected_H_shifted = H_unshifted + energy_shift * np.eye(2)
        np.testing.assert_allclose(H_shifted, expected_H_shifted, rtol=1e-14)

        # 2. U_unshifted = exp(-i*H_unshifted*t)
        expected_U_unshifted = scipy.linalg.expm(-1j * H_unshifted * timestep)
        np.testing.assert_allclose(U_unshifted, expected_U_unshifted, rtol=1e-14)

        # 3. U_shifted = exp(-i*E*t)*U_unshifted (FIXED: was +i)
        phase_factor = np.exp(-1j * energy_shift * timestep)
        expected_U_shifted = phase_factor * U_unshifted
        np.testing.assert_allclose(U_shifted, expected_U_shifted, rtol=1e-14)


class TestCaching:
    """Test that caching works correctly."""

    def test_repeated_get_uses_cache(self):
        """Test that repeated calls don't recompute."""
        H = np.array([[1.0, 0.5], [0.5, 2.0]])

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            representation='dense_matrix',
            timestep=1.0
        )

        # Get U twice
        U1 = op.get(operator_type='time_evolution')
        U2 = op.get(operator_type='time_evolution')

        # Should be the exact same object (cached)
        assert U1 is U2

    def test_different_forms_cached_independently(self):
        """Test that different forms are cached independently."""
        H = np.diag([1.0, 2.0])
        timestep = 1.0
        energy_shift = 0.5

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            representation='dense_matrix',
            timestep=timestep,
            energy_shift=energy_shift
        )

        # Get multiple different forms
        H_unshifted = op.get(operator_type='hamiltonian', energy_shifted=False)
        H_shifted = op.get(operator_type='hamiltonian', energy_shifted=True)
        U_unshifted = op.get(operator_type='time_evolution', energy_shifted=False)
        U_shifted = op.get(operator_type='time_evolution', energy_shifted=True)

        # Get them again - should be same objects
        assert op.get(operator_type='hamiltonian', energy_shifted=False) is H_unshifted
        assert op.get(operator_type='hamiltonian', energy_shifted=True) is H_shifted
        assert op.get(operator_type='time_evolution', energy_shifted=False) is U_unshifted
        assert op.get(operator_type='time_evolution', energy_shifted=True) is U_shifted


class TestEdgeCases:
    """Test edge cases and numerical stability."""

    def test_identity_hamiltonian(self):
        """Test with H = I (identity)."""
        H = np.eye(3)
        timestep = 1.0

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            representation='dense_matrix',
            timestep=timestep
        )

        U = op.get(operator_type='time_evolution')

        # U = exp(-i*I*t) = exp(-i*t)*I
        expected = np.exp(-1j * timestep) * np.eye(3)
        np.testing.assert_allclose(U, expected, rtol=1e-14)

    def test_zero_hamiltonian(self):
        """Test with H = 0 (zero matrix)."""
        H = np.zeros((2, 2))
        timestep = 1.0

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            representation='dense_matrix',
            timestep=timestep
        )

        U = op.get(operator_type='time_evolution')

        # U = exp(-i*0*t) = I
        expected = np.eye(2)
        np.testing.assert_allclose(U, expected, rtol=1e-14)

    def test_large_energy_shift(self):
        """Test with large energy shift values."""
        H = np.diag([1.0, 2.0])
        timestep = 0.5
        energy_shift = 100.0  # Large shift

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            representation='dense_matrix',
            timestep=timestep,
            energy_shift=energy_shift
        )

        # Should still work correctly (NEW CONVENTION)
        H_shifted = op.get(operator_type='hamiltonian', energy_shifted=True)
        expected = np.diag([1.0 + 100.0, 2.0 + 100.0])  # FIXED: was - (sign error)
        np.testing.assert_allclose(H_shifted, expected, rtol=1e-14)

    def test_very_small_timestep(self):
        """Test with very small timestep."""
        H = np.diag([1.0, 2.0])
        timestep = 1e-10  # Very small

        op = OperatorRepresentation(
            data=H,
            operator_type='hamiltonian',
            representation='dense_matrix',
            timestep=timestep
        )

        U = op.get(operator_type='time_evolution')

        # For small t, U ≈ I - i*H*t
        # But we should still get exact result from expm
        expected = scipy.linalg.expm(-1j * H * timestep)
        np.testing.assert_allclose(U, expected, rtol=1e-14)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
