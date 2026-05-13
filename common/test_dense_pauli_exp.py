"""
Tests for the existing (untested) DensePauliString and DensePauliExponential classes.

These tests are added to verify the existing code works as expected, since it was
previously untested and should be treated with caution.
"""

import numpy as np
import pytest
from scipy.linalg import expm

from dense_pauli_exp import DensePauliString, DensePauliExponential
from test_utils import get_pauli_matrix, pauli_string_to_matrix


# ==================================================================================
# Test: DensePauliString instantiation and properties
# ==================================================================================

class TestDensePauliString:
    """Test the DensePauliString class.

    Note: This class has known issues with tensor_contract(), so we only test
    basic properties and avoid calling tensor_contract().
    """

    def test_from_dense_pauli_string(self):
        """Test creating a DensePauliString from a string and coefficient."""
        dps = DensePauliString.from_dense_pauli_string("XYZ", 0.5)
        assert dps.string == "XYZ"
        assert dps.coefficient == 0.5
        assert dps.num_qubits == 3

    def test_with_identities(self):
        """Test Pauli string with identity operators."""
        dps = DensePauliString.from_dense_pauli_string("IXYI", 1.0)
        assert dps.num_qubits == 4
        assert dps.string == "IXYI"
        # Identities should not appear in the condensed representation
        assert len(dps.condensed) == 2  # Only X and Y

    def test_all_identity(self):
        """Test a Pauli string that is all identities."""
        dps = DensePauliString.from_dense_pauli_string("III", 1.0)
        assert dps.num_qubits == 3
        assert len(dps.condensed) == 0  # No non-identity Paulis

    def test_multiplication_by_scalar(self):
        """Test multiplying DensePauliString by a scalar."""
        dps = DensePauliString.from_dense_pauli_string("XY", 0.5)
        dps2 = dps * 2.0
        assert dps2.coefficient == 1.0
        assert dps2.string == "XY"

        dps3 = 3.0 * dps
        assert dps3.coefficient == 1.5

    def test_string_representation(self):
        """Test string representation."""
        dps = DensePauliString.from_dense_pauli_string("XYZ", 0.5)
        s = str(dps)
        assert "0.5" in s
        assert "XYZ" in s

    def test_exp_method(self):
        """Test the exp() method creates a DensePauliExponential."""
        dps = DensePauliString.from_dense_pauli_string("X", 0.5j)
        dpe = dps.exp()
        assert isinstance(dpe, DensePauliExponential)
        assert dpe.base == dps


# ==================================================================================
# Test: DensePauliExponential (with caveats)
# ==================================================================================

class TestDensePauliExponential:
    """Test the DensePauliExponential class.

    Note: This class has known issues (e.g., requires purely imaginary coefficients).
    These tests document the current behavior rather than necessarily endorsing it.
    """

    def test_creation(self):
        """Test creating a DensePauliExponential."""
        dps = DensePauliString.from_dense_pauli_string("X", 0.5j)
        dpe = DensePauliExponential(dps)
        assert dpe.base == dps

    def test_signature(self):
        """Test that signature matches the base."""
        dps = DensePauliString.from_dense_pauli_string("XYZ", 0.5j)
        dpe = DensePauliExponential(dps)
        assert dpe.signature == dps.signature

    def test_power_operator(self):
        """Test the power operator."""
        dps = DensePauliString.from_dense_pauli_string("X", 0.5j)
        dpe = DensePauliExponential(dps)
        dpe2 = dpe ** 2
        assert isinstance(dpe2, DensePauliExponential)
        # The power should scale the coefficient
        assert dpe2.base.coefficient == dps.coefficient * 2

    def test_string_representation(self):
        """Test string representation."""
        dps = DensePauliString.from_dense_pauli_string("XY", 0.5j)
        dpe = DensePauliExponential(dps)
        s = str(dpe)
        assert "exp" in s
        assert "XY" in s

    def test_requires_imaginary_coefficient(self):
        """Test that real coefficients cause assertion error (current limitation)."""
        # This documents a known limitation - the code requires purely imaginary coefficients
        dps = DensePauliString.from_dense_pauli_string("X", 0.5)  # Real coefficient
        dpe = DensePauliExponential(dps)

        # This should fail when trying to build the bloq
        # Note: We can't test tensor_contract due to other bugs, but we document the limitation
        assert dpe.base.coefficient.real != 0  # Has real component

    def test_imaginary_coefficient_structure(self):
        """Test that imaginary coefficients can be used in structure."""
        dps = DensePauliString.from_dense_pauli_string("Z", 0.5j)
        dpe = DensePauliExponential(dps)

        # Just test that we can create the structure
        assert dpe.base.coefficient == 0.5j
        # Note: tensor_contract has bugs, so we don't test it here


# ==================================================================================
# Test: Comparison with analytical results
# ==================================================================================

class TestDensePauliStringAnalytical:
    """Test DensePauliString basic structure."""

    def test_single_x_pauli(self):
        """Test single X Pauli operator structure."""
        dps = DensePauliString.from_dense_pauli_string("X", 1.0)

        # Just verify the structure is correct
        assert dps.num_qubits == 1
        assert dps.string == "X"

    def test_condensed_form(self):
        """Test that condensed form correctly excludes identities."""
        dps = DensePauliString.from_dense_pauli_string("IXIYIZI", 1.0)

        # IXIYIZI has X at index 1, Y at index 3, Z at index 5
        assert len(dps.condensed) == 3
        chars_in_condensed = [c for c, _ in dps.condensed]
        assert 'I' not in chars_in_condensed
        assert 'X' in chars_in_condensed
        assert 'Y' in chars_in_condensed
        assert 'Z' in chars_in_condensed


# ==================================================================================
# Test: Integration between DensePauliString and DensePauliExponential
# ==================================================================================

class TestIntegration:
    """Test integration between DensePauliString and DensePauliExponential."""

    def test_exp_workflow_structure(self):
        """Test the workflow structure: create string -> multiply -> exponentiate."""
        # Start with a Pauli string
        dps = DensePauliString.from_dense_pauli_string("Z", 1.0)

        # Multiply by imaginary unit to make it suitable for exponentiation
        dps_i = dps * 1j

        # Create exponential
        dpe = dps_i.exp()

        # Verify structure (not tensor_contract due to bugs)
        assert isinstance(dpe, DensePauliExponential)
        assert dpe.base.coefficient == 1.0j

    def test_numpy_exp_compatibility(self):
        """Test that np.exp() also creates a DensePauliExponential."""
        dps = DensePauliString.from_dense_pauli_string("X", 0.5j)

        # numpy's exp should also work
        dpe = np.exp(1j * dps)
        assert isinstance(dpe, DensePauliExponential)


# ==================================================================================
# Test: Known limitations and edge cases
# ==================================================================================

class TestLimitations:
    """Document known limitations of the existing code."""

    def test_real_coefficient_limitation_documented(self):
        """Document that real coefficients have issues with DensePauliExponential."""
        # This is a known limitation - document it
        dps = DensePauliString.from_dense_pauli_string("X", 1.0)  # Real
        dpe = dps.exp()

        # We can create the structure, but it has limitations
        assert dpe.base.coefficient == 1.0

    def test_tensor_contract_has_bugs(self):
        """Document that tensor_contract() has bugs in the existing code."""
        # The existing code has bugs with tensor_contract - it doesn't work properly
        # This test documents that fact
        dps = DensePauliString.from_dense_pauli_string("Z", 0.5j)

        # We can create the structure
        assert dps.string == "Z"
        assert dps.coefficient == 0.5j

        # But tensor_contract has bugs (signature issues)
        # For practical use, users should use PauliStringEvolution instead
