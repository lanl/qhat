"""
Comprehensive test suite for PauliStringEvolution Bloq

Tests both resource estimation pathway (T-complexity) and unitary generation
pathway (tensor_contract).
"""

import numpy as np
import pytest
from scipy.linalg import expm

from qhat.common.pauli_string_evolution import PauliStringEvolution
from qualtran.cirq_interop.t_complexity_protocol import t_complexity
from qhat.common.pauli_utils import get_pauli_matrix, pauli_string_to_matrix, analytical_evolution


# ==================================================================================
# Test: Basic instantiation and validation
# ==================================================================================

class TestInstantiation:
    """Test basic instantiation and input validation."""

    def test_valid_pauli_string(self):
        """Test creation with valid Pauli strings."""
        pse = PauliStringEvolution(pauli_string="XYZ", coefficient=1.0, time=0.5)
        assert pse.pauli_string == "XYZ"
        assert pse.coefficient == 1.0
        assert pse.time == 0.5
        assert pse.hbar == 1.0

    def test_with_identities(self):
        """Test Pauli strings containing identity operators."""
        pse = PauliStringEvolution(pauli_string="IXYI", coefficient=0.5, time=1.0)
        assert pse.pauli_string == "IXYI"
        assert pse.num_qubits == 4

    def test_all_identity(self):
        """Test a Pauli string that is all identities."""
        pse = PauliStringEvolution(pauli_string="III", coefficient=1.0, time=1.0)
        assert pse.num_qubits == 3

    def test_single_qubit(self):
        """Test single-qubit Pauli strings."""
        pse = PauliStringEvolution(pauli_string="X", coefficient=1.0, time=0.5)
        assert pse.num_qubits == 1

    def test_zero_coefficient(self):
        """Test with zero coefficient (should give identity)."""
        pse = PauliStringEvolution(pauli_string="XYZ", coefficient=0.0, time=1.0)
        assert pse.coefficient == 0.0

    def test_negative_time(self):
        """Test with negative time (backward evolution)."""
        pse = PauliStringEvolution(pauli_string="X", coefficient=1.0, time=-0.5)
        assert pse.time == -0.5

    def test_zero_time(self):
        """Test with zero time (should give identity)."""
        pse = PauliStringEvolution(pauli_string="X", coefficient=1.0, time=0.0)
        assert pse.time == 0.0

    def test_custom_hbar(self):
        """Test with non-default hbar."""
        pse = PauliStringEvolution(pauli_string="Z", coefficient=1.0, time=1.0, hbar=2.0)
        assert pse.hbar == 2.0

    def test_invalid_character_raises_error(self):
        """Test that invalid characters raise ValueError."""
        with pytest.raises(ValueError, match="Invalid character"):
            PauliStringEvolution(pauli_string="XAZ", coefficient=1.0, time=1.0)

    def test_empty_string_raises_error(self):
        """Test that empty string raises ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            PauliStringEvolution(pauli_string="", coefficient=1.0, time=1.0)


# ==================================================================================
# Test: Signature properties
# ==================================================================================

class TestSignature:
    """Test the Bloq signature."""

    def test_signature_num_qubits(self):
        """Test that signature matches Pauli string length."""
        pse = PauliStringEvolution(pauli_string="XYZIX", coefficient=1.0, time=1.0)
        assert pse.num_qubits == 5
        sig = pse.signature
        # The signature should have a register 'q' with 5 individual qubits (shape=(5,))
        assert 'q' in [reg.name for reg in sig]
        q_reg = [reg for reg in sig if reg.name == 'q'][0]
        # With shape=(5,), we have 5 individual qubits (bitsize=1 each)
        assert q_reg.shape == (5,)
        assert q_reg.bitsize == 1  # Each element is a single qubit

    def test_single_qubit_signature(self):
        """Test signature for single qubit."""
        pse = PauliStringEvolution(pauli_string="Z", coefficient=1.0, time=1.0)
        assert pse.num_qubits == 1
        sig = pse.signature
        # The signature should have a register 'q' with 1 qubit (shape=(1,))
        assert 'q' in [reg.name for reg in sig]
        q_reg = [reg for reg in sig if reg.name == 'q'][0]
        assert q_reg.shape == (1,)
        assert q_reg.bitsize == 1  # Single qubit


# ==================================================================================
# Test: Unitary matrix generation via tensor_contract
# ==================================================================================

class TestUnitaryGeneration:
    """Test unitary matrix generation via tensor_contract."""

    def test_single_x_gate(self):
        """Test evolution under single X gate."""
        coef = 0.5
        time = 1.0
        pse = PauliStringEvolution(pauli_string="X", coefficient=coef, time=time)
        U = pse.tensor_contract()

        # Analytical result: exp(-i c t X) = cos(ct)I - i sin(ct)X
        ct = coef * time
        expected = np.cos(ct) * np.eye(2) - 1j * np.sin(ct) * get_pauli_matrix('X')

        assert np.allclose(U, expected), f"Mismatch:\nGot:\n{U}\nExpected:\n{expected}"

    def test_single_z_gate(self):
        """Test evolution under single Z gate."""
        coef = 0.7
        time = 0.5
        pse = PauliStringEvolution(pauli_string="Z", coefficient=coef, time=time)
        U = pse.tensor_contract()

        # Analytical result
        ct = coef * time
        expected = np.cos(ct) * np.eye(2) - 1j * np.sin(ct) * get_pauli_matrix('Z')

        assert np.allclose(U, expected)

    def test_single_y_gate(self):
        """Test evolution under single Y gate."""
        coef = 0.3
        time = 2.0
        pse = PauliStringEvolution(pauli_string="Y", coefficient=coef, time=time)
        U = pse.tensor_contract()

        ct = coef * time
        expected = np.cos(ct) * np.eye(2) - 1j * np.sin(ct) * get_pauli_matrix('Y')

        assert np.allclose(U, expected)

    def test_identity_evolution(self):
        """Test evolution under identity (should be identity up to phase)."""
        pse = PauliStringEvolution(pauli_string="I", coefficient=1.0, time=1.0)
        U = pse.tensor_contract()

        # Evolution under identity is just a global phase
        # exp(-i c t I) = exp(-i c t) * I
        expected = np.exp(-1j * 1.0 * 1.0) * np.eye(2)

        assert np.allclose(U, expected)

    def test_two_qubit_iz(self):
        """Test two-qubit IZ evolution."""
        coef = 0.5
        time = 1.0
        pse = PauliStringEvolution(pauli_string="IZ", coefficient=coef, time=time)
        U = pse.tensor_contract()

        expected = analytical_evolution("IZ", coef, time)
        assert np.allclose(U, expected)

    def test_two_qubit_xy(self):
        """Test two-qubit XY evolution."""
        coef = 0.25
        time = 2.0
        pse = PauliStringEvolution(pauli_string="XY", coefficient=coef, time=time)
        U = pse.tensor_contract()

        expected = analytical_evolution("XY", coef, time)
        assert np.allclose(U, expected)

    def test_three_qubit_xyz(self):
        """Test three-qubit XYZ evolution."""
        coef = 0.1
        time = 1.0
        pse = PauliStringEvolution(pauli_string="XYZ", coefficient=coef, time=time)
        U = pse.tensor_contract()

        expected = analytical_evolution("XYZ", coef, time)
        assert np.allclose(U, expected)

    def test_unitarity(self):
        """Test that the generated matrix is unitary: U† U = I."""
        pse = PauliStringEvolution(pauli_string="XYZ", coefficient=0.5, time=1.0)
        U = pse.tensor_contract()

        # Check U† U = I
        U_dag_U = U.conj().T @ U
        identity = np.eye(U.shape[0])

        assert np.allclose(U_dag_U, identity), "Matrix is not unitary"

    def test_coefficient_scaling(self):
        """Test that doubling coefficient is equivalent to doubling time."""
        pauli_str = "XZ"
        time = 1.0

        # Version 1: coefficient = 2, time = t
        pse1 = PauliStringEvolution(pauli_string=pauli_str, coefficient=2.0, time=time)
        U1 = pse1.tensor_contract()

        # Version 2: coefficient = 1, time = 2t
        pse2 = PauliStringEvolution(pauli_string=pauli_str, coefficient=1.0, time=2*time)
        U2 = pse2.tensor_contract()

        assert np.allclose(U1, U2), "Coefficient scaling not equivalent to time scaling"

    def test_composition_property(self):
        """Test that exp(-i H t1) * exp(-i H t2) = exp(-i H (t1+t2))."""
        pauli_str = "YZ"
        coef = 0.3
        t1 = 0.5
        t2 = 0.7

        pse1 = PauliStringEvolution(pauli_string=pauli_str, coefficient=coef, time=t1)
        pse2 = PauliStringEvolution(pauli_string=pauli_str, coefficient=coef, time=t2)
        pse_combined = PauliStringEvolution(pauli_string=pauli_str, coefficient=coef, time=t1+t2)

        U1 = pse1.tensor_contract()
        U2 = pse2.tensor_contract()
        U_combined = pse_combined.tensor_contract()

        U_product = U2 @ U1  # Note: matrix multiplication order

        assert np.allclose(U_product, U_combined), "Composition property violated"

    def test_zero_coefficient_gives_identity(self):
        """Test that zero coefficient gives identity matrix."""
        pse = PauliStringEvolution(pauli_string="XYZ", coefficient=0.0, time=1.0)
        U = pse.tensor_contract()

        expected = np.eye(2**3)
        assert np.allclose(U, expected), "Zero coefficient should give identity"

    def test_zero_time_gives_identity(self):
        """Test that zero time gives identity matrix."""
        pse = PauliStringEvolution(pauli_string="XY", coefficient=0.5, time=0.0)
        U = pse.tensor_contract()

        expected = np.eye(2**2)
        assert np.allclose(U, expected), "Zero time should give identity"

    def test_custom_hbar(self):
        """Test evolution with custom hbar value."""
        pauli_str = "Z"
        coef = 1.0
        time = 1.0
        hbar = 2.0

        pse = PauliStringEvolution(pauli_string=pauli_str, coefficient=coef, time=time, hbar=hbar)
        U = pse.tensor_contract()

        expected = analytical_evolution(pauli_str, coef, time, hbar)
        assert np.allclose(U, expected)

    def test_multiple_of_2pi(self):
        """Test that evolution by 2π returns to identity (up to global phase)."""
        pauli_str = "X"
        coef = 1.0
        time = 2 * np.pi  # Full rotation

        pse = PauliStringEvolution(pauli_string=pauli_str, coefficient=coef, time=time)
        U = pse.tensor_contract()

        # Should be identity up to global phase
        # Remove global phase by dividing by U[0,0]
        U_normalized = U / U[0, 0]
        expected_normalized = np.eye(2)

        assert np.allclose(U_normalized, expected_normalized, atol=1e-10)

    def test_negative_time(self):
        """Test that negative time gives inverse evolution."""
        pauli_str = "XY"
        coef = 0.5
        time = 1.0

        pse_forward = PauliStringEvolution(pauli_string=pauli_str, coefficient=coef, time=time)
        pse_backward = PauliStringEvolution(pauli_string=pauli_str, coefficient=coef, time=-time)

        U_forward = pse_forward.tensor_contract()
        U_backward = pse_backward.tensor_contract()

        # U(-t) should be the inverse (adjoint) of U(t)
        assert np.allclose(U_backward, U_forward.conj().T), "Backward evolution should be inverse"


# ==================================================================================
# Test: Resource estimation via T-complexity
# ==================================================================================

class TestResourceEstimation:
    """Test resource estimation via _t_complexity_."""

    def test_returns_tcomplexity_object(self):
        """Test that _t_complexity_ returns a TComplexity object."""
        pse = PauliStringEvolution(pauli_string="XYZ", coefficient=1.0, time=1.0)
        tc = t_complexity(pse)

        # Should have attributes t, rotations, clifford
        assert hasattr(tc, 't')
        assert hasattr(tc, 'rotations')
        assert hasattr(tc, 'clifford')

    def test_non_negative_counts(self):
        """Test that all resource counts are non-negative."""
        pse = PauliStringEvolution(pauli_string="XY", coefficient=0.5, time=1.0)
        tc = t_complexity(pse)

        assert tc.t >= 0
        assert tc.rotations >= 0
        assert tc.clifford >= 0

    def test_scaling_with_paulis(self):
        """Test that complexity scales with number of non-identity Paulis."""
        # Single Pauli
        pse1 = PauliStringEvolution(pauli_string="X", coefficient=1.0, time=1.0)
        tc1 = t_complexity(pse1)

        # Two Paulis
        pse2 = PauliStringEvolution(pauli_string="XY", coefficient=1.0, time=1.0)
        tc2 = t_complexity(pse2)

        # Three Paulis
        pse3 = PauliStringEvolution(pauli_string="XYZ", coefficient=1.0, time=1.0)
        tc3 = t_complexity(pse3)

        # Complexity should generally increase (though not necessarily strictly)
        # Check that more Paulis don't result in simpler circuits
        # Compare consecutive pairs (transitivity implies 1 vs 3)
        assert tc2.rotations >= tc1.rotations or tc2.t >= tc1.t
        assert tc3.rotations >= tc2.rotations or tc3.t >= tc2.t

    def test_identity_has_minimal_cost(self):
        """Test that evolution under identity has minimal cost."""
        pse_identity = PauliStringEvolution(pauli_string="I", coefficient=1.0, time=1.0)
        tc_identity = t_complexity(pse_identity)

        # Compare against a non-trivial Pauli string
        pse_x = PauliStringEvolution(pauli_string="X", coefficient=1.0, time=1.0)
        tc_x = t_complexity(pse_x)

        # Identity should have cost less than or equal to a non-trivial Pauli
        # Check that all metrics are less than or equal for identity
        assert (tc_identity.t <= tc_x.t and
                tc_identity.rotations <= tc_x.rotations and
                tc_identity.clifford <= tc_x.clifford)


# ==================================================================================
# Test: Integration and compatibility
# ==================================================================================

class TestIntegration:
    """Test integration with Qualtran features."""

    def test_tensor_contract_works(self):
        """Test that tensor_contract (via add_my_tensors) works correctly."""
        pse = PauliStringEvolution(pauli_string="XY", coefficient=0.5, time=1.0)

        # Should be able to get the unitary matrix
        U = pse.tensor_contract()
        assert U.shape == (4, 4)  # 2 qubits = 2^2 x 2^2 matrix

        # Should be unitary
        assert np.allclose(U.conj().T @ U, np.eye(4))

    def test_controlled_operation_structure(self):
        """Test that .controlled() creates a controlled bloq with expected signature."""
        pse = PauliStringEvolution(pauli_string="Z", coefficient=0.5, time=1.0)

        # Create controlled version
        controlled_pse = pse.controlled()

        # Verify it's a valid Bloq
        from qualtran import Bloq
        assert isinstance(controlled_pse, Bloq)

        # Verify it has a signature
        sig = controlled_pse.signature
        assert sig is not None

        # The signature should have more registers than the original
        # (original has 'q', controlled should have 'ctrl' and 'q')
        register_names = [reg.name for reg in sig]
        assert 'ctrl' in register_names  # Control qubit(s)
        assert 'q' in register_names  # Original register

        # Original had 1 qubit (Z acts on 1 qubit)
        # Controlled version should have control + original = at least 2 qubits total
        original_qubits = pse.num_qubits
        assert original_qubits == 1

    def test_str_representation(self):
        """Test string representation."""
        pse = PauliStringEvolution(pauli_string="XYZ", coefficient=0.5, time=1.0, hbar=1.0)
        s = str(pse)
        assert "0.5" in s
        assert "XYZ" in s
        assert "1.0" in s

    def test_repr_representation(self):
        """Test repr representation."""
        pse = PauliStringEvolution(pauli_string="XYZ", coefficient=0.5, time=1.0)
        r = repr(pse)
        assert "PauliStringEvolution" in r
        assert "XYZ" in r


# ==================================================================================
# Test: Numerical precision
# ==================================================================================

class TestNumericalPrecision:
    """Test numerical precision for various parameter ranges."""

    def test_small_time_evolution(self):
        """Test evolution for small time values."""
        pse = PauliStringEvolution(pauli_string="XYZ", coefficient=1.0, time=0.001)
        U = pse.tensor_contract()

        # Should still be unitary
        assert np.allclose(U.conj().T @ U, np.eye(U.shape[0]))

        # Compare with analytical
        expected = analytical_evolution("XYZ", 1.0, 0.001)
        assert np.allclose(U, expected)

    def test_large_time_evolution(self):
        """Test evolution for large time values."""
        pse = PauliStringEvolution(pauli_string="Z", coefficient=1.0, time=100.0)
        U = pse.tensor_contract()

        # Should still be unitary
        assert np.allclose(U.conj().T @ U, np.eye(U.shape[0]))

        expected = analytical_evolution("Z", 1.0, 100.0)
        assert np.allclose(U, expected)

    def test_small_coefficient(self):
        """Test with very small coefficient."""
        pse = PauliStringEvolution(pauli_string="X", coefficient=1e-6, time=1.0)
        U = pse.tensor_contract()

        # Should be very close to identity
        expected = analytical_evolution("X", 1e-6, 1.0)
        assert np.allclose(U, expected)
