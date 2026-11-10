"""
Comprehensive test suite for PauliStringEvolution Bloq

Tests both resource estimation pathway (T-complexity) and unitary generation
pathway (tensor_contract).
"""

import numpy as np
import pytest
from scipy.linalg import expm

from qualtran.cirq_interop.t_complexity_protocol import t_complexity

from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.pauli_utils import analytical_evolution, get_pauli_matrix, pauli_string_to_matrix


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
        # Check that more Paulis don't result in simpler algorithms
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


# ==================================================================================
# Test: Exact decomposition and identity positions
# ==================================================================================

# Strings with identities in front, in the middle, at the end, all identities, and 1 qubit.
DECOMPOSITION_STRINGS = ["Z", "X", "Y", "IZ", "ZI", "IIIZ", "XIYIZ", "II", "XYZ", "ZIIX", "IYI"]


class TestExactDecomposition:
    """The decomposition equals the exact matrix, global phase included (not "up to a phase")."""

    @pytest.mark.parametrize("pauli_string", DECOMPOSITION_STRINGS)
    @pytest.mark.parametrize("coefficient, time, hbar", [(0.37, 0.91, 1.0), (-0.62, 1.37, 1.7)])
    def test_decomposition_equals_matrix(self, pauli_string, coefficient, time, hbar):
        pse = PauliStringEvolution(pauli_string, coefficient, time, hbar)
        np.testing.assert_allclose(
            pse.decompose_bloq().tensor_contract(), pse.tensor_contract(), atol=1e-10)

    @pytest.mark.parametrize("pauli_string", DECOMPOSITION_STRINGS)
    def test_decomposition_equals_analytical_matrix(self, pauli_string):
        """Independent of the class's own matrix code."""
        pse = PauliStringEvolution(pauli_string, 0.43, 0.77)
        expected = expm(-1j * 0.43 * 0.77 * pauli_string_to_matrix(pauli_string))
        np.testing.assert_allclose(pse.decompose_bloq().tensor_contract(), expected, atol=1e-10)

    def test_identity_positions_are_untouched(self):
        """Only the two active qubits of "IXIZI" are used: one CNOT pair, no gates on the rest.

        (Cirq's PauliStringPhasorGate treated I as Z, which added 2 CNOTs per identity.)
        """
        cbloq = PauliStringEvolution("IXIZI", 0.5, 0.3).decompose_bloq()
        names = sorted(type(binst.bloq).__name__ for binst in cbloq.bloq_instances)
        assert names == ["CNOT", "CNOT", "Hadamard", "Hadamard", "Rz"]


class TestTrotterContractionsAgree:
    """Trotter's fast contraction and its gate-by-gate "qualtran" contraction still agree."""

    @pytest.mark.parametrize("method", ["first order", "second order", "fourth order"])
    def test_fast_equals_qualtran(self, method):
        from qhat.common.trotter_flattened import Trotterization
        terms = [("XIZ", 0.5), ("IYI", 0.3), ("ZZI", -0.4), ("IXY", 0.2)]
        kwargs = dict(pauli_terms=terms, method=method, time=0.83, num_steps=3)
        fast = Trotterization.from_method(**kwargs).tensor_contract()
        gates = Trotterization.from_method(
            tensor_contraction_method="qualtran", **kwargs).tensor_contract()
        np.testing.assert_allclose(fast, gates, atol=1e-10)


class TestDecompositionCounts:
    """T-complexity equals a hand count of the decomposition."""

    # (string, Cliffords, rotations): H = 1 Clifford, S^dagger/S = 1 each, CNOT = 1.
    #   X: H, H = 2.  Y: S^dagger, H, H, S = 4.  Parity ladder: 2 CNOTs per extra active qubit.
    @pytest.mark.parametrize("pauli_string, cliffords, rotations", [
        ("Z", 0, 1),
        ("IIIZ", 0, 1),
        ("X", 2, 1),
        ("IYI", 4, 1),
        ("IZ", 0, 1),
        ("XIYIZ", 2 + 4 + 0 + 4, 1),      # 3 active qubits -> 4 CNOTs
        ("XYZ", 2 + 4 + 0 + 4, 1),
        ("ZZII", 2, 1),
        ("II", 0, 0),
        ("YY", 8 + 2, 1),
    ])
    def test_hand_count(self, pauli_string, cliffords, rotations):
        tc = t_complexity(PauliStringEvolution(pauli_string, 0.37, 0.91))
        assert (tc.t, tc.clifford, tc.rotations) == (0, cliffords, rotations)

    @pytest.mark.parametrize("pauli_string", DECOMPOSITION_STRINGS)
    def test_matches_counting_the_decomposition(self, pauli_string):
        """The hand-coded `_t_complexity_` agrees with Qualtran counting the decomposition."""
        pse = PauliStringEvolution(pauli_string, 0.37, 0.91)
        assert pse._t_complexity_() == t_complexity(pse.decompose_bloq())

    def test_identities_add_nothing(self):
        """Padding a string with identities does not change its cost."""
        base = t_complexity(PauliStringEvolution("XZ", 0.4, 0.8))
        padded = t_complexity(PauliStringEvolution("IXIIZI", 0.4, 0.8))
        assert base == padded


# ==================================================================================
# Test: adjoint flag and powers
# ==================================================================================

class TestAdjointAndPower:

    @pytest.mark.parametrize("pauli_string", DECOMPOSITION_STRINGS)
    def test_adjoint_is_conjugate_transpose(self, pauli_string):
        pse = PauliStringEvolution(pauli_string, 0.37, 0.91)
        expected = pse.tensor_contract().conj().T
        adj = pse.adjoint()
        assert isinstance(adj, PauliStringEvolution) and adj.is_adjoint
        np.testing.assert_allclose(adj.tensor_contract(), expected, atol=1e-10)
        np.testing.assert_allclose(adj.decompose_bloq().tensor_contract(), expected, atol=1e-10)

    def test_adjoint_keeps_time_and_coefficient(self):
        pse = PauliStringEvolution("XY", 0.37, 0.91)
        adj = pse.adjoint()
        assert (adj.time, adj.coefficient, adj.hbar) == (0.91, 0.37, 1.0)

    def test_double_adjoint_equals_original_with_same_hash(self):
        pse = PauliStringEvolution("XYZ", 0.37, 0.91)
        assert pse.adjoint() != pse
        assert pse.adjoint().adjoint() == pse
        assert hash(pse.adjoint().adjoint()) == hash(pse)

    def test_adjoint_and_negative_time_are_distinct_but_equal_matrices(self):
        pse = PauliStringEvolution("XY", 0.37, 0.91)
        negated = PauliStringEvolution("XY", 0.37, -0.91)
        assert pse.adjoint() != negated
        np.testing.assert_allclose(
            pse.adjoint().tensor_contract(), negated.tensor_contract(), atol=1e-12)

    def test_t_complexity_unchanged(self):
        pse = PauliStringEvolution("XIYZ", 0.37, 0.91)
        assert t_complexity(pse.adjoint()) == t_complexity(pse)

    def test_str_and_repr(self):
        pse = PauliStringEvolution("XY", 0.37, 0.91)
        assert "†" in str(pse.adjoint()) and "†" not in str(pse)
        assert "is_adjoint=True" in repr(pse.adjoint())

    @pytest.mark.parametrize("k", [1, 2, 3, 7])
    def test_power_is_matrix_power(self, k):
        pse = PauliStringEvolution("XYZ", 0.37, 0.91)
        np.testing.assert_allclose(
            (pse ** k).tensor_contract(),
            np.linalg.matrix_power(pse.tensor_contract(), k), atol=1e-10)
        assert (pse ** k).time == pytest.approx(0.91 * k)

    @pytest.mark.parametrize("k", [1, 2, 5])
    def test_negative_power_is_adjoint_of_positive(self, k):
        pse = PauliStringEvolution("XIZ", 0.37, 0.91)
        assert pse ** (-k) == (pse ** k).adjoint()
        np.testing.assert_allclose(
            (pse ** -k).decompose_bloq().tensor_contract(),
            np.linalg.matrix_power(pse.tensor_contract(), k).conj().T, atol=1e-10)

    def test_power_of_adjoint_keeps_flag(self):
        pse = PauliStringEvolution("XZ", 0.37, 0.91)
        assert pse.adjoint() ** 3 == (pse ** 3).adjoint()
        assert pse.adjoint() ** -3 == pse ** 3

    @pytest.mark.parametrize("power", [0, 1.5, 2.0, True])
    def test_invalid_power_raises(self, power):
        with pytest.raises(ValueError):
            PauliStringEvolution("XZ", 0.37, 0.91) ** power


# ==================================================================================
# Test: Rotation counting at special angles (plan D14)
# ==================================================================================

class TestRotationCountIgnoresAngle:
    """Every non-identity term counts as one rotation, whatever its angle."""

    @pytest.mark.parametrize("coefficient, time", [(0.0, 0.91), (np.pi / 8, 1.0), (np.pi / 2, 1.0)])
    def test_special_angles_count_as_one_rotation(self, coefficient, time):
        # Angles 0, pi/4 and pi: Qualtran would count Rz at these as a Clifford or a T.
        tc = t_complexity(PauliStringEvolution("XZ", coefficient, time))
        assert (tc.t, tc.clifford, tc.rotations) == (0, 4, 1)

    def test_fast_trotter_count_is_sum_of_term_counts(self):
        """Includes a zero-coefficient term; R1 finding 2 (13 vs 7 rotations before D14)."""
        from qualtran.cirq_interop.t_complexity_protocol import TComplexity
        from qhat.common.trotter_flattened import Trotterization
        tr = Trotterization((('XZ', .4), ('ZI', 0.), ('YX', .3)), (.5, .5), .7, 3)
        cbloq = tr.decompose_bloq().flatten(
            lambda binst: not isinstance(binst.bloq, PauliStringEvolution))
        terms = [binst.bloq for binst in cbloq.bloq_instances]
        assert all(isinstance(term, PauliStringEvolution) for term in terms)
        total = sum((t_complexity(term) for term in terms), start=TComplexity())
        assert t_complexity(tr) == total
        assert total.rotations == tr.num_operations


# ==================================================================================
# Test: Controlled evolution
# ==================================================================================

def _controlled_matrix(unitary, ctrl_value):
    """U on the block where the control equals `ctrl_value`, the identity elsewhere."""
    dim = unitary.shape[0]
    result = np.eye(2 * dim, dtype=complex)
    result[ctrl_value * dim:(ctrl_value + 1) * dim, ctrl_value * dim:(ctrl_value + 1) * dim] = unitary
    return result


class TestControlled:
    """Controlled versions are our own bloq, exact, and cheap."""

    @pytest.mark.parametrize("pauli_string", DECOMPOSITION_STRINGS)
    @pytest.mark.parametrize("ctrl_value", [1, 0])
    @pytest.mark.parametrize("is_adjoint", [False, True])
    def test_decomposition_equals_block_matrix(self, pauli_string, ctrl_value, is_adjoint):
        """Exact equality (no "up to a phase"), against a matrix built independently."""
        from qualtran import CtrlSpec
        pse = PauliStringEvolution(pauli_string, -0.62, 1.37, 1.7, is_adjoint=is_adjoint)
        forward = expm(-1j * -0.62 * 1.37 / 1.7 * pauli_string_to_matrix(pauli_string))
        expected = _controlled_matrix(forward.conj().T if is_adjoint else forward, ctrl_value)
        cpse = pse.controlled(CtrlSpec(cvs=ctrl_value))
        np.testing.assert_allclose(cpse.decompose_bloq().tensor_contract(), expected, atol=1e-10)
        np.testing.assert_allclose(cpse.tensor_contract(), expected, atol=1e-10)

    def test_returns_our_bloq(self):
        from qualtran import CtrlSpec
        pse = PauliStringEvolution("XZ", 0.37, 0.91)
        assert pse.controlled() == PauliStringEvolution("XZ", 0.37, 0.91, ctrl_values=(1,))
        assert pse.controlled(CtrlSpec(cvs=0)).ctrl_values == (0,)
        assert [reg.name for reg in pse.controlled().signature] == ['ctrl', 'q']
        assert pse.controlled().signature.n_qubits() == 3

    def test_controlled_and_adjoint_commute(self):
        from qualtran import CtrlSpec
        pse = PauliStringEvolution("XYZ", 0.37, 0.91)
        for spec in (CtrlSpec(), CtrlSpec(cvs=0)):
            assert pse.controlled(spec).adjoint() == pse.adjoint().controlled(spec)
            assert type(pse.controlled(spec).adjoint()) is PauliStringEvolution

    def test_power_keeps_controls(self):
        pse = PauliStringEvolution("XZ", 0.37, 0.91).controlled()
        assert (pse ** 3).ctrl_values == (1,) and (pse ** -2).ctrl_values == (1,)
        np.testing.assert_allclose(
            (pse ** -2).decompose_bloq().tensor_contract(),
            np.linalg.matrix_power(pse.tensor_contract(), 2).conj().T, atol=1e-10)

    @pytest.mark.parametrize("ctrl_value", [1, 0])
    def test_only_the_rotation_is_controlled(self, ctrl_value):
        """The basis changes and the ladder aren't controlled; Rz is replaced by 2 Rz + 2 CNOT."""
        from qualtran import CtrlSpec
        cbloq = PauliStringEvolution("IXIZI", 0.5, 0.3).controlled(
            CtrlSpec(cvs=ctrl_value)).decompose_bloq()
        names = sorted(type(binst.bloq).__name__ for binst in cbloq.bloq_instances)
        assert names == ["CNOT"] * 4 + ["Hadamard"] * 2 + ["Rz"] * 2

    # (string, Cliffords, rotations) for the controlled term: the uncontrolled Cliffords + 2 CNOTs,
    # and 2 rotations. All-identity: 1 rotation on the control. Negative controls cost the same.
    @pytest.mark.parametrize("pauli_string, cliffords, rotations", [
        ("Z", 2, 2),
        ("IIIZ", 2, 2),
        ("X", 4, 2),
        ("IYI", 6, 2),
        ("XIYIZ", 2 + 4 + 4 + 2, 2),
        ("YY", 8 + 2 + 2, 2),
        ("II", 0, 1),
    ])
    @pytest.mark.parametrize("ctrl_value", [1, 0])
    def test_hand_count(self, pauli_string, cliffords, rotations, ctrl_value):
        from qualtran import CtrlSpec
        cpse = PauliStringEvolution(pauli_string, 0.37, 0.91).controlled(CtrlSpec(cvs=ctrl_value))
        tc = t_complexity(cpse)
        assert (tc.t, tc.clifford, tc.rotations) == (0, cliffords, rotations)
        assert cpse._t_complexity_() == t_complexity(cpse.decompose_bloq())

    @pytest.mark.parametrize("spec", ["cvs=(1, 1)", "cvs=[1]", "cvs=2", "QUInt(2), cvs=2"])
    def test_unsupported_ctrl_spec_raises(self, spec):
        from qualtran import CtrlSpec, QUInt
        ctrl_spec = eval(f"CtrlSpec({spec})", {"CtrlSpec": CtrlSpec, "QUInt": QUInt})
        with pytest.raises(ValueError, match="single-qubit controls"):
            PauliStringEvolution("XZ", 0.37, 0.91).controlled(ctrl_spec)

    def test_second_control_prepends(self):
        from qualtran import CtrlSpec
        pse = PauliStringEvolution("XZ", 0.37, 0.91)
        assert pse.controlled().controlled() == PauliStringEvolution(
            "XZ", 0.37, 0.91, ctrl_values=(1, 1))
        assert pse.controlled().controlled(CtrlSpec(cvs=0)).ctrl_values == (0, 1)

    @pytest.mark.parametrize("ctrl_values", [(2,), (True,), (-1,)])
    def test_invalid_ctrl_values_raise(self, ctrl_values):
        with pytest.raises(ValueError):
            PauliStringEvolution("XZ", 0.37, 0.91, ctrl_values=ctrl_values)

    def test_str_and_repr(self):
        from qualtran import CtrlSpec
        pse = PauliStringEvolution("XY", 0.37, 0.91).controlled(CtrlSpec(cvs=0))
        assert str(pse).startswith("C[0](")
        assert "ctrl_values=(0,)" in repr(pse)
        assert "ctrl_values" not in repr(PauliStringEvolution("XY", 0.37, 0.91))

    def test_inside_a_composite_bloq(self):
        """The adder wires an outside qubit into `ctrl` (used by FlexibleQPE)."""
        from qualtran import BloqBuilder, CtrlSpec, QBit, Register
        pse = PauliStringEvolution("XZ", 0.37, 0.91)
        bb = BloqBuilder()
        c = bb.add_register('c', 1)
        q = bb.add_register(Register('q', QBit(), shape=(2,)))
        _, adder = pse.get_ctrl_system(CtrlSpec(cvs=0))
        (c,), (q,) = adder(bb, [c], {'q': q})
        cbloq = bb.finalize(c=c, q=q)
        np.testing.assert_allclose(
            cbloq.tensor_contract(), _controlled_matrix(pse.tensor_contract(), 0), atol=1e-10)

    def test_pyliqtr_counts_drop(self):
        """pyLIQTR decomposes the circuit; only the rotation should be controlled (plan D10).

        A controlled term costs pyLIQTR exactly two of the uncontrolled term's rotations, plus 2
        CNOTs. Qualtran's generic `Controlled` wrapper controlled every gate, i.e. 6x the T for
        "XYZI" (1,252 vs 206).
        """
        from qualtran import Controlled, CtrlSpec
        from qhat.analysis.resource_estimation import estimate_resources
        from qhat.common.logging_utils import add_custom_levels
        add_custom_levels()

        def counts(bloq):
            res = estimate_resources("pyliqtr", bloq)["pyliqtr"]
            return res["T_count"], res["Clifford_count"]

        pse = PauliStringEvolution("XYZI", 0.37, 0.91)
        t_u, _ = counts(pse)
        t_new, cliff_new = counts(pse.controlled())
        t_old, cliff_old = counts(Controlled(pse, CtrlSpec()))
        assert t_new == 2 * t_u
        assert t_new < t_old / 5 and cliff_new < cliff_old / 5


class TestCallGraph:
    """`build_call_graph` lists exactly the decomposition's gates, angles included (plan T7)."""

    STRINGS = ["Z", "X", "Y", "I", "IZ", "ZI", "II", "IIIZ", "XIYIZ", "XYZ", "ZIIX", "IYI", "YYXZI"]
    CONTROLS = [((), "ancilla"), ((1,), "ancilla"), ((0,), "ancilla"), ((1, 0), "ancilla"),
                ((0, 1, 1), "ancilla"), ((1, 0), "direct"), ((0, 0), "direct"),
                ((1, 0, 1), "direct")]

    @pytest.mark.parametrize("adjoint", [False, True])
    @pytest.mark.parametrize("ctrls", CONTROLS)
    @pytest.mark.parametrize("pauli_string", STRINGS)
    def test_equals_decomposition_call_graph(self, pauli_string, ctrls, adjoint):
        ctrl_values, method = ctrls
        pse = PauliStringEvolution(pauli_string, 0.37, 0.91, hbar=0.8, is_adjoint=adjoint,
                                   ctrl_values=ctrl_values, multi_control_method=method)
        assert dict(pse.build_call_graph(None)) == dict(pse.decompose_bloq().build_call_graph(None))

    @pytest.mark.parametrize("ctrls", CONTROLS)
    @pytest.mark.parametrize("pauli_string", STRINGS)
    def test_costs_agree_with_t_complexity(self, pauli_string, ctrls):
        """Without the AND ladder (which the two count differently), `QECGatesCost` through
        the call graph gives the same Cliffords and rotations as `t_complexity`."""
        from qualtran.resource_counting import get_cost_value, QECGatesCost
        from qhat.common import controls
        ctrl_values, method = ctrls
        pse = PauliStringEvolution(pauli_string, 0.37, 0.91, ctrl_values=ctrl_values,
                                   multi_control_method=method)
        gates = get_cost_value(pse, QECGatesCost())
        tc = t_complexity(pse)
        assert gates.rotation == tc.rotations
        if not controls.uses_and_ancilla(ctrl_values, method):
            assert (gates.t, gates.clifford, gates.and_bloq) == (0, tc.clifford, 0)

    def test_gate_list(self):
        from qualtran.bloqs.basic_gates import CNOT, Hadamard, Rz
        pse = PauliStringEvolution("IXIZI", 0.37, 0.91)
        assert dict(pse.build_call_graph(None)) == {
            Hadamard(): 2, CNOT(): 2, Rz(2 * 0.37 * 0.91): 1}

    def test_no_zero_counts(self):
        """1-qubit strings have no CNOT chain; it must not appear with count 0."""
        for ctrl_values in [(), (1,)]:
            counts = dict(PauliStringEvolution("Z", 0.3, 0.5, ctrl_values=ctrl_values)
                          .build_call_graph(None))
            assert all(n > 0 for n in counts.values())
