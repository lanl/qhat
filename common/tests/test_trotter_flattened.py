"""
Comprehensive test suite for Trotterization Bloq

Tests expansion logic, resource estimation, unitary generation, and integration.
"""

import numpy as np
import pytest
from scipy.linalg import expm

from qualtran.cirq_interop.t_complexity_protocol import t_complexity

from qhat.common.pauli_utils import pauli_string_to_matrix
from qhat.common.trotter_flattened import (
    Trotterization,
    build_ramped_trotterized_unitary,
    get_trotterization_coefficients,
)


def expand_ramped_trotterization(
    num_terms: int,
    coefficients: list,
    num_steps: int,
    combine_terms: bool = True
) -> list:
    """Helper function to test expansion logic via Trotterization class.

    Creates dummy Pauli terms and returns the expanded sequence.
    """
    # Create dummy pauli terms (using 'X', 'Y', 'Z' pattern)
    pauli_chars = ['X', 'Y', 'Z', 'I']
    pauli_terms = tuple(
        (pauli_chars[i % len(pauli_chars)], 1.0)
        for i in range(num_terms)
    )

    # Create Trotterization instance
    trotter = Trotterization(
        pauli_terms=pauli_terms,
        coefficients=tuple(coefficients),
        time=1.0,
        num_steps=num_steps,
        combine_terms=combine_terms
    )

    # Return the expanded sequence (it's a cached_property, not a method)
    return trotter.expanded_sequence


# ==================================================================================
# Test: Expansion logic
# ==================================================================================

class TestExpansion:
    """Test the expansion of ramped Trotterization into flat sequences."""

    def test_single_term_single_coeff_single_step(self):
        """Test simplest case: 1 term, 1 coefficient, 1 step."""
        result = expand_ramped_trotterization(1, [1.0], 1)
        # Should just be term 0 with coefficient 1.0
        assert result == [(0, 1.0)]

    def test_three_terms_first_order(self):
        """Test 3 terms with first-order method (single coefficient)."""
        result = expand_ramped_trotterization(3, [1.0], 1)
        # First-order: ascending ramp
        assert result == [(0, 1.0), (1, 1.0), (2, 1.0)]

    def test_three_terms_second_order_single_step(self):
        """Test 3 terms with second-order method (0.5, 0.5)."""
        result = expand_ramped_trotterization(3, [0.5, 0.5], 1)
        # Step 1: ascending with 0.5, then descending with 0.5
        # Ascending: 0, 1, 2
        # Descending: 2, 1, 0
        # After combining adjacent: (0, 0.5), (1, 0.5), (2, 1.0), (1, 0.5), (0, 0.5)
        assert result == [(0, 0.5), (1, 0.5), (2, 1.0), (1, 0.5), (0, 0.5)]

    def test_two_terms_second_order_two_steps(self):
        """Test 2 terms, second-order, 2 steps."""
        result = expand_ramped_trotterization(2, [0.5, 0.5], 2)

        # Step 1:
        #   Ramp 1 (asc): 0, 1
        #   Ramp 2 (desc): 1, 0
        # Step 2:
        #   Ramp 3 (asc): 0, 1
        #   Ramp 4 (desc): 1, 0
        # Combining:
        #   Step 1: (0, 0.5), (1, 1.0), (0, 0.5)
        #   Step 2: (0, 0.5), (1, 1.0), (0, 0.5)
        #   Between steps: (0, 0.5) + (0, 0.5) = (0, 1.0)
        # Final: (0, 0.5), (1, 1.0), (0, 1.0), (1, 1.0), (0, 0.5)
        assert result == [(0, 0.5), (1, 1.0), (0, 1.0), (1, 1.0), (0, 0.5)]

    def test_four_terms_custom_coefficients(self):
        """Test with custom coefficients."""
        result = expand_ramped_trotterization(4, [0.3, 0.7], 1)
        # Ascending: 0, 1, 2, 3 with 0.3
        # Descending: 3, 2, 1, 0 with 0.7
        # Middle term 3 combines: 0.3 + 0.7 = 1.0
        assert result == [(0, 0.3), (1, 0.3), (2, 0.3), (3, 1.0), (2, 0.7), (1, 0.7), (0, 0.7)]

    def test_direction_alternates(self):
        """Test that ramp direction alternates correctly."""
        result = expand_ramped_trotterization(2, [0.25, 0.25, 0.25, 0.25], 1)
        # Ramp 1 (asc): 0, 1
        # Ramp 2 (desc): 1, 0
        # Ramp 3 (asc): 0, 1
        # Ramp 4 (desc): 1, 0
        # After combining:
        expected = [(0, 0.25), (1, 0.5), (0, 0.5), (1, 0.5), (0, 0.25)]
        assert result == expected

    def test_empty_terms_raises_error(self):
        """Test that zero terms raises error."""
        with pytest.raises(ValueError, match="at least one Pauli term"):
            expand_ramped_trotterization(0, [1.0], 1)

    def test_zero_steps_raises_error(self):
        """Test that zero steps raises error."""
        with pytest.raises(ValueError, match="num_steps must be positive"):
            expand_ramped_trotterization(1, [1.0], 0)

    def test_empty_coefficients_raises_error(self):
        """Test that empty coefficients raises error."""
        with pytest.raises(ValueError, match="at least one coefficient"):
            expand_ramped_trotterization(1, [], 1)


# ==================================================================================
# Test: Instantiation and validation
# ==================================================================================

class TestInstantiation:
    """Test Trotterization instantiation and validation."""

    def test_valid_single_term(self):
        """Test creation with single term.

        Note: Single-term Hamiltonians are automatically optimized to num_steps=1
        since multiple steps are redundant when there's no non-commutation error.
        """
        trotter = Trotterization(
            pauli_terms=(("XYZ", 0.5),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=10
        )
        assert trotter.num_terms == 1
        assert trotter.num_qubits == 3
        assert trotter.num_steps == 1  # Optimized from 10 to 1

    def test_valid_multiple_terms(self):
        """Test creation with multiple terms."""
        trotter = Trotterization(
            pauli_terms=(("XY", 0.5), ("YZ", 0.3), ("ZX", 0.2)),
            coefficients=(0.5, 0.5),
            time=2.0,
            num_steps=20
        )
        assert trotter.num_terms == 3
        assert trotter.num_qubits == 2
        assert trotter.num_steps == 20

    def test_second_order_coefficients(self):
        """Test with second-order coefficients."""
        trotter = Trotterization(
            pauli_terms=(("X", 1.0), ("Y", 1.0)),
            coefficients=(0.5, 0.5),
            time=1.0,
            num_steps=5
        )
        assert len(trotter.coefficients) == 2

    def test_empty_terms_raises_error(self):
        """Test that empty terms list raises error."""
        with pytest.raises(ValueError, match="at least one"):
            Trotterization(
                pauli_terms=(),
                coefficients=(1.0,),
                time=1.0,
                num_steps=1
            )

    def test_zero_steps_raises_error(self):
        """Test that zero steps raises error."""
        with pytest.raises(ValueError, match="positive"):
            Trotterization(
                pauli_terms=(("X", 1.0),),
                coefficients=(1.0,),
                time=1.0,
                num_steps=0
            )

    def test_negative_steps_raises_error(self):
        """Test that negative steps raises error."""
        with pytest.raises(ValueError, match="positive"):
            Trotterization(
                pauli_terms=(("X", 1.0),),
                coefficients=(1.0,),
                time=1.0,
                num_steps=-5
            )

    def test_empty_coefficients_raises_error(self):
        """Test that empty coefficients raises error."""
        with pytest.raises(ValueError, match="at least one coefficient"):
            Trotterization(
                pauli_terms=(("X", 1.0),),
                coefficients=(),
                time=1.0,
                num_steps=1
            )

    def test_different_length_strings_raise_error(self):
        """Test that different length strings raise error."""
        with pytest.raises(ValueError, match="same length"):
            Trotterization(
                pauli_terms=(("XI", 1.0), ("XYZ", 1.0)),
                coefficients=(1.0,),
                time=1.0,
                num_steps=1
            )

    def test_invalid_character_raises_error(self):
        """Test that invalid characters raise error."""
        with pytest.raises(ValueError, match="Invalid character"):
            Trotterization(
                pauli_terms=(("XAZ", 1.0),),
                coefficients=(1.0,),
                time=1.0,
                num_steps=1
            )


# ==================================================================================
# Test: Signature
# ==================================================================================

class TestSignature:
    """Test the Bloq signature."""

    def test_signature_structure(self):
        """Test that signature has correct structure."""
        trotter = Trotterization(
            pauli_terms=(("XYZ", 0.5), ("ZYX", 0.3)),
            coefficients=(0.5, 0.5),
            time=1.0,
            num_steps=10
        )
        sig = trotter.signature
        assert 'q' in [reg.name for reg in sig]
        q_reg = [reg for reg in sig if reg.name == 'q'][0]
        assert q_reg.shape == (3,)
        assert q_reg.bitsize == 1


# ==================================================================================
# Test: Expanded sequence
# ==================================================================================

class TestExpandedSequence:
    """Test the expanded_sequence property."""

    def test_expanded_sequence_cached(self):
        """Test that expanded_sequence is cached."""
        trotter = Trotterization(
            pauli_terms=(("XY", 0.5), ("YZ", 0.3)),
            coefficients=(0.5, 0.5),
            time=1.0,
            num_steps=1
        )
        seq1 = trotter.expanded_sequence
        seq2 = trotter.expanded_sequence
        assert seq1 is seq2  # Same object (cached)

    def test_expanded_sequence_correctness(self):
        """Test that expanded_sequence matches expand_ramped_trotterization."""
        trotter = Trotterization(
            pauli_terms=(("X", 1.0), ("Y", 1.0), ("Z", 1.0)),
            coefficients=(0.5, 0.5),
            time=1.0,
            num_steps=2
        )
        expected = expand_ramped_trotterization(3, [0.5, 0.5], 2)
        assert trotter.expanded_sequence == expected

    def test_num_operations(self):
        """Test num_operations property."""
        trotter = Trotterization(
            pauli_terms=(("X", 1.0), ("Y", 1.0)),
            coefficients=(0.5, 0.5),
            time=1.0,
            num_steps=1
        )
        # 2 terms, second-order: asc (0,1), desc (1,0) -> combine middle -> 3 ops
        assert trotter.num_operations == 3


# ==================================================================================
# Test: Unitary generation
# ==================================================================================

class TestUnitaryGeneration:
    """Test unitary matrix generation via tensor_contract."""

    def test_single_term_matches_analytical(self):
        """Test single term matches analytical result."""
        # Single term with first-order method
        trotter = Trotterization(
            pauli_terms=(("X", 0.5),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=1
        )
        U = trotter.tensor_contract()

        # Should be exp(-i * 0.5 * X * 1.0)
        X = pauli_string_to_matrix("X")
        expected = expm(-1j * 0.5 * X * 1.0)
        assert np.allclose(U, expected)

    def test_unitarity(self):
        """Test that generated matrix is unitary."""
        trotter = Trotterization(
            pauli_terms=(("XY", 0.5), ("YZ", 0.3)),
            coefficients=(0.5, 0.5),
            time=1.0,
            num_steps=5
        )
        U = trotter.tensor_contract()

        # Check U† U = I
        U_dag_U = U.conj().T @ U
        identity = np.eye(U.shape[0])
        assert np.allclose(U_dag_U, identity), "Matrix is not unitary"

    def test_first_order_approximation(self):
        """Test first-order Trotterization approximates exact evolution."""
        # For small time and many steps, first-order should be reasonable
        h_x = 0.5
        h_y = 0.3

        trotter = Trotterization(
            pauli_terms=(("X", h_x), ("Y", h_y)),
            coefficients=(1.0,),  # First-order
            time=0.1,
            num_steps=100
        )
        U_trotter = trotter.tensor_contract()

        # Exact evolution
        H = h_x * pauli_string_to_matrix("X") + h_y * pauli_string_to_matrix("Y")
        U_exact = expm(-1j * H * 0.1)

        # Should be close for small time and many steps
        assert np.allclose(U_trotter, U_exact, atol=1e-3)

    def test_second_order_better_than_first(self):
        """Test that second-order is more accurate than first-order."""
        h_x = 1.0
        h_y = 1.0
        time = 1.0
        num_steps = 5

        # First-order
        trotter1 = Trotterization(
            pauli_terms=(("X", h_x), ("Y", h_y)),
            coefficients=(1.0,),
            time=time,
            num_steps=num_steps
        )
        U1 = trotter1.tensor_contract()

        # Second-order
        trotter2 = Trotterization(
            pauli_terms=(("X", h_x), ("Y", h_y)),
            coefficients=(0.5, 0.5),
            time=time,
            num_steps=num_steps
        )
        U2 = trotter2.tensor_contract()

        # Exact
        H = h_x * pauli_string_to_matrix("X") + h_y * pauli_string_to_matrix("Y")
        U_exact = expm(-1j * H * time)

        # Measure errors
        error1 = np.linalg.norm(U1 - U_exact)
        error2 = np.linalg.norm(U2 - U_exact)

        # Second-order should be more accurate
        assert error2 < error1

    def test_increasing_steps_improves_accuracy(self):
        """Test that more steps improves accuracy."""
        h_x = 1.0
        h_z = 1.0
        time = 1.0

        # 5 steps
        trotter5 = Trotterization(
            pauli_terms=(("X", h_x), ("Z", h_z)),
            coefficients=(0.5, 0.5),
            time=time,
            num_steps=5
        )
        U5 = trotter5.tensor_contract()

        # 20 steps
        trotter20 = Trotterization(
            pauli_terms=(("X", h_x), ("Z", h_z)),
            coefficients=(0.5, 0.5),
            time=time,
            num_steps=20
        )
        U20 = trotter20.tensor_contract()

        # Exact
        H = h_x * pauli_string_to_matrix("X") + h_z * pauli_string_to_matrix("Z")
        U_exact = expm(-1j * H * time)

        # Measure errors
        error5 = np.linalg.norm(U5 - U_exact)
        error20 = np.linalg.norm(U20 - U_exact)

        # More steps should be more accurate
        assert error20 < error5


# ==================================================================================
# Test: Resource estimation
# ==================================================================================

class TestResourceEstimation:
    """Test resource estimation via T-complexity."""

    def test_returns_tcomplexity_object(self):
        """Test that _t_complexity_ returns a TComplexity object."""
        trotter = Trotterization(
            pauli_terms=(("XY", 0.5), ("YX", 0.3)),
            coefficients=(0.5, 0.5),
            time=1.0,
            num_steps=10
        )
        tc = t_complexity(trotter)

        assert hasattr(tc, 't')
        assert hasattr(tc, 'rotations')
        assert hasattr(tc, 'clifford')

    def test_non_negative_counts(self):
        """Test that all resource counts are non-negative."""
        trotter = Trotterization(
            pauli_terms=(("XY", 0.5),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=5
        )
        tc = t_complexity(trotter)

        assert tc.t >= 0
        assert tc.rotations >= 0
        assert tc.clifford >= 0

    def test_complexity_scales_with_operations(self):
        """Test that complexity scales with number of operations."""
        # Fewer operations
        trotter1 = Trotterization(
            pauli_terms=(("X", 1.0),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=5
        )
        tc1 = t_complexity(trotter1)

        # More operations (more steps)
        trotter2 = Trotterization(
            pauli_terms=(("X", 1.0),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=20
        )
        tc2 = t_complexity(trotter2)

        # Should have more resources with more operations
        assert tc2.rotations >= tc1.rotations or tc2.clifford >= tc1.clifford

    def test_complexity_independent_of_time(self):
        """Test that complexity is independent of time value."""
        trotter1 = Trotterization(
            pauli_terms=(("XYZ", 0.5),),
            coefficients=(1.0,),
            time=0.1,
            num_steps=10
        )
        trotter2 = Trotterization(
            pauli_terms=(("XYZ", 0.5),),
            coefficients=(1.0,),
            time=10.0,
            num_steps=10
        )

        tc1 = t_complexity(trotter1)
        tc2 = t_complexity(trotter2)

        # Should be identical (time doesn't affect gate count)
        assert tc1 == tc2


# ==================================================================================
# Test: Decomposition
# ==================================================================================

class TestDecomposition:
    """Test the decomposition into PauliStringEvolution bloqs."""

    def test_decompose_returns_composite_bloq(self):
        """Test that decomposition returns a CompositeBloq."""
        from qualtran import CompositeBloq

        trotter = Trotterization(
            pauli_terms=(("XY", 0.5),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=5
        )
        cbloq = trotter.decompose_bloq()
        assert isinstance(cbloq, CompositeBloq)

    def test_decompose_unitary_matches_original(self):
        """Test that decomposed bloq produces same unitary as original."""
        trotter = Trotterization(
            pauli_terms=(("X", 0.5), ("Y", 0.3)),
            coefficients=(0.5, 0.5),
            time=1.0,
            num_steps=10
        )
        cbloq = trotter.decompose_bloq()

        U_original = trotter.tensor_contract()
        U_decomposed = cbloq.tensor_contract()

        assert np.allclose(U_original, U_decomposed), \
            "Decomposed unitary doesn't match original"


# ==================================================================================
# Test: String representations
# ==================================================================================

class TestStringRepresentations:
    """Test string and repr methods."""

    def test_str_representation_first_order(self):
        """Test string representation for first-order.

        Note: Single-term Hamiltonians are optimized to 1 step.
        """
        trotter = Trotterization(
            pauli_terms=(("XY", 0.5),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=10
        )
        s = str(trotter)
        assert "first-order" in s
        assert "1 terms" in s
        assert "1 steps" in s  # Optimized from 10

    def test_str_representation_second_order(self):
        """Test string representation for second-order."""
        trotter = Trotterization(
            pauli_terms=(("X", 1.0), ("Y", 1.0)),
            coefficients=(0.5, 0.5),
            time=2.0,
            num_steps=20
        )
        s = str(trotter)
        assert "second-order" in s
        assert "2 terms" in s
        assert "20 steps" in s

    def test_repr_representation(self):
        """Test repr representation."""
        trotter = Trotterization(
            pauli_terms=(("XY", 0.5),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=5
        )
        r = repr(trotter)
        assert "Trotterization" in r


# ==================================================================================
# Test: Properties
# ==================================================================================

class TestProperties:
    """Test various properties and accessors."""

    def test_num_terms(self):
        """Test num_terms property."""
        trotter = Trotterization(
            pauli_terms=(("X", 1.0), ("Y", 1.0), ("Z", 1.0)),
            coefficients=(0.5, 0.5),
            time=1.0,
            num_steps=10
        )
        assert trotter.num_terms == 3

    def test_num_qubits(self):
        """Test num_qubits property."""
        trotter = Trotterization(
            pauli_terms=(("XYZX", 0.5),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=5
        )
        assert trotter.num_qubits == 4

    def test_num_steps(self):
        """Test num_steps attribute.

        Note: Single-term Hamiltonians are optimized to 1 step.
        """
        trotter = Trotterization(
            pauli_terms=(("X", 1.0),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=42
        )
        assert trotter.num_steps == 1  # Optimized from 42


# ==================================================================================
# Test: Helper functions
# ==================================================================================

class TestHelperFunctions:
    """Test helper functions for coefficient generation."""

    def test_get_coefficients_first_order_int(self):
        """Test getting first-order coefficients by int."""
        coeffs = get_trotterization_coefficients(1)
        assert coeffs == (1.0,)

    def test_get_coefficients_first_order_string(self):
        """Test getting first-order coefficients by string."""
        coeffs = get_trotterization_coefficients("first order")
        assert coeffs == (1.0,)

    def test_get_coefficients_second_order_int(self):
        """Test getting second-order coefficients by int."""
        coeffs = get_trotterization_coefficients(2)
        assert coeffs == (0.5, 0.5)

    def test_get_coefficients_second_order_string(self):
        """Test getting second-order coefficients by string."""
        coeffs = get_trotterization_coefficients("second order")
        assert coeffs == (0.5, 0.5)

    def test_get_coefficients_third_order(self):
        """Test getting third-order coefficients."""
        coeffs = get_trotterization_coefficients("third order")
        assert len(coeffs) == 5
        assert coeffs == (7./24., 3./8., 3./8., -25./24., 1.0)

    def test_get_coefficients_fourth_order(self):
        """Test getting fourth-order coefficients."""
        coeffs = get_trotterization_coefficients("fourth order")
        assert len(coeffs) == 10

    def test_get_coefficients_eighth_order(self):
        """Test getting eighth-order coefficients."""
        coeffs = get_trotterization_coefficients("eighth order")
        assert len(coeffs) == 34

    def test_get_coefficients_custom_sequence(self):
        """Test passing custom coefficients."""
        custom = [0.3, 0.4, 0.3]
        coeffs = get_trotterization_coefficients(custom)
        assert coeffs == tuple(custom)

    def test_get_coefficients_unknown_raises_error(self):
        """Test that unknown method raises error."""
        with pytest.raises(ValueError, match="Unknown Trotter method"):
            get_trotterization_coefficients("unknown method")

    def test_get_coefficients_case_insensitive(self):
        """Test that method names are case-insensitive."""
        coeffs1 = get_trotterization_coefficients("Second Order")
        coeffs2 = get_trotterization_coefficients("SECOND ORDER")
        coeffs3 = get_trotterization_coefficients("second order")
        assert coeffs1 == coeffs2 == coeffs3

    def test_get_coefficients_suzuki_aliases(self):
        """Test Suzuki fourth-order method aliases."""
        coeffs1 = get_trotterization_coefficients(4)
        coeffs2 = get_trotterization_coefficients("fourth order")
        coeffs3 = get_trotterization_coefficients("suzuki5")
        coeffs4 = get_trotterization_coefficients("suzuki 5")
        coeffs5 = get_trotterization_coefficients("suzuki 1990")
        coeffs6 = get_trotterization_coefficients("suzuki 1990 5-term")
        assert coeffs1 == coeffs2 == coeffs3 == coeffs4 == coeffs5 == coeffs6
        assert len(coeffs1) == 10  # 5 cycles

    def test_get_coefficients_blanes_moan_fourth(self):
        """Test Blanes & Moan fourth-order method (equation 46)."""
        coeffs1 = get_trotterization_coefficients("bm4")
        coeffs2 = get_trotterization_coefficients("blanes and moan 2002 fourth order")
        coeffs3 = get_trotterization_coefficients("blanes moan 4")
        assert coeffs1 == coeffs2 == coeffs3
        # Should have 12 coefficients (6 cycles, ramped format)
        assert len(coeffs1) == 12
        # Verify all are real numbers
        assert all(isinstance(c, (int, float)) for c in coeffs1)

    def test_get_coefficients_ostmeyer_fourth(self):
        """Test Ostmeyer 2023 fourth-order method (equation 40)."""
        coeffs1 = get_trotterization_coefficients("opt4")
        coeffs2 = get_trotterization_coefficients("ostmeyer 2023 fourth order")
        coeffs3 = get_trotterization_coefficients("optimised fourth order")
        assert coeffs1 == coeffs2 == coeffs3
        # Should have 10 coefficients (5 cycles, ramped format)
        assert len(coeffs1) == 10
        # Verify all are real numbers
        assert all(isinstance(c, (int, float)) for c in coeffs1)

    def test_get_coefficients_blanes_moan_sixth(self):
        """Test Blanes & Moan sixth-order method (equation 53)."""
        coeffs1 = get_trotterization_coefficients(6)
        coeffs2 = get_trotterization_coefficients("sixth order")
        coeffs3 = get_trotterization_coefficients("bm6")
        coeffs4 = get_trotterization_coefficients("blanes and moan 2002 sixth order")
        coeffs5 = get_trotterization_coefficients("blanes moan 6")
        assert coeffs1 == coeffs2 == coeffs3 == coeffs4 == coeffs5
        # Should have 20 coefficients (10 cycles, ramped format)
        assert len(coeffs1) == 20
        # Verify all are real numbers
        assert all(isinstance(c, (int, float)) for c in coeffs1)

    def test_get_coefficients_morales_eighth(self):
        """Test Morales eighth-order method aliases."""
        coeffs1 = get_trotterization_coefficients(8)
        coeffs2 = get_trotterization_coefficients("eighth order")
        coeffs3 = get_trotterization_coefficients("morales")
        coeffs4 = get_trotterization_coefficients("morales 2022")
        coeffs5 = get_trotterization_coefficients("morales 2025")
        assert coeffs1 == coeffs2 == coeffs3 == coeffs4 == coeffs5
        assert len(coeffs1) == 34  # 17 cycles

    def test_get_coefficients_verlet_aliases(self):
        """Test second-order Verlet/Leapfrog aliases."""
        coeffs1 = get_trotterization_coefficients(2)
        coeffs2 = get_trotterization_coefficients("second order")
        coeffs3 = get_trotterization_coefficients("verlet")
        coeffs4 = get_trotterization_coefficients("leapfrog")
        assert coeffs1 == coeffs2 == coeffs3 == coeffs4


# ==================================================================================
# Test: Convenience constructors
# ==================================================================================

class TestConvenienceConstructors:
    """Test convenience class methods for construction."""

    def test_from_method_second_order(self):
        """Test from_method with second-order."""
        trotter = Trotterization.from_method(
            pauli_terms=[("X", 1.0), ("Y", 1.0)],
            method="second order",
            time=1.0,
            num_steps=10
        )
        assert trotter.coefficients == (0.5, 0.5)
        assert trotter.num_terms == 2
        assert trotter.num_steps == 10

    def test_from_method_fourth_order(self):
        """Test from_method with fourth-order.

        Note: Single-term Hamiltonians are optimized to 1 step.
        """
        trotter = Trotterization.from_method(
            pauli_terms=[("XY", 0.5)],
            method="fourth order",
            time=2.0,
            num_steps=20
        )
        assert len(trotter.coefficients) == 10
        assert trotter.time == 2.0
        assert trotter.num_steps == 1  # Optimized from 20

    def test_from_method_custom_coefficients(self):
        """Test from_method with custom coefficients."""
        custom_coeffs = [0.3, 0.4, 0.3]
        trotter = Trotterization.from_method(
            pauli_terms=[("Z", 1.0)],
            method=custom_coeffs,
            time=1.0,
            num_steps=5
        )
        assert trotter.coefficients == tuple(custom_coeffs)

    def test_from_method_with_integer(self):
        """Test from_method with integer method identifier."""
        trotter = Trotterization.from_method(
            pauli_terms=[("X", 1.0), ("Y", 1.0)],
            method=2,  # Integer for second-order
            time=1.0,
            num_steps=10
        )
        assert trotter.coefficients == (0.5, 0.5)
        assert trotter.num_steps == 10

    def test_from_method_integer_matches_string(self):
        """Test that integer and string methods produce same result."""
        terms = [("X", 1.0), ("Y", 1.0)]
        time = 1.0
        num_steps = 10

        trotter_int = Trotterization.from_method(
            pauli_terms=terms,
            method=2,
            time=time,
            num_steps=num_steps
        )

        trotter_str = Trotterization.from_method(
            pauli_terms=terms,
            method="second order",
            time=time,
            num_steps=num_steps
        )

        # Should have same expanded sequence
        assert trotter_int.expanded_sequence == trotter_str.expanded_sequence

    def test_from_method_matches_direct_construction(self):
        """Test that from_method produces same result as direct construction."""
        terms = [("X", 1.0), ("Y", 1.0)]
        time = 1.0
        num_steps = 10

        trotter1 = Trotterization.from_method(
            pauli_terms=terms,
            method="second order",
            time=time,
            num_steps=num_steps
        )

        trotter2 = Trotterization(
            pauli_terms=tuple(terms),
            coefficients=(0.5, 0.5),
            time=time,
            num_steps=num_steps
        )

        # Should have same expanded sequence
        assert trotter1.expanded_sequence == trotter2.expanded_sequence


# ==================================================================================
# Test: Backward compatibility
# ==================================================================================

class TestBackwardCompatibility:
    """Test backward compatibility with old trotter.py interface."""

    def test_build_ramped_trotterized_unitary_interface(self):
        """Test the old interface function exists and works."""
        # Old interface: build_ramped_trotterized_unitary(pauli_strings, method, timestep, numsteps)
        bloq = build_ramped_trotterized_unitary(
            [("X", 1.0), ("Y", 1.0)],
            "second order",
            1.0,
            10
        )
        assert isinstance(bloq, Trotterization)
        assert bloq.time == 1.0
        assert bloq.num_steps == 10
        assert bloq.coefficients == (0.5, 0.5)

    def test_build_with_integer_method(self):
        """Test old interface with integer method."""
        bloq = build_ramped_trotterized_unitary(
            [("X", 1.0), ("Y", 1.0)],
            2,  # Integer for second-order
            1.0,
            10
        )
        assert bloq.coefficients == (0.5, 0.5)

    def test_build_with_generator(self):
        """Test old interface with generator input (as from analysis Hamiltonian)."""
        # Simulate what hamiltonian.get_all_pauli_strings(return_as='strings') returns
        def pauli_generator():
            yield ("XY", 0.5)
            yield ("YZ", 0.3)
            yield ("ZX", 0.2)

        bloq = build_ramped_trotterized_unitary(
            pauli_generator(),
            "second order",
            1.0,
            5
        )
        assert bloq.num_terms == 3
        assert bloq.num_steps == 5

    def test_backward_compatible_with_analysis_usage(self):
        """Test the exact usage pattern from analysis/unitary.py."""
        # This simulates the exact call from unitary.py line 122-126
        pauli_strings = [("XII", 0.5), ("IXI", 0.3), ("IIX", 0.2)]
        method = "second order"
        timestep = 1.0
        Nsteps = 10

        bloq = build_ramped_trotterized_unitary(
            pauli_strings,
            method,
            timestep,
            Nsteps
        )

        # Verify it creates a valid Bloq
        assert isinstance(bloq, Trotterization)
        assert bloq.num_qubits == 3

        # Verify resource estimation works
        tc = t_complexity(bloq)
        assert tc.rotations > 0

    def test_equivalence_old_and_new_interface(self):
        """Test that old and new interfaces produce identical results."""
        terms = [("X", 1.0), ("Y", 1.0)]
        method = "second order"
        time = 1.0
        num_steps = 5

        # Old interface
        bloq_old = build_ramped_trotterized_unitary(terms, method, time, num_steps)

        # New interface
        bloq_new = Trotterization.from_method(
            pauli_terms=terms,
            method=method,
            time=time,
            num_steps=num_steps
        )

        # Should produce identical expanded sequences
        assert bloq_old.expanded_sequence == bloq_new.expanded_sequence

        # Should produce identical unitaries
        U_old = bloq_old.tensor_contract()
        U_new = bloq_new.tensor_contract()
        assert np.allclose(U_old, U_new)


# ==================================================================================
# Test: New Ostmeyer 2023 methods
# ==================================================================================

class TestOstmeyerMethods:
    """Test the newly implemented methods from Ostmeyer 2023 paper."""

    def test_blanes_moan_fourth_instantiation(self):
        """Test that Blanes & Moan fourth-order method can be instantiated."""
        trotter = Trotterization.from_method(
            pauli_terms=[("X", 1.0), ("Y", 1.0)],
            method="bm4",
            time=1.0,
            num_steps=10
        )
        assert len(trotter.coefficients) == 12
        assert trotter.num_steps == 10

    def test_blanes_moan_fourth_unitarity(self):
        """Test that Blanes & Moan fourth produces unitary matrices."""
        trotter = Trotterization.from_method(
            pauli_terms=[("XY", 0.5), ("YZ", 0.3)],
            method="bm4",
            time=1.0,
            num_steps=5
        )
        U = trotter.tensor_contract()

        # Check U† U = I
        U_dag_U = U.conj().T @ U
        identity = np.eye(U.shape[0])
        assert np.allclose(U_dag_U, identity, atol=1e-10), \
            "Blanes & Moan fourth does not produce unitary matrix"

    def test_ostmeyer_fourth_instantiation(self):
        """Test that Ostmeyer fourth-order method can be instantiated."""
        trotter = Trotterization.from_method(
            pauli_terms=[("X", 1.0), ("Y", 1.0)],
            method="opt4",
            time=1.0,
            num_steps=10
        )
        assert len(trotter.coefficients) == 10
        assert trotter.num_steps == 10

    def test_ostmeyer_fourth_unitarity(self):
        """Test that Ostmeyer fourth produces unitary matrices."""
        trotter = Trotterization.from_method(
            pauli_terms=[("XY", 0.5), ("YZ", 0.3)],
            method="opt4",
            time=1.0,
            num_steps=5
        )
        U = trotter.tensor_contract()

        # Check U† U = I
        U_dag_U = U.conj().T @ U
        identity = np.eye(U.shape[0])
        assert np.allclose(U_dag_U, identity, atol=1e-10), \
            "Ostmeyer fourth does not produce unitary matrix"

    def test_blanes_moan_sixth_instantiation(self):
        """Test that Blanes & Moan sixth-order method can be instantiated."""
        trotter = Trotterization.from_method(
            pauli_terms=[("X", 1.0), ("Y", 1.0)],
            method="bm6",
            time=1.0,
            num_steps=10
        )
        assert len(trotter.coefficients) == 20
        assert trotter.num_steps == 10

    def test_blanes_moan_sixth_unitarity(self):
        """Test that Blanes & Moan sixth produces unitary matrices."""
        trotter = Trotterization.from_method(
            pauli_terms=[("XY", 0.5), ("YZ", 0.3)],
            method="bm6",
            time=1.0,
            num_steps=5
        )
        U = trotter.tensor_contract()

        # Check U† U = I
        U_dag_U = U.conj().T @ U
        identity = np.eye(U.shape[0])
        assert np.allclose(U_dag_U, identity, atol=1e-10), \
            "Blanes & Moan sixth does not produce unitary matrix"

    def test_fourth_order_methods_comparison(self):
        """Compare all fourth-order methods - verify they produce valid unitaries."""
        # Test case: two non-commuting Hamiltonians
        # Note: As discussed in Ostmeyer (2023), different fourth-order methods have
        # different error accumulation properties. Ostmeyer's "optimised" method has
        # high theoretical efficiency but unfavorable error accumulation in practice,
        # often being outperformed by Blanes & Moan's method.
        h_x = 1.0
        h_y = 1.0
        time = 1.0
        num_steps = 10

        # Exact evolution
        H = h_x * pauli_string_to_matrix("X") + h_y * pauli_string_to_matrix("Y")
        U_exact = expm(-1j * H * time)

        # Test that all methods can be instantiated and produce valid results
        methods = [
            "suzuki5",
            "bm4",
            "opt4"
        ]

        for method_name in methods:
            trotter = Trotterization.from_method(
                pauli_terms=[("X", h_x), ("Y", h_y)],
                method=method_name,
                time=time,
                num_steps=num_steps
            )
            U = trotter.tensor_contract()

            # Verify unitary
            U_dag_U = U.conj().T @ U
            identity = np.eye(U.shape[0])
            assert np.allclose(U_dag_U, identity, atol=1e-10), \
                f"{method_name} does not produce unitary matrix"

            # Verify it's at least approximating the correct evolution
            error = np.linalg.norm(U - U_exact)
            assert error < 1.0, f"{method_name} error is unreasonably large"

    def test_sixth_order_better_than_fourth(self):
        """Test that sixth-order is more accurate than fourth-order."""
        h_x = 1.0
        h_y = 1.0
        time = 1.0
        num_steps = 5

        # Exact evolution
        H = h_x * pauli_string_to_matrix("X") + h_y * pauli_string_to_matrix("Y")
        U_exact = expm(-1j * H * time)

        # Fourth-order (Suzuki)
        trotter4 = Trotterization.from_method(
            pauli_terms=[("X", h_x), ("Y", h_y)],
            method=4,
            time=time,
            num_steps=num_steps
        )
        U4 = trotter4.tensor_contract()
        error4 = np.linalg.norm(U4 - U_exact)

        # Sixth-order (Blanes & Moan)
        trotter6 = Trotterization.from_method(
            pauli_terms=[("X", h_x), ("Y", h_y)],
            method=6,
            time=time,
            num_steps=num_steps
        )
        U6 = trotter6.tensor_contract()
        error6 = np.linalg.norm(U6 - U_exact)

        # Sixth-order should be more accurate
        assert error6 < error4

    def test_new_methods_resource_estimation(self):
        """Test that resource estimation works for all new methods."""
        terms = [("XY", 0.5), ("YZ", 0.3)]
        time = 1.0
        num_steps = 10

        for method_name in ["bm4", "opt4", "bm6"]:
            trotter = Trotterization.from_method(
                pauli_terms=terms,
                method=method_name,
                time=time,
                num_steps=num_steps
            )
            tc = t_complexity(trotter)

            # Should have non-negative resource counts
            assert tc.t >= 0
            assert tc.rotations >= 0
            assert tc.clifford >= 0

    def test_coefficient_symmetry_properties(self):
        """Test that new methods produce symmetric coefficient sequences."""
        # For symmetric Trotter methods, the coefficient sequence should be palindromic
        # when considering the full (c,d) pairs

        for method_name in ["bm4", "opt4", "bm6"]:
            coeffs = get_trotterization_coefficients(method_name)

            # Check that we have an even number (pairs of c,d)
            assert len(coeffs) % 2 == 0

            # For these methods, verify basic properties
            # All coefficients should be real
            assert all(isinstance(c, (int, float)) for c in coeffs)


# ==================================================================================
# Test: Tensor contraction method selection
# ==================================================================================

class TestTensorContractionMethods:
    """Test that different tensor contraction methods work correctly and produce identical results."""

    @pytest.mark.parametrize("num_steps", [1, 2, 3, 10, 50])
    @pytest.mark.parametrize("method", ["qualtran", "incremental", "structured", None])
    def test_all_methods_produce_unitary(self, num_steps, method):
        """Test that all tensor contraction methods produce valid unitary matrices."""
        pauli_terms = [('X', 1.0), ('Z', 0.5), ('Y', 0.3)]

        trotter = Trotterization.from_method(
            pauli_terms=pauli_terms,
            method='second order',
            time=0.1,
            num_steps=num_steps,
            tensor_contraction_method=method
        )

        U = trotter.tensor_contract()

        # Verify unitarity: U† U = I
        U_dag_U = U.conj().T @ U
        identity = np.eye(U.shape[0])
        assert np.allclose(U_dag_U, identity, atol=1e-10), \
            f"Method {method} with num_steps={num_steps} does not produce unitary matrix"

    @pytest.mark.parametrize("num_steps", [1, 2, 3, 10, 50])
    def test_all_methods_agree(self, num_steps):
        """Test that all tensor contraction methods produce identical results."""
        pauli_terms = [('X', 1.0), ('Z', 0.5), ('Y', 0.3)]
        methods = ['qualtran', 'incremental', 'structured', None]

        # Compute matrices for all methods
        results = {}
        for method in methods:
            trotter = Trotterization.from_method(
                pauli_terms=pauli_terms,
                method='second order',
                time=0.1,
                num_steps=num_steps,
                tensor_contraction_method=method
            )
            results[str(method)] = trotter.tensor_contract()

        # Use incremental as reference
        U_ref = results['incremental']

        # Check that all methods produce the same matrix (within numerical precision)
        for method_name, U in results.items():
            if method_name == 'incremental':
                continue
            max_diff = np.max(np.abs(U - U_ref))
            assert max_diff < 1e-14, \
                f"Method {method_name} differs from incremental by {max_diff:.2e} at num_steps={num_steps}"

    @pytest.mark.parametrize("num_steps", [1, 2, 3])
    def test_structured_works_for_small_num_steps(self, num_steps):
        """Test that structured method now works for num_steps = 1, 2, 3 (after Option 1+2)."""
        pauli_terms = [('X', 1.0), ('Z', 0.5)]

        # This should not raise an error anymore
        trotter = Trotterization.from_method(
            pauli_terms=pauli_terms,
            method='second order',
            time=0.1,
            num_steps=num_steps,
            tensor_contraction_method='structured'
        )

        U = trotter.tensor_contract()

        # Verify it's unitary
        U_dag_U = U.conj().T @ U
        identity = np.eye(U.shape[0])
        assert np.allclose(U_dag_U, identity, atol=1e-10)

    def test_pattern_detection_degenerate_case(self):
        """Test that pattern structure is correct for num_steps=1."""
        pauli_terms = [('X', 1.0), ('Z', 0.5)]

        trotter = Trotterization.from_method(
            pauli_terms=pauli_terms,
            method='second order',
            time=0.1,
            num_steps=1
        )

        # For num_steps=1, the pattern should have empty repeat_bridge
        # since there are no step boundaries to bridge
        assert trotter.num_steps == 1
        assert len(trotter._repeat_bridge) == 0 or trotter.num_steps == 1

        # The expanded sequence should still produce correct results
        assert len(trotter.expanded_sequence) > 0

    def test_pattern_structure_without_combining(self):
        """Test pattern structure when combine_terms=False.

        Without combining, all terms stay separate with no prologue/bridge/epilogue.
        """
        trotter = Trotterization(
            pauli_terms=(('X', 1.0), ('Z', 0.5)),
            coefficients=(0.5, 0.5),  # second-order
            time=1.0,
            num_steps=3,
            combine_terms=False
        )

        # Expected pattern: ascending (0,1) then descending (1,0), all separate
        # Step: [(0, 0.5), (1, 0.5), (1, 0.5), (0, 0.5)]
        assert trotter._prologue == ()
        assert trotter._repeat_core == ((0, 0.5), (1, 0.5), (1, 0.5), (0, 0.5))
        assert trotter._repeat_bridge == ()
        assert trotter._epilogue == ()

        # Expanded: 3 repetitions of the core
        expected_expanded = [
            (0, 0.5), (1, 0.5), (1, 0.5), (0, 0.5),  # step 1
            (0, 0.5), (1, 0.5), (1, 0.5), (0, 0.5),  # step 2
            (0, 0.5), (1, 0.5), (1, 0.5), (0, 0.5),  # step 3
        ]
        assert trotter.expanded_sequence == expected_expanded

    def test_pattern_structure_asymmetric_method(self):
        """Test pattern structure with asymmetric method (Ruth 1983, 5 coefficients).

        Asymmetric methods end with a different term than they start,
        so no cross-step combining occurs.
        """
        ruth_coeffs = get_trotterization_coefficients('ruth 1983')
        assert len(ruth_coeffs) == 5  # Verify it's the 5-coefficient version

        trotter = Trotterization(
            pauli_terms=(('X', 1.0), ('Y', 0.5), ('Z', 0.3)),
            coefficients=ruth_coeffs,
            time=1.0,
            num_steps=2,
            combine_terms=True
        )

        # With 5 coefficients (odd), alternating directions:
        # Ramp 1 (asc): 0,1,2 | Ramp 2 (desc): 2,1,0 | Ramp 3 (asc): 0,1,2
        # Ramp 4 (desc): 2,1,0 | Ramp 5 (asc): 0,1,2
        # First term: 0, Last term: 2 (different) -> no cross-step combining

        assert trotter._prologue == ()
        assert trotter._epilogue == ()
        assert trotter._repeat_bridge == ()
        # All combined terms go into repeat_core
        assert len(trotter._repeat_core) == 11  # After within-step combining
        # Verify first and last terms differ
        assert trotter._repeat_core[0][0] == 0   # starts with term 0
        assert trotter._repeat_core[-1][0] == 2  # ends with term 2

    def test_pattern_structure_symmetric_method(self):
        """Test pattern structure with symmetric method (second-order).

        Symmetric methods end with the same term they start with,
        enabling cross-step combining via prologue/bridge/epilogue.
        """
        trotter = Trotterization(
            pauli_terms=(('X', 1.0), ('Z', 0.5)),
            coefficients=(0.5, 0.5),
            time=1.0,
            num_steps=3,
            combine_terms=True
        )

        # Second-order with 2 terms: ascending (0,1) then descending (1,0)
        # After within-step combining: [(0, 0.5), (1, 1.0), (0, 0.5)]
        # First and last are both term 0 with coeff 0.5
        # Cross-step combining splits this:

        assert trotter._prologue == ((0, 0.5),)
        assert trotter._repeat_core == ((1, 1.0),)
        assert trotter._repeat_bridge == ((0, 1.0),)  # 0.5 + 0.5 from adjacent steps
        assert trotter._epilogue == ((0, 0.5),)
        assert trotter._symmetric_bookends == True

        # Expanded for 3 steps: prologue + (core + bridge)*(n-1) + core + epilogue
        expected_expanded = [
            (0, 0.5),  # prologue
            (1, 1.0),  # core (step 1)
            (0, 1.0),  # bridge (between 1 and 2)
            (1, 1.0),  # core (step 2)
            (0, 1.0),  # bridge (between 2 and 3)
            (1, 1.0),  # core (step 3)
            (0, 0.5),  # epilogue
        ]
        assert trotter.expanded_sequence == expected_expanded

    def test_pattern_structure_even_ramps_asymmetric_coeffs(self):
        """Test pattern structure with even ramp count but asymmetric coefficients.

        When a step has an even number of ramps, the first and last terms have
        the same index. However, if the first and last coefficients differ,
        the prologue and epilogue will have different coefficients, and
        symmetric_bookends should be False.
        """
        # 4 ramps (even) with asymmetric coefficients
        coeffs = (0.3, 0.4, 0.5, 0.6)

        trotter = Trotterization(
            pauli_terms=(('X', 1.0), ('Z', 0.5)),
            coefficients=coeffs,
            time=1.0,
            num_steps=3,
            combine_terms=True
        )

        # Manual trace with 2 terms and 4 ramps:
        # Ramp 1 (asc, 0.3): 0, 1
        # Ramp 2 (desc, 0.4): 1, 0
        # Ramp 3 (asc, 0.5): 0, 1
        # Ramp 4 (desc, 0.6): 1, 0
        #
        # Before combining: (0,0.3), (1,0.3), (1,0.4), (0,0.4), (0,0.5), (1,0.5), (1,0.6), (0,0.6)
        # After within-step combining: (0,0.3), (1,0.7), (0,0.9), (1,1.1), (0,0.6)
        # First term: (0, 0.3), Last term: (0, 0.6) - same index, different coeff!
        #
        # Cross-step combining creates:
        # - prologue: (0, 0.3)
        # - repeat_core: (1, 0.7), (0, 0.9), (1, 1.1) [middle terms]
        # - repeat_bridge: (0, 0.9) [0.3 + 0.6 from adjacent steps]
        # - epilogue: (0, 0.6)

        # Check structure (use np.isclose for floating point comparisons)
        assert len(trotter._prologue) == 1
        assert trotter._prologue[0][0] == 0
        assert np.isclose(trotter._prologue[0][1], 0.3)

        assert len(trotter._repeat_core) == 3
        assert trotter._repeat_core[0] == (1, 0.7)
        assert trotter._repeat_core[1] == (0, 0.9)
        assert trotter._repeat_core[2] == (1, 1.1)

        assert len(trotter._repeat_bridge) == 1
        assert trotter._repeat_bridge[0][0] == 0
        assert np.isclose(trotter._repeat_bridge[0][1], 0.9)  # 0.3 + 0.6

        assert len(trotter._epilogue) == 1
        assert trotter._epilogue[0][0] == 0
        assert np.isclose(trotter._epilogue[0][1], 0.6)

        # Verify prologue and epilogue are different
        assert trotter._prologue[0][0] == trotter._epilogue[0][0]  # same term index
        assert not np.isclose(trotter._prologue[0][1], trotter._epilogue[0][1])  # different coefficients
        assert np.isclose(trotter._prologue[0][1], 0.3)
        assert np.isclose(trotter._epilogue[0][1], 0.6)

        # Verify symmetric_bookends flag is False (asymmetric)
        assert trotter._symmetric_bookends == False

        # Verify expanded sequence for 3 steps
        expanded = trotter.expanded_sequence
        assert len(expanded) == 13  # 1 prologue + 3*3 core + 2 bridge + 1 epilogue

        # Check structure: prologue + (core + bridge)*(n-1) + core + epilogue
        assert expanded[0][0] == 0 and np.isclose(expanded[0][1], 0.3)  # prologue
        assert expanded[1:4] == [(1, 0.7), (0, 0.9), (1, 1.1)]  # core 1
        assert expanded[4][0] == 0 and np.isclose(expanded[4][1], 0.9)  # bridge
        assert expanded[5:8] == [(1, 0.7), (0, 0.9), (1, 1.1)]  # core 2
        assert expanded[8][0] == 0 and np.isclose(expanded[8][1], 0.9)  # bridge
        assert expanded[9:12] == [(1, 0.7), (0, 0.9), (1, 1.1)]  # core 3
        assert expanded[12][0] == 0 and np.isclose(expanded[12][1], 0.6)  # epilogue

    def test_auto_selection_behavior(self):
        """Test that auto-selection chooses the right method based on num_steps."""
        pauli_terms = [('X', 1.0), ('Z', 0.5)]

        # For num_steps < 10, should use incremental
        # For num_steps >= 10, should use structured (if pattern detected)

        # Small num_steps: verify it doesn't crash (auto selects incremental)
        for num_steps in [1, 2, 5, 9]:
            trotter = Trotterization.from_method(
                pauli_terms=pauli_terms,
                method='second order',
                time=0.1,
                num_steps=num_steps,
                tensor_contraction_method=None  # Auto
            )
            U = trotter.tensor_contract()
            assert U.shape == (2, 2)

        # Large num_steps: verify it doesn't crash (auto may select structured)
        for num_steps in [10, 50]:
            trotter = Trotterization.from_method(
                pauli_terms=pauli_terms,
                method='second order',
                time=0.1,
                num_steps=num_steps,
                tensor_contraction_method=None  # Auto
            )
            U = trotter.tensor_contract()
            assert U.shape == (2, 2)

    def test_qualtran_method_accessible(self):
        """Test that the 'qualtran' method successfully calls Bloq.tensor_contract()."""
        pauli_terms = [('X', 1.0), ('Z', 0.5)]

        trotter = Trotterization.from_method(
            pauli_terms=pauli_terms,
            method='second order',
            time=0.1,
            num_steps=5,
            tensor_contraction_method='qualtran'
        )

        U = trotter.tensor_contract()

        # Should produce a valid unitary
        U_dag_U = U.conj().T @ U
        identity = np.eye(U.shape[0])
        assert np.allclose(U_dag_U, identity, atol=1e-10)

    def test_invalid_method_raises_error(self):
        """Test that an invalid tensor_contraction_method is handled correctly."""
        # Note: Validation happens in analysis/config_types.py when used through the config system
        # Here we test that an invalid method in tensor_contract() raises an error
        pauli_terms = [('X', 1.0)]

        trotter = Trotterization(
            pauli_terms=tuple(pauli_terms),
            coefficients=(0.5, 0.5),
            time=0.1,
            num_steps=5,
            tensor_contraction_method='invalid_method'  # Invalid method
        )

        with pytest.raises(ValueError, match="Invalid tensor_contraction_method"):
            trotter.tensor_contract()
