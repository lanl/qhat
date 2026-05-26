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
    expand_ramped_trotterization,
    get_trotterization_coefficients,
)


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
        with pytest.raises(ValueError, match="at least one term"):
            expand_ramped_trotterization(0, [1.0], 1)

    def test_zero_steps_raises_error(self):
        """Test that zero steps raises error."""
        with pytest.raises(ValueError, match="at least one step"):
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
        """Test creation with single term."""
        trotter = Trotterization(
            pauli_terms=(("XYZ", 0.5),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=10
        )
        assert trotter.num_terms == 1
        assert trotter.num_qubits == 3
        assert trotter.num_steps == 10

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
        """Test string representation for first-order."""
        trotter = Trotterization(
            pauli_terms=(("XY", 0.5),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=10
        )
        s = str(trotter)
        assert "first-order" in s
        assert "1 terms" in s
        assert "10 steps" in s

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
        """Test num_steps attribute."""
        trotter = Trotterization(
            pauli_terms=(("X", 1.0),),
            coefficients=(1.0,),
            time=1.0,
            num_steps=42
        )
        assert trotter.num_steps == 42


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
        """Test from_method with fourth-order."""
        trotter = Trotterization.from_method(
            pauli_terms=[("XY", 0.5)],
            method="fourth order",
            time=2.0,
            num_steps=20
        )
        assert len(trotter.coefficients) == 10
        assert trotter.time == 2.0
        assert trotter.num_steps == 20

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
