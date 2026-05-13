"""
Comprehensive test suite for CommutingPauliStringEvolution Bloq

Tests commutativity checking, composition, resource estimation, and unitary generation.
"""

import numpy as np
import pytest
from scipy.linalg import expm

from commuting_pauli_string_evolution import (
    CommutingPauliStringEvolution,
    pauli_strings_commute,
    all_commute,
)
from qualtran.cirq_interop.t_complexity_protocol import t_complexity


# ==================================================================================
# Helper functions for analytical comparisons
# ==================================================================================

def get_pauli_matrix(char):
    """Return the 2x2 Pauli matrix for a given character."""
    if char == 'I':
        return np.array([[1, 0], [0, 1]], dtype=complex)
    elif char == 'X':
        return np.array([[0, 1], [1, 0]], dtype=complex)
    elif char == 'Y':
        return np.array([[0, -1j], [1j, 0]], dtype=complex)
    elif char == 'Z':
        return np.array([[1, 0], [0, -1]], dtype=complex)
    else:
        raise ValueError(f"Invalid Pauli character: {char}")


def pauli_string_to_matrix(pauli_string):
    """Convert a Pauli string to its matrix representation via tensor product."""
    result = get_pauli_matrix(pauli_string[0])
    for char in pauli_string[1:]:
        result = np.kron(result, get_pauli_matrix(char))
    return result


def analytical_commuting_evolution(pauli_terms, time, hbar=1.0):
    """Compute the exact evolution for commuting terms analytically.

    Since terms commute: exp(-i*(A+B)*t) = exp(-i*A*t) * exp(-i*B*t)
    """
    if len(pauli_terms) == 0:
        raise ValueError("Must have at least one term")

    num_qubits = len(pauli_terms[0][0])
    dim = 2 ** num_qubits

    # Start with identity
    U = np.eye(dim, dtype=complex)

    # Multiply by each individual evolution
    for pauli_string, coefficient in pauli_terms:
        P = pauli_string_to_matrix(pauli_string)
        U_term = expm(-1j * coefficient * P * time / hbar)
        U = U_term @ U

    return U


# ==================================================================================
# Test: Commutativity checking
# ==================================================================================

class TestCommutativity:
    """Test Pauli string commutativity checking."""

    def test_identical_strings_commute(self):
        """Test that identical strings commute."""
        assert pauli_strings_commute("XYZ", "XYZ")
        assert pauli_strings_commute("I", "I")
        assert pauli_strings_commute("XYZI", "XYZI")

    def test_identity_commutes_with_everything(self):
        """Test that identity commutes with all strings."""
        assert pauli_strings_commute("III", "XYZ")
        assert pauli_strings_commute("XYZ", "III")
        # When one string has I at a position, that position doesn't contribute to differences
        assert pauli_strings_commute("IXI", "III")

    def test_even_differences_commute(self):
        """Test strings with even number of differences commute."""
        # 2 differences: commute
        assert pauli_strings_commute("XY", "YX")
        assert pauli_strings_commute("XZ", "ZX")
        assert pauli_strings_commute("YZ", "ZY")

        # 0 differences (same non-identity positions): commute
        assert pauli_strings_commute("XIX", "XYX")  # Only middle differs, so effectively 0
        assert pauli_strings_commute("ZII", "ZII")

    def test_odd_differences_dont_commute(self):
        """Test strings with odd number of differences don't commute."""
        # 1 difference: don't commute
        assert not pauli_strings_commute("XY", "ZY")  # X vs Z at pos 0
        assert not pauli_strings_commute("ZI", "YI")  # Z vs Y at pos 0
        assert not pauli_strings_commute("IXI", "IYI")  # X vs Y at pos 1

        # 3 differences: don't commute
        assert not pauli_strings_commute("XYZ", "YZX")  # X vs Y, Y vs Z, Z vs X - all 3 diff

    def test_different_lengths_raise_error(self):
        """Test that different length strings raise an error."""
        with pytest.raises(ValueError, match="same length"):
            pauli_strings_commute("XY", "XYZ")

    def test_all_commute_empty(self):
        """Test all_commute with empty list."""
        assert all_commute([])

    def test_all_commute_single(self):
        """Test all_commute with single string."""
        assert all_commute(["XYZ"])

    def test_all_commute_pairwise_true(self):
        """Test all_commute when all pairs commute."""
        # All identity: commute
        assert all_commute(["III", "III", "III"])

        # Mix of commuting strings
        assert all_commute(["XII", "IXI", "IIX"])  # All commute pairwise

    def test_all_commute_pairwise_false(self):
        """Test all_commute when at least one pair doesn't commute."""
        assert not all_commute(["XI", "YI"])  # X and Y on same qubit don't commute


# ==================================================================================
# Test: Instantiation and validation
# ==================================================================================

class TestInstantiation:
    """Test CommutingPauliStringEvolution instantiation and validation."""

    def test_valid_single_term(self):
        """Test creation with single term."""
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XYZ", 0.5),),
            time=1.0
        )
        assert cpse.num_terms == 1
        assert cpse.num_qubits == 3

    def test_valid_multiple_commuting_terms(self):
        """Test creation with multiple commuting terms."""
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XII", 0.5), ("IXI", 0.3), ("IIX", 0.2)),
            time=1.0
        )
        assert cpse.num_terms == 3
        assert cpse.num_qubits == 3

    def test_from_pauli_string_iterable(self):
        """Test creation from iterable (generator format)."""
        terms = [("XZ", 0.5), ("ZX", 0.3)]
        cpse = CommutingPauliStringEvolution.from_pauli_string_iterable(
            terms, time=1.0
        )
        assert cpse.num_terms == 2

    def test_empty_terms_raises_error(self):
        """Test that empty terms list raises error."""
        with pytest.raises(ValueError, match="at least one"):
            CommutingPauliStringEvolution(pauli_terms=(), time=1.0)

    def test_non_commuting_terms_raise_error(self):
        """Test that non-commuting terms raise an error."""
        with pytest.raises(ValueError, match="Not all Pauli strings commute"):
            CommutingPauliStringEvolution(
                pauli_terms=(("XI", 0.5), ("YI", 0.3)),  # Don't commute
                time=1.0
            )

    def test_different_length_strings_raise_error(self):
        """Test that different length strings raise an error."""
        with pytest.raises(ValueError, match="same length"):
            CommutingPauliStringEvolution(
                pauli_terms=(("XI", 0.5), ("XYZ", 0.3)),
                time=1.0
            )

    def test_invalid_character_raises_error(self):
        """Test that invalid characters raise an error."""
        with pytest.raises(ValueError, match="Invalid character"):
            CommutingPauliStringEvolution(
                pauli_terms=(("XAZ", 0.5),),
                time=1.0
            )


# ==================================================================================
# Test: Signature
# ==================================================================================

class TestSignature:
    """Test the Bloq signature."""

    def test_signature_structure(self):
        """Test that signature has correct structure."""
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XYZ", 0.5), ("ZYX", 0.3)),
            time=1.0
        )
        sig = cpse.signature
        assert 'q' in [reg.name for reg in sig]
        q_reg = [reg for reg in sig if reg.name == 'q'][0]
        assert q_reg.shape == (3,)
        assert q_reg.bitsize == 1


# ==================================================================================
# Test: Unitary matrix generation
# ==================================================================================

class TestUnitaryGeneration:
    """Test unitary matrix generation via tensor_contract."""

    def test_single_term_matches_analytical(self):
        """Test single term matches analytical result."""
        terms = (("XY", 0.5),)
        cpse = CommutingPauliStringEvolution(pauli_terms=terms, time=1.0)
        U = cpse.tensor_contract()

        expected = analytical_commuting_evolution(terms, 1.0)
        assert np.allclose(U, expected)

    def test_two_commuting_terms_match_analytical(self):
        """Test two commuting terms match analytical result."""
        terms = (("XZ", 0.5), ("ZX", 0.3))  # These commute
        cpse = CommutingPauliStringEvolution(pauli_terms=terms, time=1.0)
        U = cpse.tensor_contract()

        expected = analytical_commuting_evolution(terms, 1.0)
        assert np.allclose(U, expected)

    def test_three_commuting_terms(self):
        """Test three commuting terms."""
        terms = (("XII", 0.5), ("IXI", 0.3), ("IIX", 0.2))
        cpse = CommutingPauliStringEvolution(pauli_terms=terms, time=1.0)
        U = cpse.tensor_contract()

        expected = analytical_commuting_evolution(terms, 1.0)
        assert np.allclose(U, expected)

    def test_unitarity(self):
        """Test that generated matrix is unitary."""
        terms = (("XY", 0.5), ("YX", 0.3))
        cpse = CommutingPauliStringEvolution(pauli_terms=terms, time=1.0)
        U = cpse.tensor_contract()

        # Check U† U = I
        U_dag_U = U.conj().T @ U
        identity = np.eye(U.shape[0])
        assert np.allclose(U_dag_U, identity), "Matrix is not unitary"

    def test_commuting_vs_sum_hamiltonian(self):
        """Test that commuting evolution equals evolution under sum Hamiltonian."""
        # For commuting terms: exp(-i*(A+B)*t) = exp(-i*A*t) * exp(-i*B*t)
        terms = (("ZI", 0.5), ("IZ", 0.3))  # Commute
        cpse = CommutingPauliStringEvolution(pauli_terms=terms, time=1.0)
        U_commuting = cpse.tensor_contract()

        # Compute evolution under the sum directly
        H_total = 0.5 * pauli_string_to_matrix("ZI") + 0.3 * pauli_string_to_matrix("IZ")
        U_sum = expm(-1j * H_total * 1.0)

        assert np.allclose(U_commuting, U_sum)

    def test_custom_hbar(self):
        """Test evolution with custom hbar."""
        terms = (("X", 0.5),)
        cpse = CommutingPauliStringEvolution(pauli_terms=terms, time=1.0, hbar=2.0)
        U = cpse.tensor_contract()

        expected = analytical_commuting_evolution(terms, 1.0, hbar=2.0)
        assert np.allclose(U, expected)


# ==================================================================================
# Test: Resource estimation
# ==================================================================================

class TestResourceEstimation:
    """Test resource estimation via T-complexity."""

    def test_returns_tcomplexity_object(self):
        """Test that _t_complexity_ returns a TComplexity object."""
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XY", 0.5), ("YX", 0.3)),
            time=1.0
        )
        tc = t_complexity(cpse)

        assert hasattr(tc, 't')
        assert hasattr(tc, 'rotations')
        assert hasattr(tc, 'clifford')

    def test_non_negative_counts(self):
        """Test that all resource counts are non-negative."""
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XY", 0.5),),
            time=1.0
        )
        tc = t_complexity(cpse)

        assert tc.t >= 0
        assert tc.rotations >= 0
        assert tc.clifford >= 0

    def test_complexity_scales_with_terms(self):
        """Test that complexity increases with more terms."""
        # Single term
        cpse1 = CommutingPauliStringEvolution(
            pauli_terms=(("XII", 0.5),),
            time=1.0
        )
        tc1 = t_complexity(cpse1)

        # Two terms
        cpse2 = CommutingPauliStringEvolution(
            pauli_terms=(("XII", 0.5), ("IXI", 0.3)),
            time=1.0
        )
        tc2 = t_complexity(cpse2)

        # Should generally have more resources with more terms
        assert tc2.rotations >= tc1.rotations or tc2.t >= tc1.t


# ==================================================================================
# Test: Decomposition
# ==================================================================================

class TestDecomposition:
    """Test the decomposition into PauliStringEvolution bloqs."""

    def test_decompose_single_term(self):
        """Test decomposition with single term."""
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XY", 0.5),),
            time=1.0
        )
        # Should be able to decompose
        cbloq = cpse.decompose_bloq()
        assert cbloq is not None

    def test_decompose_multiple_terms(self):
        """Test decomposition with multiple terms."""
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XII", 0.5), ("IXI", 0.3)),
            time=1.0
        )
        cbloq = cpse.decompose_bloq()
        assert cbloq is not None


# ==================================================================================
# Test: String representations
# ==================================================================================

class TestStringRepresentations:
    """Test string and repr methods."""

    def test_str_representation(self):
        """Test string representation."""
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XY", 0.5), ("YX", 0.3)),
            time=1.0
        )
        s = str(cpse)
        assert "XY" in s
        assert "YX" in s
        assert "0.5" in s
        assert "0.3" in s

    def test_repr_representation(self):
        """Test repr representation."""
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XY", 0.5),),
            time=1.0
        )
        r = repr(cpse)
        assert "CommutingPauliStringEvolution" in r


# ==================================================================================
# Test: Properties
# ==================================================================================

class TestProperties:
    """Test various properties and accessors."""

    def test_num_terms(self):
        """Test num_terms property."""
        # Use strings that all pairwise commute
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XII", 0.5), ("IXI", 0.3), ("IIX", 0.2)),  # All commute
            time=1.0
        )
        assert cpse.num_terms == 3

    def test_num_qubits(self):
        """Test num_qubits property."""
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XYZX", 0.5),),
            time=1.0
        )
        assert cpse.num_qubits == 4
