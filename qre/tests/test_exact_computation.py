"""
Tests for exact computation of Trotter error coefficients.

Includes:
- Analytical validation with hand-calculable test cases
- Edge cases and boundary conditions
- Basic function validation
"""

import numpy as np
from openfermion import QubitOperator
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from trotter_coefficients_fast import (
    trotter_error_estimator_fast,
    preprocess_pauli_terms,
    generate_all_pairs,
    generate_all_triples,
    compute_C1_exact,
    compute_C21_exact,
    compute_C22_exact
)
from .mock_config import mock_config


# =============================================================================
# Basic Function Tests
# =============================================================================

def test_generate_pairs():
    """Test deterministic pair generation."""
    pairs = generate_all_pairs(5)
    assert len(pairs) == 10, f"Expected 10 pairs, got {len(pairs)}"
    # Verify all pairs are unique
    pair_set = set(tuple(p) for p in pairs)
    assert len(pair_set) == 10, "Pairs should be unique"
    print("✓ Test passed: generate_all_pairs")


def test_generate_triples():
    """Test deterministic triple generation."""
    triples = generate_all_triples(5)
    expected = sum(range(4))  # k=0: 6, k=1: 3, k=2: 1, k=3: 0 = 10
    assert len(triples) == 10, f"Expected 10 triples, got {len(triples)}"
    print("✓ Test passed: generate_all_triples")


# =============================================================================
# Analytical Validation Tests
# =============================================================================

def test_two_anticommuting_paulis():
    """
    H = X₀ + Y₀

    Expected: C1 = 2, C21 = 0, C22 = 4
    Therefore: C1_est = 1.0, C2_est = 1/6
    """
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    terms = [H1, H2]

    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')
    c1_auto, c2_auto = trotter_error_estimator_fast(terms, 60, mock_config, mode='monte_carlo', auto_exact=True)

    expected_c1 = 1.0
    expected_c2 = 1.0/6.0

    assert abs(c1_ex - expected_c1) < 1e-10, f"Exact C1: got {c1_ex}, expected {expected_c1}"
    assert abs(c2_ex - expected_c2) < 1e-10, f"Exact C2: got {c2_ex}, expected {expected_c2}"
    assert abs(c1_auto - expected_c1) < 1e-10, f"Auto C1: got {c1_auto}, expected {expected_c1}"
    assert abs(c2_auto - expected_c2) < 1e-10, f"Auto C2: got {c2_auto}, expected {expected_c2}"

    print("✓ Test passed: Two anticommuting Paulis")


def test_three_mutually_anticommuting():
    """
    H = X₀ + Y₀ + Z₀

    Expected: C1 = 6, C21 = 0, C22 = 12
    Therefore: C1_est = 3.0, C2_est = 0.5
    """
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    H3 = QubitOperator('Z0', 1.0)
    terms = [H1, H2, H3]

    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    expected_c1 = 3.0
    expected_c2 = 0.5

    assert abs(c1_ex - expected_c1) < 1e-10, f"C1: got {c1_ex}, expected {expected_c1}"
    assert abs(c2_ex - expected_c2) < 1e-10, f"C2: got {c2_ex}, expected {expected_c2}"

    print("✓ Test passed: Three mutually anticommuting Paulis")


def test_commuting_paulis():
    """
    H = X₀ + X₁ (commuting)

    Expected: C1 = 0, C2 = 0
    """
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('X1', 1.0)
    terms = [H1, H2]

    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    assert abs(c1_ex) < 1e-10, f"C1: got {c1_ex}, expected 0.0"
    assert abs(c2_ex) < 1e-10, f"C2: got {c2_ex}, expected 0.0"

    print("✓ Test passed: Commuting Paulis")


def test_non_unit_coefficients():
    """
    H = 2X₀ + 3Y₀

    Expected: C1 = 12, C22 = 48
    Therefore: C1_est = 6.0, C2_est = 2.0
    """
    H1 = QubitOperator('X0', 2.0)
    H2 = QubitOperator('Y0', 3.0)
    terms = [H1, H2]

    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    expected_c1 = 6.0
    expected_c2 = 2.0

    assert abs(c1_ex - expected_c1) < 1e-10, f"C1: got {c1_ex}, expected {expected_c1}"
    assert abs(c2_ex - expected_c2) < 1e-10, f"C2: got {c2_ex}, expected {expected_c2}"

    print("✓ Test passed: Non-unit coefficients")


# =============================================================================
# Edge Cases
# =============================================================================

def test_single_term():
    """N=1: Single term - no pairs or triples."""
    terms = [QubitOperator('X0', 1.0)]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    assert c1 == 0.0, f"Expected C1=0, got {c1}"
    assert c2 == 0.0, f"Expected C2=0, got {c2}"
    print("✓ Test passed: N=1 (single term)")


def test_two_terms_minimum_pairs():
    """N=2: Minimum for pairs, no triples."""
    terms = [QubitOperator('X0', 1.0), QubitOperator('Y0', 1.0)]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    assert c1 > 0, f"Expected C1>0, got {c1}"
    assert c2 > 0, f"Expected C2>0, got {c2}"
    print("✓ Test passed: N=2 (minimum for pairs)")


def test_three_terms_minimum_triples():
    """N=3: Minimum for triples."""
    terms = [
        QubitOperator('X0', 1.0),
        QubitOperator('Y0', 1.0),
        QubitOperator('Z0', 1.0)
    ]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    assert c1 > 0, f"Expected C1>0, got {c1}"
    assert c2 > 0, f"Expected C2>0, got {c2}"
    print("✓ Test passed: N=3 (minimum for triples)")


def test_all_commute():
    """All terms commute - should give zeros."""
    terms = [
        QubitOperator('X0', 1.0),
        QubitOperator('X1', 1.0),
        QubitOperator('X2', 1.0),
        QubitOperator('X3', 1.0)
    ]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    assert abs(c1) < 1e-10, f"Expected C1≈0, got {c1}"
    assert abs(c2) < 1e-10, f"Expected C2≈0, got {c2}"
    print("✓ Test passed: All commuting terms")


def test_mixed_commuting_anticommuting():
    """Mix of commuting and anticommuting terms."""
    terms = [
        QubitOperator('X0', 1.0),  # anticommutes with Y0
        QubitOperator('Y0', 1.0),  # anticommutes with X0
        QubitOperator('X1', 1.0),  # commutes with above
        QubitOperator('Y1', 1.0)   # anticommutes with X1
    ]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    # Expected: 2 pairs contribute (X0,Y0) and (X1,Y1)
    # C1 = 2+2 = 4, so C1_est = 2.0
    assert abs(c1 - 2.0) < 1e-10, f"Expected C1=2.0, got {c1}"
    assert c2 > 0, f"Expected C2>0, got {c2}"
    print("✓ Test passed: Mixed commuting/anticommuting")


def test_zero_coefficients():
    """Terms with zero coefficients should not contribute."""
    terms = [
        QubitOperator('X0', 0.0),
        QubitOperator('Y0', 1.0),
        QubitOperator('Z0', 1.0)
    ]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    # Only Y0-Z0 pair contributes: C1 = 2, so C1_est = 1.0
    assert abs(c1 - 1.0) < 1e-10, f"Expected C1=1.0, got {c1}"
    print("✓ Test passed: Zero coefficients")


def test_complex_coefficients():
    """Complex coefficients should work correctly."""
    terms = [
        QubitOperator('X0', 1.0 + 1.0j),
        QubitOperator('Y0', 2.0 - 0.5j)
    ]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    # |ab| = sqrt(2.5^2 + 1.5^2) = sqrt(8.5)
    expected = np.sqrt(8.5)
    assert abs(c1 - expected) < 1e-10, f"Expected C1={expected}, got {c1}"
    print("✓ Test passed: Complex coefficients")


def test_large_coefficients():
    """Large coefficients should scale correctly."""
    terms = [
        QubitOperator('X0', 100.0),
        QubitOperator('Y0', 200.0)
    ]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    # C1 = 2|100*200| = 40000, so C1_est = 20000
    assert abs(c1 - 20000.0) < 1e-8, f"Expected C1=20000, got {c1}"
    print("✓ Test passed: Large coefficients")


# =============================================================================
# Test Runner
# =============================================================================

if __name__ == "__main__":
    print("="*70)
    print("EXACT COMPUTATION TESTS")
    print("="*70)
    print()

    # Basic function tests
    print("Basic Function Tests:")
    test_generate_pairs()
    test_generate_triples()
    print()

    # Analytical tests
    print("Analytical Validation:")
    test_two_anticommuting_paulis()
    test_three_mutually_anticommuting()
    test_commuting_paulis()
    test_non_unit_coefficients()
    print()

    # Edge cases
    print("Edge Cases:")
    test_single_term()
    test_two_terms_minimum_pairs()
    test_three_terms_minimum_triples()
    test_all_commute()
    test_mixed_commuting_anticommuting()
    test_zero_coefficients()
    test_complex_coefficients()
    test_large_coefficients()
    print()

    print("="*70)
    print("✅ ALL EXACT COMPUTATION TESTS PASSED!")
    print("="*70)
