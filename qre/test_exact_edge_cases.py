"""
Test edge cases and boundary conditions for exact computation.
"""

import numpy as np
from openfermion import QubitOperator
from trotter_coefficients_fast import trotter_error_estimator_fast


def test_single_term():
    """N=1: Single term Hamiltonian - no pairs or triples possible."""
    terms = [QubitOperator('X0', 1.0)]

    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')

    assert c1 == 0.0, f"Expected C1=0, got {c1}"
    assert c2 == 0.0, f"Expected C2=0, got {c2}"
    print("✓ Test passed: N=1 (single term)")


def test_two_terms_minimum_pairs():
    """N=2: Minimum for pairs, no triples."""
    terms = [
        QubitOperator('X0', 1.0),
        QubitOperator('Y0', 1.0)
    ]

    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')

    # Should have 1 pair, 0 triples
    assert c1 > 0, f"Expected C1>0 for anticommuting terms, got {c1}"
    assert c2 > 0, f"Expected C2>0 (from C22), got {c2}"
    print("✓ Test passed: N=2 (minimum for pairs)")


def test_three_terms_minimum_triples():
    """N=3: Minimum for triples."""
    terms = [
        QubitOperator('X0', 1.0),
        QubitOperator('Y0', 1.0),
        QubitOperator('Z0', 1.0)
    ]

    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')

    # Should have 3 pairs, 1 triple
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

    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')

    assert abs(c1) < 1e-10, f"Expected C1≈0 for commuting terms, got {c1}"
    assert abs(c2) < 1e-10, f"Expected C2≈0 for commuting terms, got {c2}"
    print("✓ Test passed: All commuting terms")


def test_mixed_commuting_anticommuting():
    """Mix of commuting and anticommuting terms."""
    terms = [
        QubitOperator('X0', 1.0),  # anticommutes with Y0
        QubitOperator('Y0', 1.0),  # anticommutes with X0
        QubitOperator('X1', 1.0),  # commutes with above
        QubitOperator('Y1', 1.0)   # anticommutes with X1, commutes with qubit 0 ops
    ]

    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')

    # Should have contributions from (X0,Y0) and (X1,Y1) pairs
    # Expected: 2 pairs contribute, C1 = 2+2 = 4, so C1_est = 2.0
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

    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')

    # Zero coefficient means its commutators have zero norm
    # Only Y0-Z0 pair contributes: ||[Y,Z]|| = 2
    # C1 = 2, so C1_est = 1.0
    assert abs(c1 - 1.0) < 1e-10, f"Expected C1=1.0, got {c1}"
    print("✓ Test passed: Zero coefficients")


def test_complex_coefficients():
    """Complex coefficients should work correctly."""
    terms = [
        QubitOperator('X0', 1.0 + 1.0j),
        QubitOperator('Y0', 2.0 - 0.5j)
    ]

    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')

    # ||[aX, bY]|| = 2|ab| where a=(1+j), b=(2-0.5j)
    # ab = (1+j)(2-0.5j) = 2 - 0.5j + 2j + 0.5 = 2.5 + 1.5j
    # |ab| = sqrt(2.5^2 + 1.5^2) = sqrt(8.5) ≈ 2.9155
    # C1 = 2|ab| ≈ 5.831, so C1_est ≈ 2.9155
    expected = np.sqrt(8.5)
    assert abs(c1 - expected) < 1e-10, f"Expected C1={expected}, got {c1}"
    print("✓ Test passed: Complex coefficients")


def test_large_coefficients():
    """Large coefficients should scale correctly."""
    terms = [
        QubitOperator('X0', 100.0),
        QubitOperator('Y0', 200.0)
    ]

    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')

    # ||[100X, 200Y]|| = 2|100*200| = 40000
    # C1 = 40000, so C1_est = 20000
    assert abs(c1 - 20000.0) < 1e-8, f"Expected C1=20000, got {c1}"
    print("✓ Test passed: Large coefficients")


if __name__ == "__main__":
    print("="*70)
    print("EDGE CASE TESTS")
    print("="*70)
    print()
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
    print("✅ ALL EDGE CASE TESTS PASSED!")
    print("="*70)
