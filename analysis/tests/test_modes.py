"""
Tests for mode selection and integration of Trotter error coefficient computation.

Includes:
- Mode selection (monte_carlo, exact, auto_exact)
- Backward compatibility
- Integration tests
- Error handling
"""

import numpy as np
from openfermion import QubitOperator
import sys
import os

from .mock_config import mock_config
from qhat.analysis.trotter_coefficients_fast import trotter_error_estimator_fast



# =============================================================================
# Mode Selection Tests
# =============================================================================

def test_backward_compatibility():
    """Default parameters should work unchanged."""
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    terms = [H1, H2]

    # Default should be monte_carlo with auto_exact=False
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config)

    expected_c1 = 1.0
    expected_c2 = 1.0/6.0

    assert abs(c1 - expected_c1) < 1e-6, f"C1: expected {expected_c1}, got {c1}"
    assert abs(c2 - expected_c2) < 1e-6, f"C2: expected {expected_c2}, got {c2}"
    print("✓ Test passed: Backward compatibility")


def test_explicit_monte_carlo():
    """Explicit monte_carlo mode works."""
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    terms = [H1, H2]

    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='monte_carlo')

    assert abs(c1 - 1.0) < 1e-6
    assert abs(c2 - 1.0/6.0) < 1e-6
    print("✓ Test passed: Explicit monte_carlo mode")


def test_exact_mode_small_system():
    """Exact mode works on small system."""
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    terms = [H1, H2]

    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')

    assert abs(c1 - 1.0) < 1e-10, f"C1: expected 1.0, got {c1}"
    assert abs(c2 - 1.0/6.0) < 1e-10, f"C2: expected {1.0/6.0}, got {c2}"
    print("✓ Test passed: Exact mode on small system")


def test_auto_exact_mode():
    """Auto-exact switches to exact on small system."""
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    terms = [H1, H2]

    c1, c2 = trotter_error_estimator_fast(terms, 60, mock_config, mode='monte_carlo', auto_exact=True)

    assert abs(c1 - 1.0) < 1e-10, f"C1: expected 1.0, got {c1}"
    assert abs(c2 - 1.0/6.0) < 1e-10, f"C2: expected {1.0/6.0}, got {c2}"
    print("✓ Test passed: Auto-exact mode")


def test_invalid_mode_raises_error():
    """Invalid mode raises ValueError."""
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    terms = [H1, H2]

    try:
        c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='invalid')
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "Invalid mode" in str(e)
        print("✓ Test passed: Invalid mode raises error")


def test_exact_mode_large_system_raises_error():
    """Exact mode raises error on infeasible system."""
    # Create 500 terms - infeasible for exact
    terms = []
    for i in range(500):
        qubit = i % 50
        pauli = ['X', 'Y', 'Z'][i % 3]
        terms.append(QubitOperator(f'{pauli}{qubit}', 1.0))

    try:
        c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "not feasible" in str(e)
        print("✓ Test passed: Exact mode rejects large system")


def test_monte_carlo_works_large_system():
    """Monte Carlo mode works on large system."""
    # Create 200 terms
    terms = []
    for i in range(200):
        qubit = i % 50
        pauli = ['X', 'Y', 'Z'][i % 3]
        terms.append(QubitOperator(f'{pauli}{qubit}', 1.0))

    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='monte_carlo')

    assert c1 > 0 and c2 >= 0
    print("✓ Test passed: Monte Carlo on large system")


def test_modes_agree_on_small_system():
    """All modes give same answer on small system."""
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    H3 = QubitOperator('Z0', 1.0)
    terms = [H1, H2, H3]

    c1_mc, c2_mc = trotter_error_estimator_fast(terms, 20, mock_config, mode='monte_carlo', auto_exact=False)
    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 20, mock_config, mode='exact')
    c1_auto, c2_auto = trotter_error_estimator_fast(terms, 60, mock_config, mode='monte_carlo', auto_exact=True)

    # Exact and auto-exact should match perfectly
    assert abs(c1_auto - c1_ex) < 1e-10, "Auto-exact C1 should match exact"
    assert abs(c2_auto - c2_ex) < 1e-10, "Auto-exact C2 should match exact"

    # Monte Carlo should be close
    rel_err_c1 = abs(c1_mc - c1_ex) / c1_ex if c1_ex > 0 else 0
    rel_err_c2 = abs(c2_mc - c2_ex) / c2_ex if c2_ex > 0 else 0
    assert rel_err_c1 < 0.1, f"MC C1 differs by {100*rel_err_c1:.1f}%"
    assert rel_err_c2 < 0.1, f"MC C2 differs by {100*rel_err_c2:.1f}%"

    print("✓ Test passed: All modes agree")


# =============================================================================
# Integration Tests
# =============================================================================

def test_complete_workflow():
    """Test complete workflow with different configurations."""
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    H3 = QubitOperator('Z0', 1.0)
    H4 = QubitOperator('X1', 1.0)
    H5 = QubitOperator('Y1', 1.0)
    terms = [H1, H2, H3, H4, H5]

    # Test all modes work
    c1_default, c2_default = trotter_error_estimator_fast(terms, 10, mock_config)
    c1_mc, c2_mc = trotter_error_estimator_fast(terms, 10, mock_config, mode='monte_carlo')
    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')
    c1_auto, c2_auto = trotter_error_estimator_fast(terms, 60, mock_config, mode='monte_carlo', auto_exact=True)

    # All should return positive values
    assert c1_default > 0 and c2_default > 0
    assert c1_mc > 0 and c2_mc > 0
    assert c1_ex > 0 and c2_ex > 0
    assert c1_auto > 0 and c2_auto > 0

    # Exact and auto-exact should match
    assert abs(c1_ex - c1_auto) < 1e-10
    assert abs(c2_ex - c2_auto) < 1e-10

    print("✓ Test passed: Complete workflow")


def test_larger_system_integration():
    """Test with larger but still feasible system."""
    np.random.seed(42)
    terms = []
    for i in range(20):
        qubit = i % 10
        pauli = ['X', 'Y', 'Z'][i % 3]
        coeff = np.random.uniform(0.5, 2.0)
        terms.append(QubitOperator(f'{pauli}{qubit}', coeff))

    # All modes should work
    c1_mc, c2_mc = trotter_error_estimator_fast(terms, 20, mock_config, mode='monte_carlo', auto_exact=False)
    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 20, mock_config, mode='exact')
    c1_auto, c2_auto = trotter_error_estimator_fast(terms, 60, mock_config, mode='monte_carlo', auto_exact=True)

    # Exact modes should match
    assert abs(c1_ex - c1_auto) < 1e-10
    assert abs(c2_ex - c2_auto) < 1e-10

    print("✓ Test passed: Larger system integration")


def test_edge_case_integration():
    """Test edge cases work in integration."""
    # N=1
    terms = [QubitOperator('X0', 1.0)]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')
    assert c1 == 0.0 and c2 == 0.0

    # N=2
    terms = [QubitOperator('X0', 1.0), QubitOperator('Y0', 1.0)]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')
    assert c1 > 0

    # All commuting
    terms = [QubitOperator('X0', 1.0), QubitOperator('X1', 1.0), QubitOperator('X2', 1.0)]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')
    assert c1 == 0.0 and c2 == 0.0

    # Complex coefficients
    terms = [QubitOperator('X0', 1.0+1.0j), QubitOperator('Y0', 1.0-1.0j)]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mock_config, mode='exact')
    assert c1 > 0

    print("✓ Test passed: Edge case integration")


# =============================================================================
# Test Runner
# =============================================================================

if __name__ == "__main__":
    print("="*70)
    print("MODE SELECTION AND INTEGRATION TESTS")
    print("="*70)
    print()

    print("Mode Selection:")
    test_backward_compatibility()
    test_explicit_monte_carlo()
    test_exact_mode_small_system()
    test_auto_exact_mode()
    test_invalid_mode_raises_error()
    test_exact_mode_large_system_raises_error()
    test_monte_carlo_works_large_system()
    test_modes_agree_on_small_system()
    print()

    print("Integration:")
    test_complete_workflow()
    test_larger_system_integration()
    test_edge_case_integration()
    print()

    print("="*70)
    print("✅ ALL MODE AND INTEGRATION TESTS PASSED!")
    print("="*70)
