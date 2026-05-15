"""
Test exact computation against analytical calculations.

These tests use simple Hamiltonians where error coefficients can be computed by hand.
"""

import numpy as np
from openfermion import QubitOperator
from trotter_coefficients_fast import trotter_error_estimator_fast


def test_analytical_two_paulis_same_qubit():
    """
    Test Case 1: Two anticommuting Paulis on same qubit

    H = X₀ + Y₀

    Analytical calculation:
    - C1 = ||[X₀, Y₀]|| = ||2iZ₀|| = 2
    - C21 = 0 (only 2 terms, no valid triples)
    - C22 = ||[X₀, [X₀, Y₀]]|| = ||[X₀, 2iZ₀]|| = 4
    - Therefore: C1_est = C1/2 = 1.0, C2_est = 0/12 + 4/24 = 1/6
    """
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    terms = [H1, H2]

    # Test all three modes
    c1_mc, c2_mc = trotter_error_estimator_fast(terms, 10, mode='monte_carlo')
    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 10, mode='exact')
    c1_auto, c2_auto = trotter_error_estimator_fast(terms, 60, mode='monte_carlo', auto_exact=True)

    # Expected values
    expected_c1 = 1.0  # C1/2 = 2/2
    expected_c2 = 1.0/6.0  # 0/12 + 4/24

    # All modes should give same answer (exact modes should be exact, MC should be close)
    assert abs(c1_ex - expected_c1) < 1e-10, f"Exact C1: got {c1_ex}, expected {expected_c1}"
    assert abs(c2_ex - expected_c2) < 1e-10, f"Exact C2: got {c2_ex}, expected {expected_c2}"
    assert abs(c1_mc - expected_c1) < 1e-6, f"Monte Carlo C1: got {c1_mc}, expected {expected_c1}"
    assert abs(c2_mc - expected_c2) < 1e-6, f"Monte Carlo C2: got {c2_mc}, expected {expected_c2}"
    assert abs(c1_auto - expected_c1) < 1e-10, f"Auto C1: got {c1_auto}, expected {expected_c1}"
    assert abs(c2_auto - expected_c2) < 1e-10, f"Auto C2: got {c2_auto}, expected {expected_c2}"

    print("✓ Test 1 passed: Two anticommuting Paulis")


def test_analytical_three_paulis():
    """
    Test Case 2: Three mutually anticommuting Paulis

    H = X₀ + Y₀ + Z₀

    Analytical calculation:
    - Pairs: (X,Y), (X,Z), (Y,Z) - all anticommute
    - C1 = 2 + 2 + 2 = 6
    - C21: All triples give 0 (inner commutators commute with outer)
    - C22: k<j
      - k=0, j=1: [X₀, [X₀, Y₀]] = [X₀, 2iZ₀] = 4
      - k=0, j=2: [X₀, [X₀, Z₀]] = [X₀, -2iY₀] = 4
      - k=1, j=2: [Y₀, [Y₀, Z₀]] = [Y₀, 2iX₀] = 4
    - C22 = 4 + 4 + 4 = 12
    - Therefore: C1_est = 6/2 = 3.0, C2_est = 0/12 + 12/24 = 0.5
    """
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    H3 = QubitOperator('Z0', 1.0)
    terms = [H1, H2, H3]

    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 10, mode='exact')

    expected_c1 = 3.0
    expected_c2 = 0.5

    assert abs(c1_ex - expected_c1) < 1e-10, f"C1: got {c1_ex}, expected {expected_c1}"
    assert abs(c2_ex - expected_c2) < 1e-10, f"C2: got {c2_ex}, expected {expected_c2}"

    print("✓ Test 2 passed: Three mutually anticommuting Paulis")


def test_analytical_commuting_paulis():
    """
    Test Case 3: Commuting Paulis (different qubits)

    H = X₀ + X₁

    All terms commute, so:
    - C1 = 0, C21 = 0, C22 = 0
    - Therefore: C1_est = 0, C2_est = 0
    """
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('X1', 1.0)
    terms = [H1, H2]

    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 10, mode='exact')

    expected_c1 = 0.0
    expected_c2 = 0.0

    assert abs(c1_ex - expected_c1) < 1e-10, f"C1: got {c1_ex}, expected {expected_c1}"
    assert abs(c2_ex - expected_c2) < 1e-10, f"C2: got {c2_ex}, expected {expected_c2}"

    print("✓ Test 3 passed: Commuting Paulis")


def test_analytical_with_coefficients():
    """
    Test Case 4: Non-unit coefficients

    H = 2·X₀ + 3·Y₀

    Analytical calculation:
    - [2X₀, 3Y₀] = 6[X₀, Y₀] = 12iZ₀
    - ||[2X₀, 3Y₀]|| = 2|2·3| = 12
    - C1 = 12
    - C22 = ||[2X₀, [2X₀, 3Y₀]]|| = ||[2X₀, 12iZ₀]|| = 24i||[X₀, Z₀]|| = 24i·2i||Y₀|| = 48
    - Therefore: C1_est = 12/2 = 6.0, C2_est = 0/12 + 48/24 = 2.0
    """
    H1 = QubitOperator('X0', 2.0)
    H2 = QubitOperator('Y0', 3.0)
    terms = [H1, H2]

    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 10, mode='exact')

    expected_c1 = 6.0
    expected_c2 = 2.0

    assert abs(c1_ex - expected_c1) < 1e-10, f"C1: got {c1_ex}, expected {expected_c1}"
    assert abs(c2_ex - expected_c2) < 1e-10, f"C2: got {c2_ex}, expected {expected_c2}"

    print("✓ Test 4 passed: Non-unit coefficients")


if __name__ == "__main__":
    print("="*70)
    print("ANALYTICAL VALIDATION TESTS")
    print("="*70)
    print()
    test_analytical_two_paulis_same_qubit()
    test_analytical_three_paulis()
    test_analytical_commuting_paulis()
    test_analytical_with_coefficients()
    print()
    print("="*70)
    print("✅ ALL ANALYTICAL TESTS PASSED!")
    print("="*70)
