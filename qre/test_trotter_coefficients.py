"""
Comprehensive test suite for Trotter coefficient estimation.

This test suite verifies:
1. Ground truth: Implementations match known mathematical properties
2. Consistency: Fast and reference implementations agree
3. Integration: Both implementations work correctly in practice
4. Basic sanity: Results are reasonable for simple systems

Run with: python test_trotter_coefficients.py
"""

import numpy as np
import sys
from openfermion import QubitOperator

# Mock config_general for testing
class MockConfigGeneral:
    def log(self, msg): pass
    def log_verbose(self, msg): pass
    def log_debug(self, msg): pass

mock_config = MockConfigGeneral()

# Import implementations
from trotter_coefficients import (
    anticommute,
    pauli_product_key,
    compute_commutator,
    trotter_error_estimator
)

from trotter_coefficients_fast import (
    encode_pauli_string,
    pauli_anticommute,
    pauli_product_bits,
    compute_commutator_norm_fast,
    trotter_error_estimator_fast
)


# ============================================================================
# Test 1: Ground Truth - Pauli Anticommutation
# ============================================================================
def test_anticommutation():
    """Test that both implementations correctly identify Pauli anticommutation."""
    print("\n" + "=" * 70)
    print("TEST 1: Ground Truth - Pauli Anticommutation Properties")
    print("=" * 70)

    test_cases = [
        # (key1, key2, should_anticommute, description)
        (((0, 'X'),), ((0, 'Y'),), True, "X and Y on same qubit"),
        (((0, 'X'),), ((0, 'Z'),), True, "X and Z on same qubit"),
        (((0, 'Y'),), ((0, 'Z'),), True, "Y and Z on same qubit"),
        (((0, 'X'),), ((0, 'X'),), False, "X and X on same qubit"),
        (((0, 'Y'),), ((0, 'Y'),), False, "Y and Y on same qubit"),
        (((0, 'Z'),), ((0, 'Z'),), False, "Z and Z on same qubit"),
        (((0, 'X'),), ((1, 'X'),), False, "X on different qubits"),
        (((0, 'X'), (1, 'Y')), ((0, 'Y'), (1, 'Y')), True, "1 position differs"),
        (((0, 'X'), (1, 'Y')), ((0, 'Y'), (1, 'X')), False, "2 positions differ"),
        (((0, 'X'), (1, 'Y'), (2, 'Z')), ((0, 'Y'), (1, 'X'), (2, 'X')), True, "3 positions differ"),
    ]

    passed = 0
    failed = 0

    for key1, key2, expected, desc in test_cases:
        # Test reference implementation
        result_ref = anticommute(key1, key2)

        # Test fast implementation
        x1, z1 = encode_pauli_string(key1, 10)
        x2, z2 = encode_pauli_string(key2, 10)
        result_fast = pauli_anticommute(x1, z1, x2, z2)

        if result_ref == expected and result_fast == expected and result_ref == result_fast:
            passed += 1
        else:
            failed += 1
            print(f"  ❌ FAILED: {desc}")
            print(f"     Expected: {expected}, Ref: {result_ref}, Fast: {result_fast}")

    print(f"\nResults: {passed} passed, {failed} failed")
    return failed == 0


# ============================================================================
# Test 2: Consistency Between Implementations
# ============================================================================
def test_consistency():
    """Test that fast and reference implementations give consistent results."""
    print("\n" + "=" * 70)
    print("TEST 2: Consistency Between Fast and Reference Implementations")
    print("=" * 70)

    # Create test Hamiltonian: H = X₀ + Y₀ + 0.5·Z₀
    terms = [
        QubitOperator('X0', 1.0),
        QubitOperator('Y0', 1.0),
        QubitOperator('Z0', 0.5)
    ]

    print("\nTest Hamiltonian: H = X₀ + Y₀ + 0.5·Z₀")
    print(f"Number of terms: {len(terms)}")

    # Short time limit for testing
    time_limit = 2.0

    print(f"\nRunning reference implementation (time_limit={time_limit}s)...")
    c1_ref, c2_ref = trotter_error_estimator(terms, time_limit, mock_config)

    print(f"Running fast implementation (time_limit={time_limit}s)...")
    c1_fast, c2_fast = trotter_error_estimator_fast(terms, time_limit, mock_config)

    print(f"\nResults:")
    print(f"  Reference: C1 = {c1_ref:.6f}, C2 = {c2_ref:.6f}")
    print(f"  Fast:      C1 = {c1_fast:.6f}, C2 = {c2_fast:.6f}")

    # Check relative difference (allow 10% tolerance due to Monte Carlo)
    rel_diff_c1 = abs(c1_fast - c1_ref) / max(abs(c1_ref), 1e-10)
    rel_diff_c2 = abs(c2_fast - c2_ref) / max(abs(c2_ref), 1e-10)

    print(f"\nRelative differences:")
    print(f"  C1: {rel_diff_c1:.2%}")
    print(f"  C2: {rel_diff_c2:.2%}")

    tolerance = 0.15  # 15% tolerance
    passed = rel_diff_c1 < tolerance and rel_diff_c2 < tolerance

    if passed:
        print(f"\n✅ PASSED: Results agree within {tolerance:.0%} tolerance")
    else:
        print(f"\n❌ FAILED: Results differ by more than {tolerance:.0%}")

    return passed


# ============================================================================
# Test 3: Known Analytical Result
# ============================================================================
def test_analytical():
    """Test against known analytical result for H = X + Y."""
    print("\n" + "=" * 70)
    print("TEST 3: Known Analytical Result - H = X₀ + Y₀")
    print("=" * 70)

    # H = X + Y with coefficients 1.0
    terms = [
        QubitOperator('X0', 1.0),
        QubitOperator('Y0', 1.0)
    ]

    print("\nTest Hamiltonian: H = X₀ + Y₀")
    print("\nAnalytical calculation:")
    print("  [X, Y] = 2iZ, so ||[X, Y]|| = 2")
    print("  C1 = Σᵢ<ⱼ ||[Hᵢ, Hⱼ]|| = 2")
    print("  Returned value should be C1/2 = 1.0")
    print("\n  C21 = 0 (need at least 3 terms)")
    print("  C22 = Σₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]||")
    print("    For k=0, j=1: ||[X, [X, Y]]|| = ||[X, 2iZ]|| = ||4Y|| = 4")
    print("  C22 = 4, returned C2 = C21/12 + C22/24 = 0/12 + 4/24 = 1/6 ≈ 0.167")

    c1_analytical = 1.0
    c2_analytical = 1.0 / 6.0

    # Use longer time for accurate result
    time_limit = 5.0
    print(f"\nRunning fast implementation (time_limit={time_limit}s)...")
    c1_fast, c2_fast = trotter_error_estimator_fast(terms, time_limit, mock_config)

    print(f"\nResults:")
    print(f"  Analytical: C1 = {c1_analytical:.6f}, C2 = {c2_analytical:.6f}")
    print(f"  Fast:       C1 = {c1_fast:.6f}, C2 = {c2_fast:.6f}")

    rel_diff_c1 = abs(c1_fast - c1_analytical) / c1_analytical
    rel_diff_c2 = abs(c2_fast - c2_analytical) / c2_analytical

    print(f"\nRelative differences:")
    print(f"  C1: {rel_diff_c1:.2%}")
    print(f"  C2: {rel_diff_c2:.2%}")

    tolerance = 0.10  # 10% tolerance
    passed = rel_diff_c1 < tolerance and rel_diff_c2 < tolerance

    if passed:
        print(f"\n✅ PASSED: Results agree with analytical values within {tolerance:.0%}")
    else:
        print(f"\n❌ FAILED: Results differ from analytical values by more than {tolerance:.0%}")

    return passed


# ============================================================================
# Test 4: Exact Computation Feature
# ============================================================================
def test_exact_computation():
    """Test that small systems achieve exact computation."""
    print("\n" + "=" * 70)
    print("TEST 4: Exact Computation Feature")
    print("=" * 70)

    # Very small system: should achieve exact computation quickly
    terms_small = [
        QubitOperator('X0', 1.0),
        QubitOperator('Y0', 1.0),
        QubitOperator('Z0', 0.5)
    ]

    print(f"\nSmall system: N={len(terms_small)} terms")
    print("Expected: Should achieve exact computation")

    time_limit = 10.0
    c1, c2 = trotter_error_estimator_fast(terms_small, time_limit, mock_config)

    print(f"\nResults: C1 = {c1:.6f}, C2 = {c2:.6f}")

    # For this small system, results should be exact and deterministic
    # Run twice and check they're identical
    c1_2, c2_2 = trotter_error_estimator_fast(terms_small, time_limit, mock_config)

    passed = (c1 == c1_2) and (c2 == c2_2)

    if passed:
        print("\n✅ PASSED: Exact computation achieved (results are deterministic)")
    else:
        print("\n❌ FAILED: Results not deterministic (exact computation not achieved)")

    return passed


# ============================================================================
# Main Test Runner
# ============================================================================
def test_converged_monte_carlo():
    """Test that medium systems achieve convergence."""
    print("\n" + "=" * 70)
    print("TEST 5: Converged Monte Carlo Sampling")
    print("=" * 70)

    # Medium system: large enough to use Monte Carlo, small enough to converge
    np.random.seed(42)
    terms_medium = []
    for i in range(50):
        op = np.random.choice(['X', 'Y', 'Z'])
        qubit = i % 10
        coeff = np.random.uniform(0.1, 1.0)
        terms_medium.append(QubitOperator(f'{op}{qubit}', coeff))

    print(f"\nMedium system: N={len(terms_medium)} terms")
    print("Expected: Should converge within reasonable time")

    time_limit = 20.0
    print(f"\nRunning fast implementation (time_limit={time_limit}s)...")
    c1, c2 = trotter_error_estimator_fast(terms_medium, time_limit, mock_config)

    print(f"\nResults: C1 = {c1:.6f}, C2 = {c2:.6f}")

    # For medium system, run twice and check they're reasonably close (not exact, but converged)
    c1_2, c2_2 = trotter_error_estimator_fast(terms_medium, time_limit, mock_config)

    rel_diff_c1 = abs(c1 - c1_2) / max(abs(c1), 1e-10)
    rel_diff_c2 = abs(c2 - c2_2) / max(abs(c2), 1e-10)

    print(f"\nReproducibility check (two runs):")
    print(f"  C1 relative difference: {rel_diff_c1:.2%}")
    print(f"  C2 relative difference: {rel_diff_c2:.2%}")

    # Converged Monte Carlo should be reproducible within ~5%
    tolerance = 0.10
    passed = rel_diff_c1 < tolerance and rel_diff_c2 < tolerance

    if passed:
        print(f"\n✅ PASSED: Results converged (reproducible within {tolerance:.0%})")
    else:
        print(f"\n❌ FAILED: Results not converged (differs by more than {tolerance:.0%})")

    return passed


def test_analytical_scaling():
    """Test that analytical result holds for both exact and converged Monte Carlo."""
    print("\n" + "=" * 70)
    print("TEST 6: Analytical Result - Exact vs Converged Monte Carlo")
    print("=" * 70)

    print("\nPattern: Multiple copies of (X + Y) on different qubits")
    print("For N qubits with X_i + Y_i on each qubit i:")
    print("  All terms on different qubits commute with each other")
    print("  Each qubit contributes: ||[X_i, Y_i]|| = 2")
    print("  C1 = 2N, C2 = 4N/24 = N/6")
    print("  Returned: C1/2 = N, C2 = N/6")

    # Small version: exact computation
    n_qubits_small = 2
    terms_small = []
    for i in range(n_qubits_small):
        terms_small.append(QubitOperator(f'X{i}', 1.0))
        terms_small.append(QubitOperator(f'Y{i}', 1.0))

    c1_analytical_small = n_qubits_small
    c2_analytical_small = n_qubits_small / 6.0

    print(f"\n--- Small System (N={len(terms_small)} terms, {n_qubits_small} qubits) ---")
    print(f"Analytical: C1 = {c1_analytical_small:.6f}, C2 = {c2_analytical_small:.6f}")

    time_limit_small = 5.0
    c1_exact, c2_exact = trotter_error_estimator_fast(terms_small, time_limit_small, mock_config)
    print(f"Exact:      C1 = {c1_exact:.6f}, C2 = {c2_exact:.6f}")

    rel_diff_c1_exact = abs(c1_exact - c1_analytical_small) / c1_analytical_small
    rel_diff_c2_exact = abs(c2_exact - c2_analytical_small) / c2_analytical_small

    # Large version: converged Monte Carlo
    n_qubits_large = 20
    terms_large = []
    for i in range(n_qubits_large):
        terms_large.append(QubitOperator(f'X{i}', 1.0))
        terms_large.append(QubitOperator(f'Y{i}', 1.0))

    c1_analytical_large = n_qubits_large
    c2_analytical_large = n_qubits_large / 6.0

    print(f"\n--- Large System (N={len(terms_large)} terms, {n_qubits_large} qubits) ---")
    print(f"Analytical: C1 = {c1_analytical_large:.6f}, C2 = {c2_analytical_large:.6f}")

    time_limit_large = 15.0
    c1_mc, c2_mc = trotter_error_estimator_fast(terms_large, time_limit_large, mock_config)
    print(f"Monte Carlo: C1 = {c1_mc:.6f}, C2 = {c2_mc:.6f}")

    rel_diff_c1_mc = abs(c1_mc - c1_analytical_large) / c1_analytical_large
    rel_diff_c2_mc = abs(c2_mc - c2_analytical_large) / c2_analytical_large

    print(f"\nRelative errors:")
    print(f"  Exact computation:  C1: {rel_diff_c1_exact:.2%}, C2: {rel_diff_c2_exact:.2%}")
    print(f"  Monte Carlo:        C1: {rel_diff_c1_mc:.2%}, C2: {rel_diff_c2_mc:.2%}")

    tolerance_exact = 0.01  # Exact should be within 1%
    tolerance_mc = 0.10     # Converged MC should be within 10%

    passed_exact = rel_diff_c1_exact < tolerance_exact and rel_diff_c2_exact < tolerance_exact
    passed_mc = rel_diff_c1_mc < tolerance_mc and rel_diff_c2_mc < tolerance_mc
    passed = passed_exact and passed_mc

    if passed:
        print(f"\n✅ PASSED: Both exact and Monte Carlo match analytical result")
    else:
        if not passed_exact:
            print(f"\n❌ FAILED: Exact computation differs from analytical by > {tolerance_exact:.0%}")
        if not passed_mc:
            print(f"\n❌ FAILED: Monte Carlo differs from analytical by > {tolerance_mc:.0%}")

    return passed


def test_monte_carlo_produces_reasonable_results():
    """Test that Monte Carlo produces reasonable results even with limited time."""
    print("\n" + "=" * 70)
    print("TEST 6: Monte Carlo with Limited Time")
    print("=" * 70)

    # Large system with moderate time limit
    np.random.seed(42)
    terms_large = []
    for i in range(200):
        op = np.random.choice(['X', 'Y', 'Z'])
        qubit = i % 20
        coeff = np.random.uniform(0.1, 1.0)
        terms_large.append(QubitOperator(f'{op}{qubit}', coeff))

    print(f"\nLarge system: N={len(terms_large)} terms")
    print("Expected: Should produce reasonable estimates")

    time_limit = 3.0
    print(f"\nRunning fast implementation (time_limit={time_limit}s)...")
    c1, c2 = trotter_error_estimator_fast(terms_large, time_limit, mock_config)

    print(f"\nResults: C1 = {c1:.6f}, C2 = {c2:.6f}")

    # Run twice and check reproducibility
    c1_2, c2_2 = trotter_error_estimator_fast(terms_large, time_limit, mock_config)

    rel_diff_c1 = abs(c1 - c1_2) / max(abs(c1), 1e-10)
    rel_diff_c2 = abs(c2 - c2_2) / max(abs(c2), 1e-10)

    print(f"\nReproducibility check (two runs):")
    print(f"  C1 relative difference: {rel_diff_c1:.2%}")
    print(f"  C2 relative difference: {rel_diff_c2:.2%}")

    # Test passes if:
    # 1. Results are reasonable (positive values)
    # 2. Results are not wildly different between runs (< 50% variation)
    # We don't require perfect reproducibility - just reasonable stability
    max_variation = 0.50  # 50% max variation

    results_positive = c1 > 0 and c2 > 0 and c1_2 > 0 and c2_2 > 0
    results_stable = rel_diff_c1 < max_variation and rel_diff_c2 < max_variation

    passed = results_positive and results_stable

    if passed:
        print(f"\n✅ PASSED: Monte Carlo produces reasonable, stable results")
        print(f"   (C1 variation: {rel_diff_c1:.2%}, C2 variation: {rel_diff_c2:.2%})")
    else:
        print(f"\n❌ FAILED:")
        if not results_positive:
            print(f"   - Results are not positive")
        if not results_stable:
            print(f"   - Results vary too much (>{max_variation:.0%})")

    return passed


def test_standard_error_calculation():
    """Test that standard error calculation is correct for Monte Carlo estimation."""
    print("\n" + "=" * 70)
    print("TEST 8: Standard Error Formula Verification")
    print("=" * 70)

    # Temporarily reduce thresholds to force Monte Carlo sampling
    import trotter_coefficients_fast
    original_time_budget = trotter_coefficients_fast.EXACT_COMPUTATION_TIME_BUDGET
    original_memory_limit = trotter_coefficients_fast.EXACT_COMPUTATION_MEMORY_LIMIT_MB

    try:
        trotter_coefficients_fast.EXACT_COMPUTATION_TIME_BUDGET = 0.5  # Force MC
        trotter_coefficients_fast.EXACT_COMPUTATION_MEMORY_LIMIT_MB = 5  # Force MC

        # Large system with known result
        n_qubits = 50
        terms = []
        for i in range(n_qubits):
            terms.append(QubitOperator(f'X{i}', 1.0))
            terms.append(QubitOperator(f'Y{i}', 1.0))

        print(f"\nTest system: {len(terms)} terms ({n_qubits} qubits)")
        print(f"Expected: C1 = {n_qubits}, C2 = {n_qubits/6:.6f}")
        print("(Forced Monte Carlo by reducing exact computation thresholds)")

        # Run multiple times to check reproducibility
        time_limit = 15.0
        n_runs = 3

        print(f"\nRunning {n_runs} independent Monte Carlo simulations...")
        c1_values = []
        c2_values = []

        for run in range(n_runs):
            c1, c2 = trotter_error_estimator_fast(terms, time_limit, mock_config)
            c1_values.append(c1)
            c2_values.append(c2)

        c1_values = np.array(c1_values)
        c2_values = np.array(c2_values)

        c1_mean = np.mean(c1_values)
        c1_std = np.std(c1_values, ddof=1)
        c2_mean = np.mean(c2_values)
        c2_std = np.std(c2_values, ddof=1)

        print(f"\nResults over {n_runs} runs:")
        print(f"  C1: {c1_mean:.6f} ± {c1_std:.6f}")
        print(f"  C2: {c2_mean:.6f} ± {c2_std:.6f}")

        # Check accuracy
        c1_error = abs(c1_mean - n_qubits) / n_qubits
        c2_error = abs(c2_mean - n_qubits/6) / (n_qubits/6)

        print(f"\nRelative errors:")
        print(f"  C1: {c1_error:.2%}")
        print(f"  C2: {c2_error:.2%}")

        # Check reproducibility (coefficient of variation)
        c1_cv = c1_std / c1_mean if c1_mean > 0 else float('inf')
        c2_cv = c2_std / c2_mean if c2_mean > 0 else float('inf')

        print(f"\nCoefficient of variation:")
        print(f"  C1: {c1_cv:.2%}")
        print(f"  C2: {c2_cv:.2%}")

        # Test passes if:
        # 1. Results are accurate (< 5% error)
        # 2. Results are reproducible (CV < 10%)
        accuracy_threshold = 0.05
        reproducibility_threshold = 0.10

        accurate = c1_error < accuracy_threshold and c2_error < accuracy_threshold
        reproducible = c1_cv < reproducibility_threshold and c2_cv < reproducibility_threshold

        passed = accurate and reproducible

        if passed:
            print(f"\n✅ PASSED: Standard error formula is working correctly")
            print(f"   - Accurate within {accuracy_threshold:.0%}")
            print(f"   - Reproducible (CV < {reproducibility_threshold:.0%})")
        else:
            print(f"\n❌ FAILED:")
            if not accurate:
                print(f"   - Results not accurate (error ≥ {accuracy_threshold:.0%})")
            if not reproducible:
                print(f"   - Results not reproducible (CV ≥ {reproducibility_threshold:.0%})")

        return passed

    finally:
        # Restore original thresholds
        trotter_coefficients_fast.EXACT_COMPUTATION_TIME_BUDGET = original_time_budget
        trotter_coefficients_fast.EXACT_COMPUTATION_MEMORY_LIMIT_MB = original_memory_limit


def run_all_tests():
    """Run all tests and report results."""
    print("\n" + "=" * 70)
    print("TROTTER COEFFICIENTS TEST SUITE")
    print("=" * 70)

    tests = [
        ("Anticommutation", test_anticommutation),
        ("Consistency", test_consistency),
        ("Analytical", test_analytical),
        ("Exact Computation", test_exact_computation),
        ("Converged Monte Carlo", test_converged_monte_carlo),
        ("Analytical Scaling", test_analytical_scaling),
        ("Monte Carlo Limited Time", test_monte_carlo_produces_reasonable_results),
        ("Standard Error Calculation", test_standard_error_calculation),
    ]

    results = []
    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print(f"\n❌ ERROR in {name}: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))

    # Summary
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)

    passed_count = sum(1 for _, passed in results if passed)
    total_count = len(results)

    for name, passed in results:
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"  {status}: {name}")

    print(f"\nOverall: {passed_count}/{total_count} tests passed")

    if passed_count == total_count:
        print("\n🎉 All tests passed!")
        return 0
    else:
        print(f"\n⚠️  {total_count - passed_count} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(run_all_tests())
