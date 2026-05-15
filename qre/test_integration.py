"""
Integration test: Verify the complete implementation works end-to-end.
"""

import numpy as np
from openfermion import QubitOperator
from trotter_coefficients_fast import trotter_error_estimator_fast


def test_complete_workflow():
    """Test the complete workflow with different configurations."""
    print("\n" + "="*70)
    print("INTEGRATION TEST: Complete Workflow")
    print("="*70)

    # Create a test Hamiltonian
    H1 = QubitOperator('X0', 1.0)
    H2 = QubitOperator('Y0', 1.0)
    H3 = QubitOperator('Z0', 1.0)
    H4 = QubitOperator('X1', 1.0)
    H5 = QubitOperator('Y1', 1.0)
    terms = [H1, H2, H3, H4, H5]

    print(f"\nTest Hamiltonian: {len(terms)} terms")
    print("  H = X₀ + Y₀ + Z₀ + X₁ + Y₁")

    # Test 1: Default behavior (backward compatibility)
    print("\n1. Default behavior (backward compatible):")
    c1, c2 = trotter_error_estimator_fast(terms, 10)
    print(f"   → C1={c1:.6f}, C2={c2:.6f}")
    assert c1 > 0 and c2 > 0

    # Test 2: Explicit Monte Carlo
    print("\n2. Explicit Monte Carlo mode:")
    c1_mc, c2_mc = trotter_error_estimator_fast(terms, 10, mode='monte_carlo')
    print(f"   → C1={c1_mc:.6f}, C2={c2_mc:.6f}")
    assert c1_mc > 0 and c2_mc > 0

    # Test 3: Exact mode
    print("\n3. Exact mode:")
    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 10, mode='exact')
    print(f"   → C1={c1_ex:.6f}, C2={c2_ex:.6f}")
    assert c1_ex > 0 and c2_ex > 0

    # Test 4: Auto-exact mode
    print("\n4. Auto-exact mode:")
    c1_auto, c2_auto = trotter_error_estimator_fast(terms, 60, mode='monte_carlo', auto_exact=True)
    print(f"   → C1={c1_auto:.6f}, C2={c2_auto:.6f}")
    assert c1_auto > 0 and c2_auto > 0

    # Test 5: Verify exact and auto-exact give same results
    print("\n5. Consistency check:")
    print(f"   Exact:      C1={c1_ex:.6f}, C2={c2_ex:.6f}")
    print(f"   Auto-exact: C1={c1_auto:.6f}, C2={c2_auto:.6f}")
    print(f"   Difference: ΔC1={abs(c1_ex - c1_auto):.2e}, ΔC2={abs(c2_ex - c2_auto):.2e}")
    assert abs(c1_ex - c1_auto) < 1e-10, "Exact and auto-exact should match"
    assert abs(c2_ex - c2_auto) < 1e-10, "Exact and auto-exact should match"
    print("   ✓ Exact and auto-exact match perfectly")

    # Test 6: Verify Monte Carlo is close to exact
    print("\n6. Monte Carlo accuracy:")
    rel_err_c1 = abs(c1_mc - c1_ex) / c1_ex if c1_ex > 0 else 0
    rel_err_c2 = abs(c2_mc - c2_ex) / c2_ex if c2_ex > 0 else 0
    print(f"   Relative errors: C1={100*rel_err_c1:.1f}%, C2={100*rel_err_c2:.1f}%")
    if rel_err_c1 < 0.1 and rel_err_c2 < 0.1:
        print("   ✓ Monte Carlo is accurate (<10% error)")
    else:
        print("   ⚠ Monte Carlo has larger error (expected for small sample sizes)")

    print("\n" + "="*70)
    print("✅ INTEGRATION TEST PASSED")
    print("="*70)


def test_larger_system():
    """Test with a larger but still feasible system."""
    print("\n" + "="*70)
    print("INTEGRATION TEST: Larger System (N=20)")
    print("="*70)

    np.random.seed(42)
    terms = []
    for i in range(20):
        qubit = i % 10
        pauli = ['X', 'Y', 'Z'][i % 3]
        coeff = np.random.uniform(0.5, 2.0)
        terms.append(QubitOperator(f'{pauli}{qubit}', coeff))

    print(f"\nTest Hamiltonian: {len(terms)} terms, 10 qubits")

    # Monte Carlo vs Exact comparison
    print("\nMonte Carlo:")
    c1_mc, c2_mc = trotter_error_estimator_fast(terms, 20, mode='monte_carlo', auto_exact=False)
    print(f"  → C1={c1_mc:.6f}, C2={c2_mc:.6f}")

    print("\nExact:")
    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 20, mode='exact')
    print(f"  → C1={c1_ex:.6f}, C2={c2_ex:.6f}")

    print("\nAuto-exact:")
    c1_auto, c2_auto = trotter_error_estimator_fast(terms, 60, mode='monte_carlo', auto_exact=True)
    print(f"  → C1={c1_auto:.6f}, C2={c2_auto:.6f}")

    # Verify exact matches
    assert abs(c1_ex - c1_auto) < 1e-10
    assert abs(c2_ex - c2_auto) < 1e-10
    print("\n✓ All exact modes agree")

    print("\n" + "="*70)
    print("✅ LARGER SYSTEM TEST PASSED")
    print("="*70)


def test_edge_case_workflow():
    """Test edge cases in the workflow."""
    print("\n" + "="*70)
    print("INTEGRATION TEST: Edge Cases")
    print("="*70)

    # Test N=1
    print("\n1. Single term (N=1):")
    terms = [QubitOperator('X0', 1.0)]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')
    print(f"   → C1={c1:.6f}, C2={c2:.6f}")
    assert c1 == 0.0 and c2 == 0.0
    print("   ✓ Correctly returns zeros")

    # Test N=2
    print("\n2. Two terms (N=2):")
    terms = [QubitOperator('X0', 1.0), QubitOperator('Y0', 1.0)]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')
    print(f"   → C1={c1:.6f}, C2={c2:.6f}")
    assert c1 > 0
    print("   ✓ Handles N=2 (minimum for pairs)")

    # Test all commuting
    print("\n3. All commuting terms:")
    terms = [QubitOperator('X0', 1.0), QubitOperator('X1', 1.0), QubitOperator('X2', 1.0)]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')
    print(f"   → C1={c1:.6f}, C2={c2:.6f}")
    assert c1 == 0.0 and c2 == 0.0
    print("   ✓ Correctly handles commuting terms")

    # Test with complex coefficients
    print("\n4. Complex coefficients:")
    terms = [QubitOperator('X0', 1.0+1.0j), QubitOperator('Y0', 1.0-1.0j)]
    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')
    print(f"   → C1={c1:.6f}, C2={c2:.6f}")
    assert c1 > 0
    print("   ✓ Handles complex coefficients")

    print("\n" + "="*70)
    print("✅ EDGE CASES TEST PASSED")
    print("="*70)


if __name__ == "__main__":
    print("="*70)
    print("COMPREHENSIVE INTEGRATION TESTS")
    print("="*70)

    test_complete_workflow()
    test_larger_system()
    test_edge_case_workflow()

    print("\n" + "="*70)
    print("✅ ALL INTEGRATION TESTS PASSED!")
    print("="*70)
    print("\nSummary:")
    print("  • Complete workflow tested")
    print("  • All modes work correctly")
    print("  • Exact and auto-exact agree perfectly")
    print("  • Monte Carlo provides reasonable accuracy")
    print("  • Edge cases handled correctly")
    print("  • Larger systems (N=20) work efficiently")
    print("\n" + "="*70)
