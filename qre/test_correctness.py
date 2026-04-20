"""
Test script to verify correctness of both implementations.

This script tests:
1. Ground truth: Both implementations against known mathematical properties
2. Consistency: The optimized implementation matches the original
"""

import numpy as np
from openfermion import QubitOperator

# Import original functions
from trotter_coefficients import anticommute, pauli_product_key, compute_commutator

# Import optimized functions
from trotter_coefficients_fast import (
    encode_pauli_string,
    pauli_anticommute,
    pauli_product_bits,
    compute_commutator_norm_fast
)


def test_ground_truth_anticommutation():
    """Test that both implementations correctly identify Pauli anticommutation."""
    print("\n" + "=" * 70)
    print("TEST: Ground Truth - Pauli Anticommutation Properties")
    print("=" * 70)

    # Known mathematical facts about Pauli matrices:
    # - Different Paulis on same qubit anticommute
    # - Same Paulis on same qubit commute
    # - Paulis on different qubits commute
    test_cases = [
        # (key1, key2, should_anticommute, description)
        (((0, 'X'),), ((0, 'Y'),), True, "X and Y on same qubit"),
        (((0, 'X'),), ((0, 'Z'),), True, "X and Z on same qubit"),
        (((0, 'Y'),), ((0, 'Z'),), True, "Y and Z on same qubit"),
        (((0, 'X'),), ((0, 'X'),), False, "X and X on same qubit"),
        (((0, 'Y'),), ((0, 'Y'),), False, "Y and Y on same qubit"),
        (((0, 'Z'),), ((0, 'Z'),), False, "Z and Z on same qubit"),
        (((0, 'X'),), ((1, 'X'),), False, "X on different qubits"),
        (((0, 'X'),), ((1, 'Y'),), False, "X and Y on different qubits"),
        # Multi-qubit: odd number of anticommuting positions = anticommute
        (((0, 'X'), (1, 'Y')), ((0, 'Y'), (1, 'Y')), True, "1 position differs"),
        (((0, 'X'), (1, 'Y')), ((0, 'Y'), (1, 'X')), False, "2 positions differ"),
        (((0, 'X'), (1, 'Y'), (2, 'Z')), ((0, 'Y'), (1, 'X'), (2, 'Z')), False, "2 positions differ"),
        (((0, 'X'), (1, 'Y'), (2, 'Z')), ((0, 'Y'), (1, 'X'), (2, 'X')), True, "3 positions differ (X0Y1Z2 vs Y0X1X2)"),
    ]

    all_passed = True
    for key1, key2, expected, description in test_cases:
        # Test original implementation
        result_orig = anticommute(key1, key2)

        # Test fast implementation
        x1, z1 = encode_pauli_string(key1, 10)
        x2, z2 = encode_pauli_string(key2, 10)
        result_fast = pauli_anticommute(x1, z1, x2, z2)

        # Both should match expected
        passed = (result_orig == expected) and (result_fast == expected)

        status = "✓" if passed else "❌"
        print(f"  {status} {description}")
        if not passed:
            print(f"     Expected: {expected}, Original: {result_orig}, Fast: {result_fast}")
            all_passed = False

    if all_passed:
        print("\n  ✅ All ground truth anticommutation tests PASSED")
    else:
        print("\n  ❌ Some ground truth anticommutation tests FAILED")

    return all_passed


def test_ground_truth_commutator_norms():
    """Test that both implementations compute correct commutator norms."""
    print("\n" + "=" * 70)
    print("TEST: Ground Truth - Commutator Norms")
    print("=" * 70)

    # [A, B] = AB - BA
    # If A and B commute: [A, B] = 0, norm = 0
    # If A and B anticommute: AB = -BA, so [A, B] = 2AB, norm = 2|coeff_A * coeff_B|
    test_cases = [
        # (op1, op2, expected_norm, description)
        (QubitOperator('X0', 1.0), QubitOperator('Y0', 1.0), 2.0, "[X, Y] with unit coefficients"),
        (QubitOperator('X0', 2.0), QubitOperator('Y0', 3.0), 12.0, "[X, Y] with coefficients 2 and 3"),
        (QubitOperator('X0', 1.0), QubitOperator('X0', 1.0), 0.0, "[X, X] = 0"),
        (QubitOperator('X0', 5.0), QubitOperator('X0', 7.0), 0.0, "[X, X] = 0 regardless of coefficients"),
        (QubitOperator('X0', 1.0), QubitOperator('X1', 1.0), 0.0, "Different qubits commute"),
        (QubitOperator('X0 Y1', 1.0), QubitOperator('Y0 Y1', 1.0), 2.0, "Anticommute at 1 position"),
        (QubitOperator('X0 Y1', 2.0), QubitOperator('Y0 X1', 3.0), 0.0, "Commute (2 position differences)"),
        # Complex coefficients
        (QubitOperator('X0', 1.0j), QubitOperator('Y0', 1.0), 2.0, "[X, Y] with imaginary coefficient"),
        (QubitOperator('Z0', 1+1j), QubitOperator('X0', 1-1j), 2*np.sqrt(2)*np.sqrt(2), "[Z, X] with complex coefficients"),
    ]

    all_passed = True
    for op1, op2, expected_norm, description in test_cases:
        # Test original implementation
        _, norm_orig = compute_commutator(op1, op2)

        # Test fast implementation
        key1 = list(op1.terms.keys())[0]
        key2 = list(op2.terms.keys())[0]
        c1 = list(op1.terms.values())[0]
        c2 = list(op2.terms.values())[0]

        x1, z1 = encode_pauli_string(key1, 10)
        x2, z2 = encode_pauli_string(key2, 10)
        norm_fast = compute_commutator_norm_fast(x1, z1, c1, x2, z2, c2)

        # Check both match expected
        diff_orig = abs(norm_orig - expected_norm)
        diff_fast = abs(norm_fast - expected_norm)
        passed = (diff_orig < 1e-10) and (diff_fast < 1e-10)

        status = "✓" if passed else "❌"
        print(f"  {status} {description}")
        if not passed:
            print(f"     Expected: {expected_norm:.6f}")
            print(f"     Original: {norm_orig:.6f} (diff: {diff_orig:.2e})")
            print(f"     Fast:     {norm_fast:.6f} (diff: {diff_fast:.2e})")
            all_passed = False

    if all_passed:
        print("\n  ✅ All ground truth commutator norm tests PASSED")
    else:
        print("\n  ❌ Some ground truth commutator norm tests FAILED")

    return all_passed


def test_ground_truth_pauli_products():
    """Test that both implementations correctly compute Pauli products."""
    print("\n" + "=" * 70)
    print("TEST: Ground Truth - Pauli Product Rules")
    print("=" * 70)

    # Pauli product rules (ignoring phase factors):
    # X*X = I, Y*Y = I, Z*Z = I
    # X*Y = Z, Y*Z = X, Z*X = Y (cyclic)
    # Y*X = Z, Z*Y = X, X*Z = Y (reverse)
    test_cases = [
        # (key1, key2, expected_product_key, description)
        (((0, 'X'),), ((0, 'X'),), (), "X*X = I"),
        (((0, 'Y'),), ((0, 'Y'),), (), "Y*Y = I"),
        (((0, 'Z'),), ((0, 'Z'),), (), "Z*Z = I"),
        (((0, 'X'),), ((0, 'Y'),), ((0, 'Z'),), "X*Y = Z"),
        (((0, 'Y'),), ((0, 'Z'),), ((0, 'X'),), "Y*Z = X"),
        (((0, 'Z'),), ((0, 'X'),), ((0, 'Y'),), "Z*X = Y"),
        (((0, 'Y'),), ((0, 'X'),), ((0, 'Z'),), "Y*X = Z"),
        (((0, 'Z'),), ((0, 'Y'),), ((0, 'X'),), "Z*Y = X"),
        (((0, 'X'),), ((0, 'Z'),), ((0, 'Y'),), "X*Z = Y"),
        # Different qubits: product keeps both
        (((0, 'X'),), ((1, 'Y'),), ((0, 'X'), (1, 'Y')), "X0 * Y1 = X0 Y1"),
        # Multi-qubit products
        (((0, 'X'), (1, 'Y')), ((0, 'Y'), (1, 'Y')), ((0, 'Z'),), "X0Y1 * Y0Y1 = Z0"),
    ]

    all_passed = True
    for key1, key2, expected_product, description in test_cases:
        # Test original implementation
        product_orig = pauli_product_key(key1, key2)

        # Test fast implementation
        x1, z1 = encode_pauli_string(key1, 10)
        x2, z2 = encode_pauli_string(key2, 10)
        x_prod, z_prod = pauli_product_bits(x1, z1, x2, z2)

        # Convert fast result back to key format
        product_fast = []
        for qubit in range(10):
            x_bit = (x_prod >> qubit) & 1
            z_bit = (z_prod >> qubit) & 1
            if x_bit or z_bit:
                if x_bit and z_bit:
                    product_fast.append((qubit, 'Y'))
                elif x_bit:
                    product_fast.append((qubit, 'X'))
                elif z_bit:
                    product_fast.append((qubit, 'Z'))
        product_fast = tuple(product_fast)

        # Both should match expected
        passed = (product_orig == expected_product) and (product_fast == expected_product)

        status = "✓" if passed else "❌"
        print(f"  {status} {description}")
        if not passed:
            print(f"     Expected: {expected_product}")
            print(f"     Original: {product_orig}")
            print(f"     Fast:     {product_fast}")
            all_passed = False

    if all_passed:
        print("\n  ✅ All ground truth Pauli product tests PASSED")
    else:
        print("\n  ❌ Some ground truth Pauli product tests FAILED")

    return all_passed


def test_encoding():
    """Test that Pauli string encoding is correct."""
    print("\n" + "=" * 70)
    print("TEST 1: Pauli String Encoding")
    print("=" * 70)

    test_cases = [
        ((), "Identity"),
        (((0, 'X'),), "X on qubit 0"),
        (((1, 'Y'),), "Y on qubit 1"),
        (((2, 'Z'),), "Z on qubit 2"),
        (((0, 'X'), (1, 'Y'), (2, 'Z')), "X0 Y1 Z2"),
        (((0, 'X'), (3, 'Y'), (7, 'Z')), "X0 Y3 Z7"),
    ]

    all_passed = True
    for key, description in test_cases:
        x_bits, z_bits = encode_pauli_string(key, 10)

        # Verify encoding
        print(f"\n  {description}: {key}")
        print(f"    x_bits: {bin(x_bits)}")
        print(f"    z_bits: {bin(z_bits)}")

        # Check that encoding is correct
        for qubit, op in key:
            if op == 'I':
                expected_x = 0
                expected_z = 0
            elif op == 'X':
                expected_x = 1
                expected_z = 0
            elif op == 'Y':
                expected_x = 1
                expected_z = 1
            elif op == 'Z':
                expected_x = 0
                expected_z = 1

            actual_x = (x_bits >> qubit) & 1
            actual_z = (z_bits >> qubit) & 1

            if actual_x != expected_x or actual_z != expected_z:
                print(f"    ❌ FAILED for qubit {qubit}")
                all_passed = False
            else:
                print(f"    ✓ Qubit {qubit}: correct")

    if all_passed:
        print("\n  ✅ All encoding tests PASSED")
    else:
        print("\n  ❌ Some encoding tests FAILED")

    return all_passed


def test_anticommutation():
    """Test that anticommutation check matches original."""
    print("\n" + "=" * 70)
    print("TEST 2: Anticommutation Check")
    print("=" * 70)

    test_pairs = [
        # (key1, key2, expected_anticommute)
        (((0, 'X'),), ((0, 'Y'),), True),   # X and Y anticommute
        (((0, 'X'),), ((0, 'Z'),), True),   # X and Z anticommute
        (((0, 'Y'),), ((0, 'Z'),), True),   # Y and Z anticommute
        (((0, 'X'),), ((0, 'X'),), False),  # X and X commute
        (((0, 'X'),), ((1, 'X'),), False),  # Different qubits commute
        (((0, 'X'), (1, 'Y')), ((0, 'Y'), (1, 'X')), False),  # Even number of differences
        (((0, 'X'), (1, 'Y')), ((0, 'Y'), (1, 'Y')), True),   # Odd number of differences
    ]

    all_passed = True
    for key1, key2, expected in test_pairs:
        # Original implementation
        result_orig = anticommute(key1, key2)

        # Fast implementation
        x1, z1 = encode_pauli_string(key1, 10)
        x2, z2 = encode_pauli_string(key2, 10)
        result_fast = pauli_anticommute(x1, z1, x2, z2)

        passed = (result_orig == result_fast == expected)

        status = "✓" if passed else "❌"
        print(f"\n  {status} {key1} vs {key2}")
        print(f"     Original: {result_orig}, Fast: {result_fast}, Expected: {expected}")

        if not passed:
            all_passed = False

    if all_passed:
        print("\n  ✅ All anticommutation tests PASSED")
    else:
        print("\n  ❌ Some anticommutation tests FAILED")

    return all_passed


def test_commutator_norm():
    """Test that commutator norm calculation matches original."""
    print("\n" + "=" * 70)
    print("TEST 3: Commutator Norm Calculation")
    print("=" * 70)

    test_cases = [
        (QubitOperator('X0', 1.0), QubitOperator('Y0', 2.0)),
        (QubitOperator('X0', 1.0), QubitOperator('X0', 1.0)),
        (QubitOperator('X0 Y1', 0.5), QubitOperator('Y0 X1', 1.5)),
        (QubitOperator('Z0 Z1', 1.0 + 1.0j), QubitOperator('X0 X1', 2.0)),
    ]

    all_passed = True
    for op1, op2 in test_cases:
        # Original implementation
        _, norm_orig = compute_commutator(op1, op2)

        # Fast implementation
        key1 = list(op1.terms.keys())[0]
        key2 = list(op2.terms.keys())[0]
        c1 = list(op1.terms.values())[0]
        c2 = list(op2.terms.values())[0]

        x1, z1 = encode_pauli_string(key1, 10)
        x2, z2 = encode_pauli_string(key2, 10)

        norm_fast = compute_commutator_norm_fast(x1, z1, c1, x2, z2, c2)

        # Check if they match (within numerical precision)
        diff = abs(norm_orig - norm_fast)
        passed = diff < 1e-10

        status = "✓" if passed else "❌"
        print(f"\n  {status} [{op1}, {op2}]")
        print(f"     Original norm: {norm_orig:.6f}")
        print(f"     Fast norm:     {norm_fast:.6f}")
        print(f"     Difference:    {diff:.2e}")

        if not passed:
            all_passed = False

    if all_passed:
        print("\n  ✅ All commutator norm tests PASSED")
    else:
        print("\n  ❌ Some commutator norm tests FAILED")

    return all_passed


def test_statistical_agreement():
    """Test that Monte Carlo estimates agree statistically."""
    print("\n" + "=" * 70)
    print("TEST 4: Statistical Agreement (Monte Carlo)")
    print("=" * 70)

    from trotter_coefficients import trotter_error_estimator
    from trotter_coefficients_fast import trotter_error_estimator_fast

    # Generate small test case
    print("\n  Generating random Pauli terms...")
    np.random.seed(42)  # For reproducibility
    pauli_ops = ['X', 'Y', 'Z']
    n_terms = 50
    n_qubits = 10

    terms = []
    for _ in range(n_terms):
        weight = np.random.randint(1, 4)
        qubits = np.random.choice(n_qubits, size=weight, replace=False)
        ops = np.random.choice(pauli_ops, size=weight)
        pauli_string = ' '.join([f'{op}{q}' for q, op in zip(qubits, ops)])
        coeff = np.random.randn() + 1j * np.random.randn()
        terms.append(QubitOperator(pauli_string, coeff))

    print(f"  Generated {len(terms)} terms")

    # Run both implementations with same time limit
    time_limit = 3.0

    print(f"\n  Running original implementation (time_limit={time_limit}s)...")
    C1_orig, C2_orig = trotter_error_estimator(terms, time_limit)

    print(f"  Running fast implementation (time_limit={time_limit}s)...")
    C1_fast, C2_fast = trotter_error_estimator_fast(terms, time_limit, batch_size=1000)

    # Check if results are statistically close (within ~10% is reasonable for MC)
    C1_diff_pct = abs(C1_fast - C1_orig) / abs(C1_orig) * 100 if C1_orig != 0 else 0
    C2_diff_pct = abs(C2_fast - C2_orig) / abs(C2_orig) * 100 if C2_orig != 0 else 0

    print(f"\n  Results:")
    print(f"    C1 original: {C1_orig:.6f}")
    print(f"    C1 fast:     {C1_fast:.6f}")
    print(f"    Difference:  {C1_diff_pct:.2f}%")
    print()
    print(f"    C2 original: {C2_orig:.6f}")
    print(f"    C2 fast:     {C2_fast:.6f}")
    print(f"    Difference:  {C2_diff_pct:.2f}%")

    # Monte Carlo estimates can vary by ~5-10%, so we allow 15% tolerance
    passed = C1_diff_pct < 15.0 and C2_diff_pct < 15.0

    if passed:
        print("\n  ✅ Statistical agreement test PASSED")
        print("     (Both implementations produce similar Monte Carlo estimates)")
    else:
        print("\n  ⚠️  Statistical agreement test: larger than expected difference")
        print("     (This can happen with Monte Carlo - try increasing time_limit)")

    return passed


def run_all_tests():
    """Run all correctness tests."""
    print("\n" + "╔" + "═" * 68 + "╗")
    print("║" + " " * 16 + "CORRECTNESS TEST SUITE" + " " * 30 + "║")
    print("╚" + "═" * 68 + "╝")

    results = []

    # Ground truth tests (verify both implementations against known results)
    print("\n" + "┌" + "─" * 68 + "┐")
    print("│" + " " * 15 + "PART 1: GROUND TRUTH TESTS" + " " * 27 + "│")
    print("└" + "─" * 68 + "┘")
    results.append(("Ground Truth: Anticommutation", test_ground_truth_anticommutation()))
    results.append(("Ground Truth: Commutator Norms", test_ground_truth_commutator_norms()))
    results.append(("Ground Truth: Pauli Products", test_ground_truth_pauli_products()))

    # Consistency tests (compare fast vs original implementation)
    print("\n" + "┌" + "─" * 68 + "┐")
    print("│" + " " * 13 + "PART 2: IMPLEMENTATION CONSISTENCY" + " " * 21 + "│")
    print("└" + "─" * 68 + "┘")
    results.append(("Pauli String Encoding", test_encoding()))
    results.append(("Anticommutation Check", test_anticommutation()))
    results.append(("Commutator Norm", test_commutator_norm()))
    results.append(("Statistical Agreement", test_statistical_agreement()))

    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)

    all_passed = True
    for name, passed in results:
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"  {status}: {name}")
        if not passed:
            all_passed = False

    print()
    if all_passed:
        print("  🎉 All tests PASSED!")
        print("     - Both implementations satisfy ground truth properties")
        print("     - Fast and original implementations are consistent")
    else:
        print("  ⚠️  Some tests FAILED. Please review the implementations.")

    print("=" * 70)


if __name__ == "__main__":
    run_all_tests()
