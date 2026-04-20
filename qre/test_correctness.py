"""
Test script to verify correctness of the optimized implementation.

This script tests the low-level functions to ensure they produce
correct results compared to the original implementation.
"""

import numpy as np
from openfermion import QubitOperator

# Import original functions
from jkg_utils import anticommute, pauli_product_key, compute_commutator

# Import optimized functions
from jkg_utils_fast import (
    encode_pauli_string,
    pauli_anticommute,
    pauli_product_bits,
    compute_commutator_norm_fast
)


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

    from jkg_utils import trotter_error_estimator
    from jkg_utils_fast import trotter_error_estimator_fast

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
        print("  🎉 All tests PASSED! The optimized implementation is correct.")
    else:
        print("  ⚠️  Some tests FAILED. Please review the implementation.")

    print("=" * 70)


if __name__ == "__main__":
    run_all_tests()
