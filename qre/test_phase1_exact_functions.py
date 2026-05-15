"""
Test Phase 1: Test the exact computation helper functions standalone.
"""

import numpy as np
from openfermion import QubitOperator
from trotter_coefficients_fast import (
    preprocess_pauli_terms,
    generate_all_pairs,
    generate_all_triples,
    compute_C1_exact,
    compute_C21_exact,
    compute_C22_exact
)

print("="*70)
print("PHASE 1 TESTING: Exact Computation Functions")
print("="*70)

# Test 1: Generate pairs and triples for small N
print("\nTest 1: Index generation for N=5")
pairs = generate_all_pairs(5)
print(f"  Generated {len(pairs)} pairs (expected {5*4//2} = 10)")
print(f"  Pairs: {pairs.tolist()}")
assert len(pairs) == 10

triples = generate_all_triples(5)
print(f"  Generated {len(triples)} triples")
print(f"  First 5 triples: {triples[:5].tolist()}")

# Test 2: Simple analytical case - X + Y on same qubit
print("\nTest 2: Analytical test with H = X₀ + Y₀")
H1 = QubitOperator('X0', 1.0)
H2 = QubitOperator('Y0', 1.0)
terms = [H1, H2]

x_bits, z_bits, coeffs, n_qubits = preprocess_pauli_terms(terms)
print(f"  Preprocessed: N=2, n_qubits={n_qubits}")

# Expected: C1 = ||[X, Y]|| = 2, so C1_exact should return 2.0
C1_exact = compute_C1_exact(x_bits, z_bits, coeffs, 2)
print(f"  C1_exact = {C1_exact:.6f} (expected 2.0)")
assert abs(C1_exact - 2.0) < 1e-10, f"Expected 2.0, got {C1_exact}"

# For N=2, no triples, so C21 = 0
C21_exact = compute_C21_exact(x_bits, z_bits, coeffs, 2)
print(f"  C21_exact = {C21_exact:.6f} (expected 0.0)")
assert abs(C21_exact - 0.0) < 1e-10, f"Expected 0.0, got {C21_exact}"

# For N=2, C22 = ||[H_0, [H_0, H_1]]||
# [X₀, Y₀] = 2iZ₀
# [X₀, 2iZ₀] = 2i[X₀, Z₀] = 2i(-2iY₀) = 4Y₀
# ||4Y₀|| = 4
C22_exact = compute_C22_exact(x_bits, z_bits, coeffs, 2)
print(f"  C22_exact = {C22_exact:.6f} (expected 4.0)")
assert abs(C22_exact - 4.0) < 1e-10, f"Expected 4.0, got {C22_exact}"

print("  ✅ Test 2 passed!")

# Test 3: Three Paulis on same qubit
print("\nTest 3: Analytical test with H = X₀ + Y₀ + Z₀")
H1 = QubitOperator('X0', 1.0)
H2 = QubitOperator('Y0', 1.0)
H3 = QubitOperator('Z0', 1.0)
terms = [H1, H2, H3]

x_bits, z_bits, coeffs, n_qubits = preprocess_pauli_terms(terms)

# Expected: C1 = 2 + 2 + 2 = 6
C1_exact = compute_C1_exact(x_bits, z_bits, coeffs, 3)
print(f"  C1_exact = {C1_exact:.6f} (expected 6.0)")
assert abs(C1_exact - 6.0) < 1e-10, f"Expected 6.0, got {C1_exact}"

C21_exact = compute_C21_exact(x_bits, z_bits, coeffs, 3)
print(f"  C21_exact = {C21_exact:.6f}")

C22_exact = compute_C22_exact(x_bits, z_bits, coeffs, 3)
print(f"  C22_exact = {C22_exact:.6f} (expected 12.0)")
assert abs(C22_exact - 12.0) < 1e-10, f"Expected 12.0, got {C22_exact}"

print("  ✅ Test 3 passed!")

# Test 4: Commuting operators (different qubits)
print("\nTest 4: Commuting operators H = X₀ + X₁")
H1 = QubitOperator('X0', 1.0)
H2 = QubitOperator('X1', 1.0)
terms = [H1, H2]

x_bits, z_bits, coeffs, n_qubits = preprocess_pauli_terms(terms)

# Expected: all commute, so C1 = 0, C21 = 0, C22 = 0
C1_exact = compute_C1_exact(x_bits, z_bits, coeffs, 2)
print(f"  C1_exact = {C1_exact:.6f} (expected 0.0)")
assert abs(C1_exact - 0.0) < 1e-10, f"Expected 0.0, got {C1_exact}"

C21_exact = compute_C21_exact(x_bits, z_bits, coeffs, 2)
print(f"  C21_exact = {C21_exact:.6f} (expected 0.0)")
assert abs(C21_exact - 0.0) < 1e-10, f"Expected 0.0, got {C21_exact}"

C22_exact = compute_C22_exact(x_bits, z_bits, coeffs, 2)
print(f"  C22_exact = {C22_exact:.6f} (expected 0.0)")
assert abs(C22_exact - 0.0) < 1e-10, f"Expected 0.0, got {C22_exact}"

print("  ✅ Test 4 passed!")

print("\n" + "="*70)
print("✅ ALL PHASE 1 TESTS PASSED!")
print("="*70)
