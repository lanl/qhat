"""
Test Phase 2 & 3: Test mode selection and backward compatibility.
"""

import numpy as np
from openfermion import QubitOperator
from trotter_coefficients_fast import trotter_error_estimator_fast

print("="*70)
print("PHASE 2 & 3 TESTING: Mode Selection and Backward Compatibility")
print("="*70)

# Test 1: Backward compatibility - default behavior unchanged
print("\nTest 1: Backward compatibility (default parameters)")
H1 = QubitOperator('X0', 1.0)
H2 = QubitOperator('Y0', 1.0)
terms = [H1, H2]

# Default should be monte_carlo with auto_exact=False
c1, c2 = trotter_error_estimator_fast(terms, 10)
print(f"  Result: C1={c1:.6f}, C2={c2:.6f}")
# Expected: C1 = 2/2 = 1.0, C2 = 0/12 + 4/24 = 1/6 ≈ 0.166667
print(f"  Expected: C1=1.0, C2=0.166667")
assert abs(c1 - 1.0) < 1e-6, f"C1: expected 1.0, got {c1}"
assert abs(c2 - 1.0/6.0) < 1e-6, f"C2: expected {1.0/6.0}, got {c2}"
print("  ✅ Test 1 passed!")

# Test 2: Explicit monte_carlo mode
print("\nTest 2: Explicit mode='monte_carlo'")
c1, c2 = trotter_error_estimator_fast(terms, 10, mode='monte_carlo')
print(f"  Result: C1={c1:.6f}, C2={c2:.6f}")
assert abs(c1 - 1.0) < 1e-6
assert abs(c2 - 1.0/6.0) < 1e-6
print("  ✅ Test 2 passed!")

# Test 3: Exact mode on small system
print("\nTest 3: mode='exact' on small system")
c1, c2 = trotter_error_estimator_fast(terms, 10, mode='exact')
print(f"  Result: C1={c1:.6f}, C2={c2:.6f}")
assert abs(c1 - 1.0) < 1e-10, f"C1: expected 1.0, got {c1}"
assert abs(c2 - 1.0/6.0) < 1e-10, f"C2: expected {1.0/6.0}, got {c2}"
print("  ✅ Test 3 passed!")

# Test 4: Auto-exact mode on small system
print("\nTest 4: mode='monte_carlo' with auto_exact=True")
c1, c2 = trotter_error_estimator_fast(terms, 60, mode='monte_carlo', auto_exact=True)
print(f"  Result: C1={c1:.6f}, C2={c2:.6f}")
assert abs(c1 - 1.0) < 1e-10, f"C1: expected 1.0, got {c1}"
assert abs(c2 - 1.0/6.0) < 1e-10, f"C2: expected {1.0/6.0}, got {c2}"
print("  ✅ Test 4 passed!")

# Test 5: Invalid mode raises error
print("\nTest 5: Invalid mode raises ValueError")
try:
    c1, c2 = trotter_error_estimator_fast(terms, 10, mode='invalid')
    print("  ❌ Should have raised ValueError!")
    assert False
except ValueError as e:
    print(f"  Caught expected error: {e}")
    print("  ✅ Test 5 passed!")

# Test 6: Exact mode on large system raises error
print("\nTest 6: mode='exact' on large system raises ValueError")
# Create a larger system (500 terms - infeasible due to time/memory)
large_terms = []
for i in range(500):
    qubit = i % 50
    pauli = ['X', 'Y', 'Z'][i % 3]
    large_terms.append(QubitOperator(f'{pauli}{qubit}', 1.0))

try:
    c1, c2 = trotter_error_estimator_fast(large_terms, 10, mode='exact')
    print("  ❌ Should have raised ValueError!")
    assert False
except ValueError as e:
    print(f"  Caught expected error: {str(e)[:80]}...")
    print("  ✅ Test 6 passed!")

# Test 7: Monte Carlo mode works on large system
print("\nTest 7: mode='monte_carlo' on large system")
# Use 200 terms for Monte Carlo test
medium_terms = []
for i in range(200):
    qubit = i % 50
    pauli = ['X', 'Y', 'Z'][i % 3]
    medium_terms.append(QubitOperator(f'{pauli}{qubit}', 1.0))
c1, c2 = trotter_error_estimator_fast(medium_terms, 10, mode='monte_carlo')
print(f"  Result: C1={c1:.6f}, C2={c2:.6f}")
assert c1 > 0 and c2 >= 0
print("  ✅ Test 7 passed!")

# Test 8: All modes agree on small system
print("\nTest 8: All modes agree on small analytical system")
H1 = QubitOperator('X0', 1.0)
H2 = QubitOperator('Y0', 1.0)
H3 = QubitOperator('Z0', 1.0)
terms = [H1, H2, H3]

c1_mc, c2_mc = trotter_error_estimator_fast(terms, 20, mode='monte_carlo', auto_exact=False)
c1_ex, c2_ex = trotter_error_estimator_fast(terms, 20, mode='exact')
c1_auto, c2_auto = trotter_error_estimator_fast(terms, 60, mode='monte_carlo', auto_exact=True)

print(f"  Monte Carlo: C1={c1_mc:.6f}, C2={c2_mc:.6f}")
print(f"  Exact:       C1={c1_ex:.6f}, C2={c2_ex:.6f}")
print(f"  Auto-exact:  C1={c1_auto:.6f}, C2={c2_auto:.6f}")

# Expected: C1 = 3.0, C2 = 0.5
assert abs(c1_ex - 3.0) < 1e-10, f"Exact C1: expected 3.0, got {c1_ex}"
assert abs(c2_ex - 0.5) < 1e-10, f"Exact C2: expected 0.5, got {c2_ex}"
assert abs(c1_auto - 3.0) < 1e-10, f"Auto C1: expected 3.0, got {c1_auto}"
assert abs(c2_auto - 0.5) < 1e-10, f"Auto C2: expected 0.5, got {c2_auto}"
# Monte Carlo should be close
assert abs(c1_mc - 3.0) / 3.0 < 0.1, f"MC C1 differs >10% from exact"
assert abs(c2_mc - 0.5) / 0.5 < 0.1, f"MC C2 differs >10% from exact"

print("  ✅ Test 8 passed!")

print("\n" + "="*70)
print("✅ ALL PHASE 2 & 3 TESTS PASSED!")
print("="*70)
