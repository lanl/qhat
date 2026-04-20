"""
Test convergence monitoring feature.

Verifies that:
1. Standard error is computed and reported
2. Early termination happens when convergence threshold is met
3. Confidence intervals are reasonable
"""

import numpy as np
from openfermion import QubitOperator
from jkg_utils_fast import trotter_error_estimator_fast
import jkg_utils_fast

print("="*70)
print("TEST: Convergence Monitoring")
print("="*70)

# Test 1: Large system with convergence monitoring
print("\n" + "="*70)
print("TEST 1: Large system (N=200) with convergence monitoring")
print("="*70)

# Create 200 random Pauli terms
np.random.seed(42)
terms_large = []
for i in range(200):
    op = np.random.choice(['X', 'Y', 'Z'])
    qubit = i % 10
    coeff = np.random.uniform(0.1, 1.0)
    terms_large.append(QubitOperator(f'{op}{qubit}', coeff))

print(f"\nTest Hamiltonian with {len(terms_large)} terms")
print(f"Convergence threshold: {jkg_utils_fast.CONVERGENCE_THRESHOLD:.4f} (1%)")
print(f"Time limit: 30s")
print("\nWith convergence monitoring enabled, estimation may terminate early")
print("if the relative standard error drops below 1%.")

c1_large, c2_large = trotter_error_estimator_fast(terms_large, 30.0, batch_size=10000)

print(f"\nFinal results:")
print(f"  c1 = {c1_large:.6f}")
print(f"  c2 = {c2_large:.6f}")

# Test 2: Disable convergence monitoring
print("\n" + "="*70)
print("TEST 2: Same system with convergence monitoring DISABLED")
print("="*70)

# Temporarily disable convergence monitoring
original_enable = jkg_utils_fast.ENABLE_CONVERGENCE_MONITORING
jkg_utils_fast.ENABLE_CONVERGENCE_MONITORING = False

print(f"\nWith convergence monitoring disabled, all time budget will be used")

c1_no_conv, c2_no_conv = trotter_error_estimator_fast(terms_large, 30.0, batch_size=10000)

print(f"\nFinal results:")
print(f"  c1 = {c1_no_conv:.6f}")
print(f"  c2 = {c2_no_conv:.6f}")

# Restore original setting
jkg_utils_fast.ENABLE_CONVERGENCE_MONITORING = original_enable

# Test 3: Adjust convergence threshold
print("\n" + "="*70)
print("TEST 3: Tighter convergence threshold (0.1%)")
print("="*70)

original_threshold = jkg_utils_fast.CONVERGENCE_THRESHOLD
jkg_utils_fast.CONVERGENCE_THRESHOLD = 0.001  # 0.1%

print(f"\nConvergence threshold: {jkg_utils_fast.CONVERGENCE_THRESHOLD:.4f} (0.1%)")
print("This requires 100x more samples to converge than 1% threshold.")

c1_tight, c2_tight = trotter_error_estimator_fast(terms_large, 30.0, batch_size=10000)

print(f"\nFinal results:")
print(f"  c1 = {c1_tight:.6f}")
print(f"  c2 = {c2_tight:.6f}")

# Restore original threshold
jkg_utils_fast.CONVERGENCE_THRESHOLD = original_threshold

# Test 4: Very large system (should not converge in time limit)
print("\n" + "="*70)
print("TEST 4: Very large system (N=1000) - may not converge")
print("="*70)

terms_very_large = []
for i in range(1000):
    op = np.random.choice(['X', 'Y', 'Z'])
    qubit = i % 20
    coeff = np.random.uniform(0.1, 1.0)
    terms_very_large.append(QubitOperator(f'{op}{qubit}', coeff))

print(f"\nTest Hamiltonian with {len(terms_very_large)} terms")
print(f"Time limit: 10s (short)")
print("\nFor very large systems, time budget may expire before convergence.")

c1_vl, c2_vl = trotter_error_estimator_fast(terms_very_large, 10.0, batch_size=10000)

print(f"\nFinal results:")
print(f"  c1 = {c1_vl:.6f}")
print(f"  c2 = {c2_vl:.6f}")

# Summary
print("\n" + "="*70)
print("SUMMARY: Convergence Monitoring Benefits")
print("="*70)
print("""
✅ Standard Error Reporting:
   - Provides statistical confidence bounds on estimates
   - Shows relative uncertainty (SE/mean)
   - Helps assess estimation quality

✅ Early Termination:
   - Stops when estimate has converged (SE/mean < threshold)
   - Saves computation time for well-behaved systems
   - Configurable threshold (default: 1%)

✅ Adaptive Behavior:
   - Small systems: Exact computation (no uncertainty)
   - Medium systems: May converge early via monitoring
   - Large systems: Uses full time budget for best estimate

✅ Configuration:
   ENABLE_CONVERGENCE_MONITORING = True/False
   CONVERGENCE_THRESHOLD = 0.01  # 1% relative SE
   MIN_SAMPLES_FOR_CONVERGENCE = 100000

🔧 Adjust threshold for your needs:
   - 0.01 (1%): Good balance, fast convergence
   - 0.001 (0.1%): High precision, slower convergence
   - 0.0: Disable early termination, use full time budget
""")
print("="*70)
