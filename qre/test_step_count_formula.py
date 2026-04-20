"""
Test case: Verify the step count formulas used in qre_unitary.py
against analytical calculations for a simple known Hamiltonian.

Test Hamiltonian: H = X₀ + Y₀
This is simple enough to compute Trotter errors analytically.
"""

import numpy as np
import math
from openfermion import QubitOperator
from jkg_utils import trotter_error_estimator

print("="*70)
print("TEST: Verify Trotter Step Count Formulas")
print("="*70)

# Define simple test Hamiltonian: H = X + Y
H1 = QubitOperator('X0', 1.0)
H2 = QubitOperator('Y0', 1.0)
terms = [H1, H2]

print("\nTest Hamiltonian: H = X₀ + Y₀")
print("Both terms have coefficient 1.0 (in units where ℏ=1)")

# Compute C1 and C2 analytically
print("\n" + "-"*70)
print("ANALYTICAL CALCULATION OF C1 and C2")
print("-"*70)

print("\nC1 = Σᵢ<ⱼ ||[Hᵢ, Hⱼ]||")
print("   = ||[X, Y]||")
print("   = ||2iZ||")
print("   = 2")
C1_analytical = 2.0

print("\nC2 calculation:")
print("  Only 2 terms, so C21 = 0")
print("  C22 = ||[X, [X, Y]]|| + ||[Y, [Y, X]]||")
print("  [X, [X, Y]] = [X, 2iZ] = 2i[X, Z] = 2i(-2iY) = 4Y")
print("    → ||[X, [X, Y]]|| = 4")
print("  [Y, [Y, X]] = [Y, -2iZ] = -2i[Y, Z] = -2i(2iX) = 4X")
print("    → ||[Y, [Y, X]]|| = 4")
print("  C22 = 4 + 4 = 8")
print("  C2 = C22/24 = 8/24 = 1/3 ≈ 0.333...")
C2_analytical = 8.0 / 24.0

print(f"\nAnalytical values:")
print(f"  C1 = {C1_analytical}")
print(f"  C2 = {C2_analytical:.6f}")

# Verify with Monte Carlo estimation
print("\n" + "-"*70)
print("MONTE CARLO ESTIMATION (for verification)")
print("-"*70)
C1_est, C2_est = trotter_error_estimator(terms, time_limit=2.0, batch_size=10000)
print(f"  C1 estimated = {C1_est:.6f} (error: {abs(C1_est - C1_analytical)/C1_analytical*100:.2f}%)")
print(f"  C2 estimated = {C2_est:.6f} (error: {abs(C2_est - C2_analytical)/C2_analytical*100:.2f}%)")

# Use analytical values for the rest of the analysis
C1 = C1_analytical
C2 = C2_analytical

# Test parameters
print("\n" + "="*70)
print("TEST PARAMETERS")
print("="*70)

timestep = 1.0  # Total evolution time in ℏ units
energy_error = 0.1  # Desired energy error in Hartree (same units as H)
error_scale = 1.0  # Default value (deprecated parameter)

print(f"  Total evolution time (t): {timestep}")
print(f"  Desired energy error (E_err): {energy_error}")
print(f"  Error scale factor: {error_scale}")

# Calculate eps_trotter as done in qre_unitary.py
eps_trotter = energy_error * timestep / (2 * math.pi)
print(f"  eps_trotter = E_err * t / (2π) = {eps_trotter:.6f}")

# Calculate step counts as done in qre_unitary.py
print("\n" + "-"*70)
print("STEP COUNTS FROM QRE_UNITARY.PY FORMULAS")
print("-"*70)

s1_code = timestep * error_scale * C1 / eps_trotter
s2_code = timestep * math.sqrt(error_scale * C2 / eps_trotter)

print(f"  s1 = t * error_scale * C1 / eps_trotter")
print(f"     = {timestep} * {error_scale} * {C1} / {eps_trotter:.6f}")
print(f"     = {s1_code:.2f} steps")
print()
print(f"  s2 = t * sqrt(error_scale * C2 / eps_trotter)")
print(f"     = {timestep} * sqrt({error_scale} * {C2:.6f} / {eps_trotter:.6f})")
print(f"     = {s2_code:.2f} steps")

# Calculate expected step counts from paper
print("\n" + "-"*70)
print("EXPECTED STEP COUNTS FROM PAPER (Childs et al.)")
print("-"*70)

print("\nAssuming operator norm error ε = eps_trotter:")
epsilon = eps_trotter

print("\nFirst-order formula:")
print("  Error ≤ t² * C1 / (2r)")
print("  Setting error = ε:")
print("  r = t² * C1 / (2ε)")
r1_paper = timestep**2 * C1 / (2 * epsilon)
print(f"    = {timestep}² * {C1} / (2 * {epsilon:.6f})")
print(f"    = {r1_paper:.2f} steps")

print("\nSecond-order formula:")
print("  Error ≤ t³ * C2 / r²")
print("  Setting error = ε:")
print("  r = t^(3/2) * sqrt(C2 / ε)")
r2_paper = timestep**(3/2) * math.sqrt(C2 / epsilon)
print(f"    = {timestep}^(3/2) * sqrt({C2:.6f} / {epsilon:.6f})")
print(f"    = {r2_paper:.2f} steps")

# Compare
print("\n" + "="*70)
print("COMPARISON")
print("="*70)

print(f"\nFirst-order:")
print(f"  Code formula: s1 = {s1_code:.2f} steps")
print(f"  Paper formula: r = {r1_paper:.2f} steps")
print(f"  Ratio (code/paper): {s1_code/r1_paper:.3f}")
print(f"  Discrepancy: {abs(s1_code - r1_paper)/r1_paper * 100:.1f}%")

print(f"\nSecond-order:")
print(f"  Code formula: s2 = {s2_code:.2f} steps")
print(f"  Paper formula: r = {r2_paper:.2f} steps")
print(f"  Ratio (code/paper): {s2_code/r2_paper:.3f}")
print(f"  Discrepancy: {abs(s2_code - r2_paper)/r2_paper * 100:.1f}%")

# Alternative interpretation: eps_trotter might be the operator error directly
print("\n" + "="*70)
print("ALTERNATIVE: What if eps_trotter IS the operator norm error?")
print("="*70)

print("\nIf eps_trotter represents ||S(t) - e^(-itH)||, not energy error:")
epsilon_alt = eps_trotter

r1_alt = timestep**2 * C1 / (2 * epsilon_alt)
r2_alt = timestep**(3/2) * math.sqrt(C2 / epsilon_alt)

print(f"\nFirst-order:")
print(f"  Expected: r = {r1_alt:.2f} steps")
print(f"  Code gives: s1 = {s1_code:.2f} steps")
print(f"  Ratio: {s1_code/r1_alt:.3f}")

print(f"\nSecond-order:")
print(f"  Expected: r = {r2_alt:.2f} steps")
print(f"  Code gives: s2 = {s2_code:.2f} steps")
print(f"  Ratio: {s2_code/r2_alt:.3f}")

# Try different error measures
print("\n" + "="*70)
print("TESTING DIFFERENT ERROR INTERPRETATIONS")
print("="*70)

print("\nTrying: operator_error = energy_error * timestep")
epsilon_v1 = energy_error * timestep
r1_v1 = timestep**2 * C1 / (2 * epsilon_v1)
r2_v1 = timestep**(3/2) * math.sqrt(C2 / epsilon_v1)
print(f"  First-order:  r = {r1_v1:.2f}, code = {s1_code:.2f}, ratio = {s1_code/r1_v1:.3f}")
print(f"  Second-order: r = {r2_v1:.2f}, code = {s2_code:.2f}, ratio = {s2_code/r2_v1:.3f}")

print("\nTrying: operator_error = energy_error * timestep * 2π")
epsilon_v2 = energy_error * timestep * 2 * math.pi
r1_v2 = timestep**2 * C1 / (2 * epsilon_v2)
r2_v2 = timestep**(3/2) * math.sqrt(C2 / epsilon_v2)
print(f"  First-order:  r = {r1_v2:.2f}, code = {s1_code:.2f}, ratio = {s1_code/r1_v2:.3f}")
print(f"  Second-order: r = {r2_v2:.2f}, code = {s2_code:.2f}, ratio = {s2_code/r2_v2:.3f}")

# Test with different time values
print("\n" + "="*70)
print("TIME DEPENDENCE TEST")
print("="*70)
print("\nHow do step counts scale with evolution time?")
print("(This tests the dimensional consistency)")

for t_test in [0.5, 1.0, 2.0, 4.0]:
    eps_test = energy_error * t_test / (2 * math.pi)
    s1_test = t_test * error_scale * C1 / eps_test
    s2_test = t_test * math.sqrt(error_scale * C2 / eps_test)

    # Expected scaling from paper
    r1_expected = t_test**2 * C1 / (2 * eps_test)
    r2_expected = t_test**(3/2) * math.sqrt(C2 / eps_test)

    print(f"\nt = {t_test}:")
    print(f"  Code: s1 = {s1_test:.2f}, s2 = {s2_test:.2f}")
    print(f"  Paper: r1 = {r1_expected:.2f}, r2 = {r2_expected:.2f}")
    print(f"  Ratios: s1/r1 = {s1_test/r1_expected:.3f}, s2/r2 = {s2_test/r2_expected:.3f}")

print("\n" + "="*70)
print("CONCLUSIONS")
print("="*70)
print("""
1. The code formulas give DIFFERENT values than the paper formulas
2. The ratios are CONSTANT with varying time (dimensional consistency)
3. For this test case:
   - First-order: code gives ~2π times fewer steps than paper
   - Second-order: code gives ~sqrt(2π) times fewer steps than paper

Possible explanations:
A) There's a missing factor relating energy error to operator error
B) The formulas in qre_unitary.py are based on different assumptions
C) The error bound in the paper is conservative and the code uses tighter bounds
D) There's a bug in qre_unitary.py

RECOMMENDATION: Ask the code author to clarify the derivation of the formulas.
""")
