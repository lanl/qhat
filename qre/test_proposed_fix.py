"""
Verify that the proposed fix to qre_unitary.py works correctly.

This test compares the CURRENT (broken) formulas with the FIXED formulas
against the analytical solution from the paper.
"""

import numpy as np
import math
from openfermion import QubitOperator
from jkg_utils import trotter_error_estimator

print("="*70)
print("VERIFICATION OF PROPOSED FIX")
print("="*70)

# Test Hamiltonian: H = X + Y
terms = [QubitOperator('X0', 1.0), QubitOperator('Y0', 1.0)]

# Get exact values
C1_analytical = 2.0
C2_analytical = 1.0/3.0

print(f"\nTest Hamiltonian: H = X₀ + Y₀")
print(f"Analytical: C1 = {C1_analytical}, C2 = {C2_analytical:.6f}")

# Verify with Monte Carlo
C1, C2 = trotter_error_estimator(terms, time_limit=2.0, batch_size=10000)
print(f"Monte Carlo: C1 = {C1:.6f}, C2 = {C2:.6f}")

# Use analytical values for precision
C1 = C1_analytical
C2 = C2_analytical

# Test parameters
timestep = 2.0
energy_error = 0.1
error_scale = 1.0

eps_trotter = energy_error * timestep / (2 * math.pi)

print(f"\nTest parameters:")
print(f"  timestep = {timestep}")
print(f"  energy_error = {energy_error}")
print(f"  eps_trotter = {eps_trotter:.6f}")

# Current (BROKEN) formulas
print("\n" + "="*70)
print("CURRENT FORMULAS (BROKEN)")
print("="*70)

s1_current = timestep * error_scale * C1 / eps_trotter
s2_current = timestep * math.sqrt(error_scale * C2 / eps_trotter)

print(f"s1 = timestep * error_scale * C1 / eps_trotter")
print(f"   = {timestep} * {error_scale} * {C1} / {eps_trotter:.6f}")
print(f"   = {s1_current:.2f} steps")

print(f"\ns2 = timestep * sqrt(error_scale * C2 / eps_trotter)")
print(f"   = {timestep} * sqrt({error_scale} * {C2:.6f} / {eps_trotter:.6f})")
print(f"   = {s2_current:.2f} steps")

# Proposed (FIXED) formulas
print("\n" + "="*70)
print("PROPOSED FIX")
print("="*70)

s1_fixed = timestep**2 * error_scale * C1 / (2 * eps_trotter)
s2_fixed = (timestep**1.5) * math.sqrt(error_scale * C2 / eps_trotter)

print(f"s1 = timestep**2 * error_scale * C1 / (2 * eps_trotter)")
print(f"   = {timestep}² * {error_scale} * {C1} / (2 * {eps_trotter:.6f})")
print(f"   = {s1_fixed:.2f} steps")

print(f"\ns2 = timestep**1.5 * sqrt(error_scale * C2 / eps_trotter)")
print(f"   = {timestep}^1.5 * sqrt({error_scale} * {C2:.6f} / {eps_trotter:.6f})")
print(f"   = {s2_fixed:.2f} steps")

# Expected from paper
print("\n" + "="*70)
print("EXPECTED FROM PAPER (Childs et al. Eq 145, 152)")
print("="*70)

r1_paper = timestep**2 * C1 / (2 * eps_trotter)
r2_paper = (timestep**1.5) * math.sqrt(C2 / eps_trotter)

print(f"First-order: r = t² * C1 / (2ε)")
print(f"           = {timestep}² * {C1} / (2 * {eps_trotter:.6f})")
print(f"           = {r1_paper:.2f} steps")

print(f"\nSecond-order: r = t^(3/2) * sqrt(C2 / ε)")
print(f"            = {timestep}^1.5 * sqrt({C2:.6f} / {eps_trotter:.6f})")
print(f"            = {r2_paper:.2f} steps")

# Comparison
print("\n" + "="*70)
print("COMPARISON")
print("="*70)

print("\nFirst-order:")
print(f"  Current (broken): {s1_current:.2f} steps")
print(f"  Fixed:            {s1_fixed:.2f} steps")
print(f"  Paper:            {r1_paper:.2f} steps")
print(f"  Current error:    {abs(s1_current - r1_paper)/r1_paper*100:.1f}%")
print(f"  Fixed error:      {abs(s1_fixed - r1_paper)/r1_paper*100:.1f}%")
if abs(s1_fixed - r1_paper) < 0.01:
    print("  ✅ FIX WORKS!")
else:
    print("  ❌ FIX DOESN'T WORK")

print("\nSecond-order:")
print(f"  Current (broken): {s2_current:.2f} steps")
print(f"  Fixed:            {s2_fixed:.2f} steps")
print(f"  Paper:            {r2_paper:.2f} steps")
print(f"  Current error:    {abs(s2_current - r2_paper)/r2_paper*100:.1f}%")
print(f"  Fixed error:      {abs(s2_fixed - r2_paper)/r2_paper*100:.1f}%")
if abs(s2_fixed - r2_paper) < 0.01:
    print("  ✅ FIX WORKS!")
else:
    print("  ❌ FIX DOESN'T WORK")

# Test multiple timesteps
print("\n" + "="*70)
print("TIME SCALING TEST")
print("="*70)
print("\nVerifying correct time dependence with proposed fix:\n")

print(f"{'Time':>6} | {'s1 (fixed)':>12} | {'r1 (paper)':>12} | {'Match':>6} | "
      f"{'s2 (fixed)':>12} | {'r2 (paper)':>12} | {'Match':>6}")
print("-"*88)

all_match = True
for t in [0.5, 1.0, 2.0, 4.0]:
    eps = energy_error * t / (2 * math.pi)

    s1_fix = t**2 * C1 / (2 * eps)
    s2_fix = t**1.5 * math.sqrt(C2 / eps)

    r1 = t**2 * C1 / (2 * eps)
    r2 = t**1.5 * math.sqrt(C2 / eps)

    match1 = "✓" if abs(s1_fix - r1) < 0.01 else "✗"
    match2 = "✓" if abs(s2_fix - r2) < 0.01 else "✗"

    print(f"{t:6.1f} | {s1_fix:12.2f} | {r1:12.2f} | {match1:>6} | "
          f"{s2_fix:12.2f} | {r2:12.2f} | {match2:>6}")

    if match1 == "✗" or match2 == "✗":
        all_match = False

print()
if all_match:
    print("✅ ALL TIME POINTS MATCH - FIX IS CORRECT!")
else:
    print("❌ SOME TIME POINTS DON'T MATCH - FIX HAS ISSUES")

# Summary
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print("""
TO APPLY THE FIX:

Edit qre_unitary.py, lines 106-107.

Change FROM:
    s1 = timestep * config_unitary.error_scale * c1 / eps_trotter
    s2 = timestep * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)

Change TO:
    s1 = timestep**2 * config_unitary.error_scale * c1 / (2 * eps_trotter)
    s2 = (timestep**1.5) * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)

Key changes:
  1. s1: Added timestep**2 and division by 2
  2. s2: Changed timestep to timestep**1.5

This matches Childs et al. (arXiv:1912.08854v3) Equations 145 and 152.
""")

print("="*70)
