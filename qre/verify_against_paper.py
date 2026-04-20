"""
Verification script to help compare implementation against Childs et al paper.

This script computes Trotter error coefficients for a small, hand-calculable case
to help verify the formulas are correct.
"""

import numpy as np
from openfermion import QubitOperator
from trotter_coefficients import compute_commutator, compute_nested_commutator_norm

print("="*70)
print("VERIFICATION: Simple Test Case for Manual Calculation")
print("="*70)

# Simple test case: 3 Pauli terms
# H = X₀ + Y₀ + Z₀
H1 = QubitOperator('X0', 1.0)
H2 = QubitOperator('Y0', 1.0)
H3 = QubitOperator('Z0', 1.0)

terms = [H1, H2, H3]
N = len(terms)

print(f"\nHamiltonian: H = X₀ + Y₀ + Z₀")
print(f"Number of terms: N = {N}")

# ============================================================================
# Compute C1 = ∑ᵢ<ⱼ ||[Hᵢ, Hⱼ]||
# ============================================================================
print("\n" + "-"*70)
print("C1 = ∑ᵢ<ⱼ ||[Hᵢ, Hⱼ]||")
print("-"*70)

C1_exact = 0.0
for i in range(N):
    for j in range(i+1, N):
        _, norm = compute_commutator(terms[i], terms[j])
        print(f"  ||[H_{i}, H_{j}]|| = {norm:.4f}")
        C1_exact += norm

print(f"\nC1_exact = {C1_exact:.4f}")
print(f"\nNote: All pairs of X, Y, Z anticommute, so:")
print(f"  [X, Y] = 2iZ  →  ||[X, Y]|| = 2")
print(f"  [Y, Z] = 2iX  →  ||[Y, Z]|| = 2")
print(f"  [Z, X] = 2iY  →  ||[Z, X]|| = 2")
print(f"Expected: C1 = 3 × 2 = 6")

# ============================================================================
# Compute C21 = ∑ₖ<ⱼ,ₖ<ᵢ ||[Hᵢ, [Hⱼ, Hₖ]]||
# ============================================================================
print("\n" + "-"*70)
print("C21 = ∑ₖ<ⱼ,ₖ<ᵢ ||[Hᵢ, [Hⱼ, Hₖ]]||")
print("-"*70)

C21_exact = 0.0
count_c21 = 0
for k in range(N):
    for i in range(k+1, N):
        for j in range(k+1, N):
            if i == j:
                continue
            # Compute [Hⱼ, Hₖ]
            inner_comm, _ = compute_commutator(terms[j], terms[k])
            # Compute ||[Hᵢ, [Hⱼ, Hₖ]]||
            norm = compute_nested_commutator_norm(terms[i], inner_comm)
            if norm > 1e-10:
                print(f"  ||[H_{i}, [H_{j}, H_{k}]]|| = {norm:.4f}")
            C21_exact += norm
            count_c21 += 1

print(f"\nC21_exact = {C21_exact:.4f}")
print(f"Number of (i,j,k) triples with k<i,k<j: {count_c21}")

# ============================================================================
# Compute C22 = ∑ₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]||
# ============================================================================
print("\n" + "-"*70)
print("C22 = ∑ₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]||")
print("-"*70)

C22_exact = 0.0
for k in range(N):
    for j in range(k+1, N):
        # Compute [Hₖ, Hⱼ]
        inner_comm, _ = compute_commutator(terms[k], terms[j])
        # Compute ||[Hₖ, [Hₖ, Hⱼ]]||
        norm = compute_nested_commutator_norm(terms[k], inner_comm)
        if norm > 1e-10:
            print(f"  ||[H_{k}, [H_{k}, H_{j}]]|| = {norm:.4f}")
        C22_exact += norm

print(f"\nC22_exact = {C22_exact:.4f}")

# ============================================================================
# Compute what the code returns
# ============================================================================
print("\n" + "="*70)
print("WHAT THE CODE RETURNS")
print("="*70)

c1_code = C1_exact
c2_code = C21_exact / 12 + C22_exact / 24

print(f"c1 = {c1_code:.4f}")
print(f"c2 = C21/12 + C22/24 = {C21_exact:.4f}/12 + {C22_exact:.4f}/24 = {c2_code:.4f}")

# ============================================================================
# Analysis
# ============================================================================
print("\n" + "="*70)
print("ANALYSIS")
print("="*70)

print("\nFor the Hamiltonian H = X₀ + Y₀ + Z₀:")
print("  • All three Pauli operators are on the same qubit")
print("  • All pairs anticommute: [X,Y], [Y,Z], [Z,X]")
print("  • Products: XY=iZ, YZ=iX, ZX=iY")
print()
print("Expected commutators:")
print("  [X, Y] = 2iZ  →  ||[X, Y]|| = 2")
print("  [Y, Z] = 2iX  →  ||[Y, Z]|| = 2")
print("  [Z, X] = 2iY  →  ||[Z, X]|| = 2")
print()
print("Nested commutators:")
print("  [Z, [X, Y]] = [Z, 2iZ] = 0  (Z commutes with itself)")
print("  Similarly, most nested commutators are 0")
print()
print("TO VERIFY AGAINST PAPER:")
print("  1. Check if Equation 145 defines C₁ as the sum over i<j or i≠j")
print("  2. Verify the coefficients 1/12 and 1/24 for C21 and C22")
print("  3. Check dimensional consistency with error bounds")

# ============================================================================
# Run the actual code
# ============================================================================
print("\n" + "="*70)
print("RUNNING ACTUAL CODE")
print("="*70)

from trotter_coefficients import trotter_error_estimator
c1_est, c2_est = trotter_error_estimator(terms, time_limit=2.0, batch_size=10000)

print(f"\nMonte Carlo estimates:")
print(f"  c1_est = {c1_est:.4f}")
print(f"  c2_est = {c2_est:.4f}")
print(f"\nExact values:")
print(f"  c1_exact = {c1_code:.4f}")
print(f"  c2_exact = {c2_code:.4f}")
print(f"\nRelative errors:")
print(f"  c1: {abs(c1_est - c1_code) / c1_code * 100:.2f}%")
print(f"  c2: {abs(c2_est - c2_code) / c2_code * 100 if c2_code > 0 else 0:.2f}%")

print("\n" + "="*70)
print("INSTRUCTIONS FOR MANUAL VERIFICATION")
print("="*70)
print("""
To verify against the paper:

1. Open the PDF: /vast/home/bkkrueger/.claude/projects/.../webfetch-***.pdf
2. Look at Equation 145 (around page 30-35 typically)
3. Check:
   a) Is C₁ defined as ∑ᵢ<ⱼ or (1/2)∑ᵢ≠ⱼ ?
   b) What are the exact coefficients for the nested commutator terms?
   c) How does the error bound scale with these coefficients?

4. If C₁ in the paper uses (1/2)∑ᵢ≠ⱼ notation, then:
   - The code should return C1_est/2 to match the paper's definition
   - Update line 209 in trotter_coefficients.py to: return C1_est/2, C21_est/12 + C22_est/24

5. Verify the 1/12 and 1/24 coefficients match the paper
""")
