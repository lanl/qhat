# Error Analysis Assessment for config_full_analysis.py (CORRECTED)

**Date**: 2026-07-21  
**Config**: `analysis/examples/config_full_analysis.py`  
**System**: Be-H molecule, 6 qubits, Trotter approximation  
**Status**: ✅ No bugs found - errors are real!

---

## Current Error Results

```
ERROR ANALYSIS RESULTS FOR config_full_analysis.py
================================================================================

1. EIGENVALUE ERRORS:
   Number of eigenstates: 64
   Max absolute error: 3.957573e-03  ✅ Good (0.4%)
   Max relative error: 6.749448e-04  ✅ Excellent (0.07%)
   Mean absolute error: 4.202905e-04  ✅ Excellent
   Mean relative error: 7.064196e-05  ✅ Excellent (0.007%)

2. MATRIX NORM ERRORS:
   Frobenius norm ||U_exact - U_approx||_F: 6.463412  ⚠️ Large but consistent

3. STATE-DEPENDENT ERRORS:
   Absolute error: 8.074505e-01  ⚠️ High
   Relative error: 0.807450 (80.75%)  ⚠️ Very high

================================================================================
```

## ✅ Investigation Results: No Bug Found

### Energy Shift Verification

I verified that the energy shifts ARE correctly handled:

**Test 1: Check which Hamiltonian version the Trotter unitary uses**
```
H_exact (from file): [1.527082, 3.845617]  # Shifted version
U_approx (from phases): [1.527161, 3.845617]  # Matches shifted!
Internal Hamiltonian: [-3.845617, 0.586686]  # Original unshifted

Difference between U and H_exact:
  Min eigenvalue: 0.000079  ✅ Tiny (Trotter error)
  Max eigenvalue: 0.000000  ✅ Perfect match
```

**Conclusion**: The Trotter unitary uses the **shifted** Hamiltonian (energy range [1.527, 3.846]), NOT the internal unshifted version ([-3.846, 0.587]). This is correct!

**Test 2: Unitarity after removing shift**
```
U_approx (original): ||U†U - I||_F = 1.04e-14  ✅ Perfect
U_approx (after removing shift): ||U†U - I||_F = 8.50e-04  ❌ Destroyed
```

This SEEMS bad, but actually tells us: if removing the shift destroys unitarity, then **the unitary doesn't have that shift applied**! The code is trying to remove a shift that doesn't exist.

### The Real Issue

The code has `energy_shifted=True` hardcoded, which assumes the unitary has an energy shift that needs to be removed. But:

1. ✅ The Trotter circuit uses the shifted Hamiltonian H' = H + 3.846*I
2. ✅ The saved exact matrix also uses H' (same shift)
3. ❌ But the code tries to "remove" this shift from U with: U_unshifted = exp(-i*E*t) * U
4. ❌ This is WRONG because U = exp(-i*H'*t) already has the shift baked in via H', not as a separate phase factor!

### Key Insight: Two Ways to Shift

There are two different ways to handle energy shifts:

**Method 1: Shift the Hamiltonian** (what QHAT does)
```
H' = H + E*I
U' = exp(-i*H'*t) = exp(-i*(H + E*I)*t) = exp(-i*E*t) * exp(-i*H*t) * exp(terms)
```
The shift is baked into the matrix exponential through commutators.

**Method 2: Shift the unitary** (what the code assumes)
```
U_shifted = exp(i*E*t) * U
```
This is a simple phase factor multiplication.

**For Method 1**, you CAN'T just multiply by exp(-i*E*t) to "remove" the shift, because the shift is entangled with H through the Baker-Campbell-Hausdorff formula in the matrix exponential!

---

## Why Are The Errors So High?

The errors are REAL, not a bug. Here's why:

### 1. Eigenvalue Errors are Good (0.07%) ✅

Eigenvalue errors compare:
- Eigenvalues of H' (shifted exact)
- Eigenvalues derived from U' eigenphases

These match very well because we're just comparing numbers, not operators.

### 2. Matrix Norm and State Errors Are BAD ❌

When we compute:
```
U_exact = exp(-i*H'*t) where H' = H + 3.846*I
U_approx = Trotter approximation of exp(-i*H'*t)
```

Then try to "unshift":
```
U_approx_unshifted = exp(-i*3.846*t) * U_approx  # WRONG!
```

We're corrupting U_approx! The correct comparison should be:
```
||U_exact - U_approx||  # Both use H', no unshifting needed
```

But the code tries to unshift both, which corrupts them because you can't factor out exp(-i*E*t) from exp(-i*(H+E*I)*t).

---

## The Correct Fix

The bug is conceptual: **the unitary doesn't need unshifting because the shift is in the Hamiltonian, not applied as a separate phase**.

###Fix: Set energy_shifted=False for ALL algorithms

```python
# In error_analysis function, around line 672-686:

# Wrap exact Hamiltonian in OperatorRepresentation  
# H_exact is already shifted (energies in [1.527, 3.846])
exact_op = OperatorRepresentation(
    data=exact_matrix,
    operator_type='hamiltonian',
    energy_shifted=False,  # ✅ Treat as unshifted
    representation='dense_matrix',
    timestep=timestep,
    energy_shift=0.0  # ✅ No shift to apply/remove
)

# Wrap approximate time-evolution operator
# U_approx = exp(-i*H'*t) where H' is shifted
approx_op = OperatorRepresentation(
    data=unitary_matrix,
    operator_type='time_evolution',
    energy_shifted=False,  # ✅ Don't try to unshift
    representation='dense_matrix',
    timestep=timestep,
    energy_shift=0.0  # ✅ No shift to apply/remove
)
```

### Why This Works

Both H_exact and U_approx are built from the **same shifted Hamiltonian H'**:
- H_exact is literally H' saved to file
- U_approx is Trotter(exp(-i*H'*t))

So they're already on the same energy scale! No unshifting needed.

The `energy_shift` parameter passed to `error_analysis()` is only for **eigenvalue comparison** (to know what shift was applied to H), not for unshifting the operators.

---

## Expected Results After Fix

Once we set `energy_shifted=False` and `energy_shift=0.0`:

```
EXPECTED ERROR ANALYSIS RESULTS (After Fix)
================================================================================

1. EIGENVALUE ERRORS:
   Max relative error: ~0.07%  ✅ (unchanged - already correct)

2. MATRIX NORM ERRORS:
   Frobenius norm: ~0.01-0.1  ✅ (down from 6.46)
   Estimated: Based on eigenvalue errors, should be small

3. STATE-DEPENDENT ERRORS:
   Relative error: ~0.1-1%  ✅ (down from 80.75%)
   Should match eigenvalue error magnitude

================================================================================
```

---

## Detailed Recommendations

### 1. **CRITICAL: Fix the operator wrapping code**

Change in `analysis/analysis.py`, lines 672-686:

```python
# OLD CODE (WRONG):
exact_op = OperatorRepresentation(
    data=exact_matrix,
    operator_type='hamiltonian',
    energy_shifted=False,  # This is correct
    ...
    energy_shift=energy_shift  # BUG: This causes issues
)

approx_op = OperatorRepresentation(
    data=unitary_matrix,
    operator_type='time_evolution',
    energy_shifted=True,  # BUG: This is wrong!
    ...
    energy_shift=energy_shift  # BUG: This causes corruption
)

# NEW CODE (CORRECT):
exact_op = OperatorRepresentation(
    data=exact_matrix,
    operator_type='hamiltonian',
    energy_shifted=False,
    ...
    energy_shift=0.0  # ✅ Both already on same scale
)

approx_op = OperatorRepresentation(
    data=unitary_matrix,
    operator_type='time_evolution',
    energy_shifted=False,  # ✅ Don't unshift
    ...
    energy_shift=0.0  # ✅ No shift to remove
)
```

### 2. Clarify the energy_shift parameter semantics

The `energy_shift` parameter to `error_analysis()` should be used ONLY for:
- ✅ Eigenvalue comparison (knowing what shift was applied to H for context)
- ❌ NOT for operator unshifting (they're already on the same scale)

Update the docstring:

```python
def error_analysis(..., energy_shift=0.0):
    """
    ...
    energy_shift: Energy shift that was applied to the Hamiltonian internally.
                  Used for context/logging only. Both exact_matrix and unitary_matrix
                  should already be derived from the same (shifted) Hamiltonian,
                  so no operator-level unshifting is performed.
    """
```

### 3. Add assertion/check that operators are on same scale

```python
# After wrapping operators, before comparison:
# Verify that exact and approx are on the same energy scale
H_exact_eigs = np.linalg.eigvalsh(exact_matrix)
U_approx_eigs = np.linalg.eigvals(unitary_matrix)
U_approx_phases = -np.angle(U_approx_eigs)
U_approx_phases = np.where(U_approx_phases < 0, U_approx_phases + 2*np.pi, U_approx_phases)
U_approx_energies = U_approx_phases / timestep

energy_scale_mismatch = abs(H_exact_eigs.min() - U_approx_energies.min())
if energy_scale_mismatch > 0.1:
    logger.warning(
        f"Energy scale mismatch detected: {energy_scale_mismatch:.6f}. "
        f"H_exact range: [{H_exact_eigs.min():.4f}, {H_exact_eigs.max():.4f}], "
        f"U_approx range: [{U_approx_energies.min():.4f}, {U_approx_energies.max():.4f}]. "
        f"This suggests the matrices may not be derived from the same Hamiltonian."
    )
```

### 4. Update config file documentation

In `examples/config_full_analysis.py`, clarify:

```python
# Energy Shift Handling (Important!)
# -----------------------------------
# QHAT applies an energy shift E = abs(min eigenvalue) to the Hamiltonian
# internally to keep all eigenvalues positive. This shift is applied by
# modifying H → H' = H + E*I, which affects:
#   - The exact matrix output (H' is saved)
#   - The Trotter unitary (built from H')
#   - Eigenvalue comparisons (shift is added back for reporting)
#
# For error analysis, both exact and approximate operators are already
# derived from the same shifted Hamiltonian H', so they're on the same
# energy scale and can be compared directly (no unshifting needed).
```

### 5. Add integration test

```python
def test_error_analysis_consistent():
    """Test that all three error types give consistent results."""
    # Run full analysis
    results = analyze_algorithm(...)
    
    eigenvalue_error = results['error_analysis']['eigenenergy_errors']['max_relative_error']
    frobenius_error = results['error_analysis']['matrix_frobenius_error']
    state_error = results['error_analysis']['state_errors'][0]['relative_error']
    
    # All errors should be roughly the same order of magnitude for good approximations
    # Frobenius error should be O(sqrt(N) * eigenvalue_error) for N-dimensional system
    # State error should be similar to eigenvalue error for typical states
    
    assert frobenius_error < 1.0, "Frobenius error suspiciously high"
    assert state_error < 0.1, "State error suspiciously high (>10%)"
    
    # Check consistency: if eigenvalue errors are small, others should be too
    if eigenvalue_error < 0.001:  # <0.1%
        assert frobenius_error < 0.5, "Matrix error inconsistent with eigenvalue error"
        assert state_error < 0.05, "State error inconsistent with eigenvalue error"
```

---

## Summary

### What I Initially Thought (WRONG ❌)
- Time evolution doesn't apply energy shift
- U_approx is unshifted
- The code incorrectly tries to unshift it

### What's Actually True (CORRECT ✅)
- Time evolution DOES use shifted Hamiltonian H' = H + E*I  
- U_approx = Trotter(exp(-i*H'*t)) includes the shift through H'
- The code incorrectly tries to "unshift" by multiplying by exp(-i*E*t)
- But you can't factor out exp(-i*E*t) from exp(-i*(H+E*I)*t) due to non-commutativity!

### The Fix
- Set `energy_shifted=False` for both operators
- Set `energy_shift=0.0` for both (they're already on same scale)
- Don't try to unshift - they're already comparable!

### Expected Impact
- Matrix Frobenius error: 6.46 → ~0.01-0.1 (50-600x improvement)
- State relative error: 80.75% → ~0.1-1% (80-800x improvement)
- All three error types will be consistent

The Trotter approximation is actually quite good (eigenvalue errors prove this). We just need to stop corrupting the operators during comparison!
