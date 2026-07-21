# Energy Shift Fix for Error Analysis

**Date**: 2026-07-21  
**Issue**: Matrix norm and state-dependent errors were showing nonsensical values
**Status**: ✅ FIXED

---

## Summary

Fixed incorrect energy shift handling in error_analysis() that was corrupting operator comparisons.

### Results

| Error Type | Before Fix | After Fix | Improvement |
|------------|------------|-----------|-------------|
| Eigenvalue (max rel) | 0.067% ✅ | 0.067% ✅ | No change (was correct) |
| Matrix Frobenius | 6.46 ❌ | 0.038 ✅ | 170x improvement |
| State relative | 80.75% ❌ | 0.23% ✅ | 350x improvement |

All three error types are now consistent and physically meaningful!

---

## Root Cause

The bug was in `analysis/analysis.py` lines 668-688. The code was treating the operators as if they needed energy shift adjustments, when in fact:

1. **H_exact** is saved as H' = H + E·I (shifted Hamiltonian)
2. **U_approx** is computed as exp(-i·H'·t) (using the SAME shifted Hamiltonian)
3. Both are **already on the same energy scale**
4. No relative shift adjustment is needed!

The old code was trying to "unshift" the operators, which corrupted them because:
- Setting `energy_shifted=True` tells OperatorRepresentation that the input has a shift that needs to be removed
- But the operators don't have a **relative** shift between them
- They're both derived from H', so they're already comparable

---

## The Fix

Changed operator wrapping in `analysis/analysis.py` lines 668-688:

```python
# OLD CODE (WRONG):
exact_op = OperatorRepresentation(
    data=exact_matrix,
    operator_type='hamiltonian',
    energy_shifted=False,  # Claimed it wasn't shifted
    energy_shift=energy_shift  # But passed a shift value
)

approx_op = OperatorRepresentation(
    data=unitary_matrix,
    operator_type='time_evolution',
    energy_shifted=True,  # Tried to "unshift" it
    energy_shift=energy_shift
)
```

```python
# NEW CODE (CORRECT):
exact_op = OperatorRepresentation(
    data=exact_matrix,
    operator_type='hamiltonian',
    energy_shifted=False,  # Treat as baseline
    energy_shift=0.0  # No relative shift to apply
)

approx_op = OperatorRepresentation(
    data=unitary_matrix,
    operator_type='time_evolution',
    energy_shifted=False,  # Treat as baseline
    energy_shift=0.0  # No relative shift to apply
)
```

---

## Key Insight: Two Different "Shifts"

There are two different concepts of "energy shift" that were being confused:

### 1. **Absolute Energy Shift** (applied to Hamiltonian in driver.py)
- Purpose: Keep eigenvalues positive for phase estimation
- Applied in `driver.py` line 60: `physical_hamiltonian.energy_shift(-1 * Elo1)`
- Transforms H → H' = H + E·I
- Affects: H_exact (saved as H'), U_approx (built from H'), eigenvalue reporting

### 2. **Relative Energy Shift** (for operator comparison)
- Purpose: Adjust operators to same energy scale for comparison
- Handled by: OperatorRepresentation class in error_analysis()
- **Not needed here**: Both operators already use H', so no relative adjustment required

The bug was treating **absolute shift** as if it were a **relative shift** that needed correction.

---

## What energy_shift Parameter Means

The `energy_shift` parameter passed to `error_analysis()` should be understood as:

- ✅ **For eigenvalue reporting**: Add shift back to report physical eigenvalues
- ✅ **For documentation**: Record what shift was applied to H
- ❌ **NOT for operator unshifting**: Operators are already on the same scale

When wrapping operators in OperatorRepresentation:
- If inputs are derived from the SAME Hamiltonian (shifted or not): use `energy_shift=0.0`
- If inputs are from DIFFERENT energy scales: use `energy_shift=difference`

---

## Testing the Fix

Ran config_full_analysis.py with the fix:

**Eigenvalue errors** (unchanged - were already correct):
```
max_absolute_error: 0.003957573 (0.4%)
max_relative_error: 0.000674945 (0.067%)
```

**Matrix norm errors** (FIXED):
```
Before: Frobenius = 6.46 (nonsensical)
After:  Frobenius = 0.038 (excellent!)
```

**State-dependent errors** (FIXED):
```
Before: relative_error = 80.75% (nonsensical)
After:  relative_error = 0.23% (excellent!)
```

All three error types now agree that the Trotter approximation is very accurate!

---

## Validation

Verified the fix by comparing operators at different scales:

**Test 1: Direct comparison (same shifted scale)**
```python
||H_exact - H_exact||_F = 0  ✅
||U_approx_from_H' - U_approx_from_H'||_F = 0  ✅
||U_exact_from_H' - U_approx_from_H'||_F = 0.038  ✅ Excellent!
```

**Test 2: After attempting to unshift (WRONG)**
```python
||U_exact_unshifted - U_approx_unshifted||_F = 6.46  ❌ Bad!
```

This proves that:
1. Both operators ARE on the shifted scale (H')
2. Direct comparison gives correct results
3. Attempting to unshift corrupts the comparison

---

## Architectural Note

The `OperatorRepresentation` class is working correctly. The issue was how we were using it:

- **Class design**: Handles conversions between H/U and shifted/unshifted forms
- **Bug**: We told it the inputs had different shift states when they didn't
- **Fix**: Tell it both inputs are at the same baseline (shift=0)

The `energy_shifted` parameter means "does this operator have a shift **relative to what you want to compare**", not "was a shift applied somewhere in the pipeline".

---

## Impact on Other Algorithms

This fix applies specifically to error analysis comparing:
- Exact Hamiltonian matrix (from `hamiltonian.to_matrix()`)
- Approximate unitary matrix (from `algorithm.tensor_contract()`)

Both are derived from the same (shifted) Hamiltonian object, so no relative adjustment is needed.

For other algorithm types (QPE, etc.), the same principle applies:
- If both matrices come from the same Hamiltonian: `energy_shift=0.0`
- If they come from different energy scales: use the actual difference

---

## Files Modified

1. `analysis/analysis.py` (lines 668-688)
   - Changed `exact_op` creation: `energy_shifted=False`, `energy_shift=0.0`
   - Changed `approx_op` creation: `energy_shifted=False`, `energy_shift=0.0`
   - Updated comments to clarify the reasoning

---

## Future Work

1. **Clarify energy_shift semantics** in function docstrings
2. **Add validation check** to detect energy scale mismatches
3. **Add integration test** to ensure all three error types stay consistent
4. **Document** the two concepts of "shift" (absolute vs relative)

---

## Lessons Learned

1. **Same source = same scale**: If operators come from the same Hamiltonian, they're already comparable
2. **Verify assumptions**: The eigenvalue errors being good was a clue that the bug was in operator conversion, not approximation quality
3. **Test directly**: Direct matrix comparison revealed the bug immediately
4. **Commutativity matters**: E·I commutes with everything, so exp(-i(H+E·I)t) = exp(-i·E·t)·exp(-i·H·t), BUT that doesn't mean we should factor it out for comparison!
