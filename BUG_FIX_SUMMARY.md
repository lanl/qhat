# Bug Fix Summary: Error Analysis Corrections

**Date**: 2026-07-21  
**Branch**: bkk_phase2_operator_framework  
**Status**: ✅ FIXED

---

## Overview

Fixed two critical bugs in the error analysis framework that were causing incorrect error measurements. After fixes, all three error types (eigenvalue, matrix norm, state-dependent) are now consistent and physically meaningful.

---

## Bug #1: Sign Convention Mismatch

### Problem
Opposite sign conventions between driver.py and OperatorRepresentation for energy shifts:

| Component | Convention | Formula |
|-----------|-----------|---------|
| Driver (hamiltonian.py) | `H' = H + E·I` | Adds positive shift |
| OperatorRepresentation (old) | `H_shifted = H - E·I` | Subtracts shift |

This caused incorrect operator conversions when `energy_shifted=True`.

### Solution
Changed OperatorRepresentation to match driver convention:
- `_apply_energy_shift`: Now adds E to Hamiltonian eigenvalues (was: subtract)
- `_remove_energy_shift`: Now subtracts E (was: add)
- Phase factors: Now use exp(-i·E·t) for shifting (was: exp(+i·E·t))

### Files Modified
- `analysis/operators.py`: Fixed sign in shift/unshift methods
- `analysis/tests/test_operators.py`: Updated 7 tests to match new convention
- `analysis/analysis.py`: Updated to use `energy_shifted=True` correctly

### Commit
`ff76979` - "Fix sign convention mismatch in OperatorRepresentation"

---

## Bug #2: Eigenvector Reconstruction Error

### Problem  
`_get_eigendecomposition()` was using `scipy.linalg.eig()` for shifted Hamiltonians instead of `scipy.linalg.eigh()`.

**Why this matters:**
- `eig()`: General eigenvalue solver, returns eigenvectors that are NOT properly orthonormal for Hermitian matrices
- `eigh()`: Hermitian-specific solver, returns properly orthonormal eigenvectors

**The corruption:**
When modifying eigenvalues and reconstructing via `V @ diag(λ) @ V†`:
- With `eigh()` eigenvectors: Reconstruction is exact ✓
- With `eig()` eigenvectors: Reconstruction has large errors ✗

**Example:**
```
Input: H' with eigenvalues [1.527, 3.846] (shifted)
Goal: Get H with eigenvalues [-2.319, 0] (unshifted)

With eigh() eigenvectors:
  Subtract E=3.846 → [-2.319, 0]
  Reconstruct → matrix with eigenvalues [-2.319, 0] ✓

With eig() eigenvectors:
  Subtract E=3.846 → [-2.319, 0]  
  Reconstruct → matrix with eigenvalues [-3.489, 0] ✗ WRONG!
```

### Solution
Changed condition in `_get_eigendecomposition()`:

**Before (buggy):**
```python
if operator_type == 'hamiltonian' and not energy_shifted:
    eigenvalues, eigenvectors = scipy.linalg.eigh(data)
else:
    eigenvalues, eigenvectors = scipy.linalg.eig(data)
```

**After (fixed):**
```python
if operator_type == 'hamiltonian':
    # Use eigh for ALL Hamiltonians (shifted or unshifted)
    eigenvalues, eigenvectors = scipy.linalg.eigh(data)
else:
    # Time-evolution operators use eig
    eigenvalues, eigenvectors = scipy.linalg.eig(data)
```

### Files Modified
- `analysis/operators.py`: Fixed eigendecomposition method selection

### Diagnostic Test Added
Created `test_operator_conversions.py` - comprehensive test that:
- Tests all four operator forms (H/U × shifted/unshifted)
- Compares OperatorRepresentation vs direct computation  
- Checks unitarity and eigenvalue consistency
- **Helped identify the eigenvector reconstruction bug**

### Commits
`7f0e92c` - "Fix eigenvector reconstruction bug in OperatorRepresentation"
`cb42777` - "Fix matrix reconstruction to use proper inverse for time-evolution operators"

---

## Results: Before vs After

### config_full_analysis.py (Be-H molecule, 6 qubits, Trotter)

| Error Type | Before Fixes | After Fixes | Improvement |
|------------|-------------|-------------|-------------|
| **Eigenvalue (max rel)** | 0.067% | 0.067% | No change (was correct) |
| **Matrix Frobenius** | 6.46 | 0.038 | **170x better** ✅ |
| **State relative** | 80.75% | 0.23% | **350x better** ✅ |

### Physical Interpretation
All three error types now agree: **the Trotter approximation is excellent**, with errors under 0.3%.

The previous high errors (6.46 and 80.75%) were artifacts of the bugs, not real approximation errors.

---

## Technical Details

### Energy Shift Flow

**Physical System:**
- Original Hamiltonian H: eigenvalues can be negative (e.g., [-2.32, 0])
- For phase estimation, need positive eigenvalues

**Driver applies shift:**
```python
# driver.py
E = abs(min_eigenvalue)  # E = 3.846
physical_hamiltonian.energy_shift(E)  # H → H' = H + 3.846·I
```

**Result:**
- H' (shifted): eigenvalues [1.53, 3.85] (all positive)
- Saved to file as "exact_hamiltonian.npz"
- U' = exp(-i·H'·t) computed via Trotter

**For error analysis:**
- Both H' and U' are on the shifted scale
- To compare on physical scale, must unshift back to H
- OperatorRepresentation handles this conversion

### Operator Conversion Path

With the fixes, the conversion works correctly:

1. **Input**: H' (dense matrix) with eigenvalues [1.53, 3.85]
2. **Eigendecompose**: Use `eigh()` to get λ=[1.53, 3.85], V (orthonormal)
3. **Unshift eigenvalues**: λ_physical = λ - E = [-2.32, 0]
4. **Reconstruct**: H_physical = V @ diag([-2.32, 0]) @ V† ✓
5. **Convert to U**: U_physical = exp(-i·H_physical·t)

For U' → U_physical:
1. **Input**: U' with eigenphases encoding energies [1.53, 3.85]
2. **Unshift**: Multiply by exp(+i·E·t) phase factor
3. **Result**: U_physical with eigenphases encoding [-2.32, 0] ✓

---

## Testing

### Unit Tests
- All 298 tests pass
- Updated 7 tests in `test_operators.py` to match new sign convention
- All tests verify the fixes work correctly

### Integration Test
- `test_operator_conversions.py`: Comprehensive diagnostic  
- Tests all four operator forms
- Verifies consistency between OperatorRepresentation and direct computation
- Confirms Frobenius error = 0.038 ✓

### Real-World Test  
- `config_full_analysis.py`: Be-H molecule analysis
- All three error types now consistent
- Physically meaningful results

---

## Lessons Learned

1. **Convention matters**: Even with correct math, opposite conventions cause bugs
2. **Use the right tool**: `eigh()` for Hermitian, not `eig()`
3. **Eigenvector quality matters**: Reconstruction errors compound
4. **Comprehensive testing**: The diagnostic test was crucial for finding Bug #2
5. **Verify physically**: Error measures should be self-consistent

---

## Files Changed

### Core Fixes
- `analysis/operators.py`: Fixed sign convention + eigenvector bug (2 commits)
- `analysis/analysis.py`: Updated to use corrected convention (1 commit)
- `analysis/tests/test_operators.py`: Updated tests (1 commit)

### Documentation
- `SIGN_CONVENTION_MISMATCH.md`: Documents the sign convention issue
- `BUG_FIX_SUMMARY.md`: This file
- `test_operator_conversions.py`: Comprehensive diagnostic test

### Commits
1. `2320ae7` - Initial (incorrect) fix attempt
2. `ff76979` - Fix sign convention mismatch ✅
3. `7f0e92c` - Fix eigenvector reconstruction bug (eigh vs eig) ✅
4. `cb42777` - Fix matrix reconstruction (V† vs V^(-1)) ✅
5. `95e9927` - Add documentation and diagnostic test ✅

---

## Status

✅ **Both bugs fixed**  
✅ **All tests pass**  
✅ **Error measures consistent and physically meaningful**  
✅ **Ready for integration**

The error analysis framework now correctly compares operators on the physical energy scale, producing accurate and consistent error measurements.
