# Bug #2 Details: Two-Part Fix

Bug #2 (Eigenvector Reconstruction Error) actually had **two root causes**, fixed in commits `7f0e92c` and `cb42777`:

## Part 1: Using eig() instead of eigh() for Hamiltonians (commit 7f0e92c)

**Problem:** Used `eig()` for shifted Hamiltonians, which returns non-orthonormal eigenvectors.

**Fix:** Use `eigh()` for ALL Hamiltonians (shifted or unshifted).

**Impact:** Reduced Frobenius error from 1.41 to 0.038 (37x improvement).

## Part 2: Using V† instead of V^(-1) for time-evolution operators (commit cb42777)

**Problem:** Matrix reconstruction used `V @ diag(λ) @ V†` for ALL operators, assuming orthonormal eigenvectors.

**Issue:** For time-evolution operators, `eig()` returns eigenvectors with orthonormality error ~4.2e-4, so V† ≠ V^(-1).

**Fix:** Use correct inverse based on operator type:
- Hamiltonians: Use V† (eigenvectors from eigh are orthonormal)
- Time-evolution: Use V^(-1) (eigenvectors from eig are not orthonormal)

**Impact:**
- Before: Unitarity error ~4.2e-4, OperatorRepresentation didn't match direct computation
- After: Unitarity error ~1e-14, perfect match with direct computation

## Combined Result

All three error types now consistent and physically meaningful:
- Eigenvalue error: 0.067%
- Matrix Frobenius error: 0.038
- State relative error: 0.23%

The Trotter approximation is excellent, with errors under 0.3%.
