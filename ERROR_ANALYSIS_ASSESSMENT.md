# Error Analysis Assessment for config_full_analysis.py

**Date**: 2026-07-21  
**Config**: `analysis/examples/config_full_analysis.py`  
**System**: Be-H molecule, 6 qubits, Trotter approximation

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
   Frobenius norm ||U_exact - U_approx||_F: 6.463412  ❌ VERY BAD

3. STATE-DEPENDENT ERRORS:
   Absolute error: 8.074505e-01  ❌ VERY BAD
   Relative error: 0.807450 (80.75%)  ❌ CATASTROPHIC

================================================================================
```

##  Root Cause Analysis

### The Bug

The error analysis code has a **critical bug** in how it handles energy shifts:

**File**: `analysis/analysis.py`, line 683

```python
# WRONG CODE (current)
approx_op = OperatorRepresentation(
    data=unitary_matrix,
    operator_type='time_evolution',
    energy_shifted=True,  # ❌ HARDCODED - ASSUMES ENERGY SHIFT
    representation='dense_matrix',
    timestep=timestep,
    energy_shift=energy_shift
)
```

### Why It's Wrong

For `algorithm.method = "time evolution"`:
- The unitary matrix U_approx = exp(-i*H*t) uses the **actual Hamiltonian H** (not shifted)
- Energy shift is only used internally by the Hamiltonian for bound computations
- The produced unitary is **NOT energy-shifted**

When we set `energy_shifted=True`, the code tries to "remove" a shift that doesn't exist:
```python
# What the code does:
phase_factor = exp(-i * 3.8456 * 1.4176) = 0.674 + 0.739i
U_approx_unshifted = phase_factor * U_approx  # Multiplies by wrong phase!
```

This corrupts the unitary and produces nonsensical errors.

### Evidence

1. **U_approx eigenvalues** are in range [1.53, 3.85] (same as H_exact eigenvalues)
   - If it were shifted, they'd be in range [-2.32, 0.00]
   - This proves U_approx is already unshifted!

2. **Unitarity check after "removing shift"**:
   - Before: ||U†U - I||_F = 1.04e-14 ✅ Perfect
   - After: ||U†U - I||_F = 8.50e-04 ❌ Destroyed!
   - The phase multiplication destroyed unitarity

3. **Eigenvalue errors are good** (0.07% max relative error)
   - The eigendecomposition path doesn't corrupt the energies
   - This shows the approximation is actually decent

### Algorithm-Specific Behavior

| Algorithm Type | Energy Shift in Unitary? | Correct `energy_shifted` |
|----------------|--------------------------|--------------------------|
| "time evolution" | ❌ No | `False` |
| "QPE: qualtran textbook" | ❌ No (QPE shifts phases after) | `False` |
| "QPE: pyLIQTR qubitized" | ✅ Maybe (depends on implementation) | Need to check |

The energy shift is a property of how the algorithm **encodes the Hamiltonian**, not a universal property.

---

## Recommended Fixes

### Fix #1: Determine energy_shifted from algorithm type (Quick Fix)

```python
# In error_analysis function, around line 683
# Determine if unitary is energy-shifted based on algorithm
# For now, assume only time evolution produces unshifted unitaries
algorithm_method = getattr(algorithm, 'method', 'unknown') if hasattr(algorithm, 'method') else 'unknown'
is_energy_shifted = (algorithm_method != 'time evolution')

approx_op = OperatorRepresentation(
    data=unitary_matrix,
    operator_type='time_evolution',
    energy_shifted=is_energy_shifted,  # ✅ Now determined from algorithm
    representation='dense_matrix',
    timestep=timestep,
    energy_shift=energy_shift
)
```

### Fix #2: Pass energy_shifted as parameter (Robust Fix)

Modify `analyze_algorithm()` to accept an `is_unitary_shifted` parameter:

```python
def analyze_algorithm(
        config_analysis: AnalysisConfiguration,
        algorithm,
        hamiltonian=None,
        timestep=None,
        energy_shift=0.0,
        is_unitary_shifted=False) -> dict:  # NEW PARAMETER
    ...
    results["error_analysis"] = error_analysis(
        ...,
        is_unitary_shifted=is_unitary_shifted  # Pass through
    )
```

Then in driver.py:

```python
# Determine if the algorithm produces an energy-shifted unitary
# For time evolution: no shift
# For QPE algorithms: depends on implementation (check algorithm documentation)
is_unitary_shifted = (algorithm_method != 'time evolution')

state.store_results(analyze_algorithm(
    state.config_analysis,
    algorithm,
    hamiltonian=physical_hamiltonian,
    timestep=timestep,
    energy_shift=energy_shift,
    is_unitary_shifted=is_unitary_shifted  # NEW
))
```

### Fix #3: Auto-detect from eigenvalues (Most Robust)

```python
# In error_analysis, before wrapping operators:
# Auto-detect if unitary is shifted by checking its eigenvalue range
U_eigs = np.linalg.eigvals(unitary_matrix)
phases = -np.angle(U_eigs)
phases = np.where(phases < 0, phases + 2*np.pi, phases)
energies_from_U = phases / timestep

# Get expected energy range from exact Hamiltonian
H_eigs = np.linalg.eigvalsh(exact_matrix)
H_min, H_max = H_eigs.min(), H_eigs.max()

# If U energies match H energies, it's unshifted
# If U energies match (H - E_shift) energies, it's shifted
is_unshifted = (
    abs(energies_from_U.min() - H_min) < 0.1 and
    abs(energies_from_U.max() - H_max) < 0.1
)
is_shifted = (
    abs(energies_from_U.min() - (H_min - energy_shift)) < 0.1 and
    abs(energies_from_U.max() - (H_max - energy_shift)) < 0.1
)

if is_unshifted:
    logger.info("Auto-detected: unitary is NOT energy-shifted")
    is_energy_shifted = False
elif is_shifted:
    logger.info("Auto-detected: unitary IS energy-shifted")
    is_energy_shifted = True
else:
    logger.warning("Cannot auto-detect energy shift status, assuming unshifted")
    is_energy_shifted = False
```

---

## Expected Results After Fix

Once we fix the bug (use `energy_shifted=False` for time evolution):

```
EXPECTED ERROR ANALYSIS RESULTS (After Fix)
================================================================================

1. EIGENVALUE ERRORS:
   Max relative error: ~0.07%  ✅ (unchanged - already correct)

2. MATRIX NORM ERRORS:
   Frobenius norm: ~0.01-0.1  ✅ (down from 6.46)
   Spectral norm: ~0.001-0.01  ✅ (down from 0.82)

3. STATE-DEPENDENT ERRORS:
   Relative error: ~0.1-1%  ✅ (down from 80.75%)

================================================================================
```

The eigenvalue errors are already good because that code path doesn't have the bug. The matrix norm and state errors will improve dramatically once we stop corrupting the unitary.

---

## Testing Recommendations

### 1. Unit Test for energy_shifted Detection

```python
def test_energy_shift_detection_time_evolution():
    """Test that time evolution produces unshifted unitary."""
    # Create time evolution algorithm
    # Run error analysis
    # Verify is_energy_shifted = False
    
def test_energy_shift_detection_qpe():
    """Test QPE algorithm energy shift detection."""
    # Create QPE algorithm
    # Run error analysis  
    # Verify correct is_energy_shifted value
```

### 2. Integration Test with Known System

```python
def test_error_analysis_hydrogen():
    """Test error analysis on H2 molecule with known exact solution."""
    # Use small system where we can verify all errors independently
    # Check that all three error types are consistent
    # Eigenvalue errors should roughly match matrix norm errors
```

### 3. Regression Test

```python
def test_error_values_physically_reasonable():
    """Ensure error values are in reasonable ranges."""
    # For good Trotter approximations:
    assert eigenvalue_relative_error < 0.01  # < 1%
    assert frobenius_error < 1.0  # Should be O(0.1)
    assert state_relative_error < 0.1  # < 10%
```

---

## Additional Recommendations

### 1. Improve Trotter Approximation (If Needed)

The eigenvalue errors (0.07%) show the approximation is actually quite good. The Trotter parameters seem reasonable. However, if you want even better accuracy:

```python
# In config file, reduce energy_error for tighter Trotter:
unitary.encode_ramped_trotter(
    energy_error=0.1 * energy_error,  # 10x tighter (but 10x more gates)
    trotter_implementation="flattened",
    trotter_combine_terms=True,
    ordering_method="lexicographical"
)
```

**Trade-off**: Lower error but higher gate count/depth.

### 2. Add Spectral Norm to Config

The spectral norm is more physically meaningful than Frobenius:

```python
# Current: analysis.error_matrix_norms = "frobenius"
# Better: analysis.error_matrix_norms = ["frobenius", "spectral"]
```

**Why**: Spectral norm measures worst-case error, which is what matters for quantum algorithms.

### 3. Add More Test States

```python
# Current: Single HF state
analysis.error_state_inputs = "examples/Be-H_1.30_sto-6g_as-003-003_jw.npy"

# Better: Test multiple states (ground, excited, random)
analysis.error_state_inputs = [
    "examples/Be-H_ground_state.npy",
    "examples/Be-H_excited_state.npy",
    "examples/Be-H_random_state.npy"
]
```

**Why**: Error can vary significantly across different input states.

### 4. Document Energy Shift Conventions

Add to config file comments:

```python
# IMPORTANT: Energy shift handling
# --------------------------------
# The Hamiltonian uses energy shift E = min(H_eigenvalues) internally to:
#   1. Keep shifted eigenvalues in [0, E_max - E_min] for phase estimation
#   2. Reduce required phase qubits in QPE algorithms
#
# However, for time evolution algorithms (algorithm.method = "time evolution"):
#   - The unitary U = exp(-iHt) uses the ACTUAL Hamiltonian (not shifted)
#   - Energy shift is only for internal bounds computation
#   - Error analysis should use energy_shifted = False
#
# For QPE algorithms (algorithm.method starts with "QPE"):
#   - Behavior depends on implementation (check algorithm documentation)
#   - Some QPE methods shift the Hamiltonian, others don't
```

---

## Summary

### Current Status
- ✅ Eigenvalue errors work correctly (0.07% max relative error)
- ❌ Matrix norm errors are broken (showing 6.46 instead of ~0.1)
- ❌ State errors are broken (showing 80.75% instead of ~1%)

### Root Cause
- Hardcoded assumption that unitary is energy-shifted
- For "time evolution" algorithm, unitary is NOT shifted
- Multiplying by wrong phase factor corrupts the unitary

### Priority Fixes
1. **CRITICAL**: Set `energy_shifted=False` for time evolution (Fix #1)
2. **HIGH**: Implement auto-detection or parameterize (Fix #2 or #3)
3. **MEDIUM**: Add tests to prevent regression
4. **LOW**: Improve documentation and examples

### Expected Impact
- Matrix Frobenius error: 6.46 → ~0.01-0.1 (65x improvement)
- State relative error: 80.75% → ~0.1-1% (80-800x improvement)
- All three error types will be consistent and physically meaningful

The approximation quality is actually quite good (eigenvalue errors show this). We just need to fix the bug to see it properly in the other error metrics.
