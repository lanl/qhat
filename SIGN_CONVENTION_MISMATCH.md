# Sign Convention Mismatch Between Driver and OperatorRepresentation

**Date**: 2026-07-21  
**Issue**: Opposite sign conventions for energy shifts
**Status**: Documented (workaround in place)

---

## The Mismatch

There are **opposite sign conventions** in two parts of the codebase:

### Driver Convention (hamiltonian.py line 234)
```python
def energy_shift(self, dE):
    t0 = self._H.constant + dE  # ADDS dE
```

**Result**: `H' = H + E·I` (shift makes eigenvalues MORE POSITIVE)

When driver.py calls `energy_shift(+3.846)`, it adds 3.846 to all eigenvalues:
- Original H: eigenvalues ≈ [-2.32, 0]
- After shift: H' with eigenvalues [1.53, 3.85]

### OperatorRepresentation Convention (operators.py line 359)
```python
# Comment in code:
# For Hamiltonians: H_shifted = H - E*I, so λ_shifted = λ - E
```

**Expectation**: `H_shifted = H - E·I` (shift makes eigenvalues LESS POSITIVE)

When OperatorRepresentation sees `energy_shifted=True` with `energy_shift=3.846`:
- Assumes input has eigenvalues 3.846 LOWER than physical H
- To "unshift": adds 3.846 back

---

## Why This Causes Problems

The driver produces: H' = H + 3.846·I (eigenvalues [1.53, 3.85])

If we tell OperatorRepresentation:
- `energy_shifted=True` (input is shifted)
- `energy_shift=3.846` (shift value)

It thinks:
- Input is H_shifted = H_physical - 3.846·I  
- To get H_physical, add 3.846
- Result: H_physical with eigenvalues [5.37, 7.69]

But actually:
- Input H' with eigenvalues [1.53, 3.85] IS the physical Hamiltonian we want!
- There's no "original" H with higher eigenvalues that we need to recover

---

## Why the Conventions Differ

### Driver's Perspective
The shift is applied to make eigenvalues positive for phase estimation:
```python
# driver.py lines 57-60
Elo1, Ehi1 = physical_hamiltonian.compute_initial_energy_bounds(...)
physical_hamiltonian.energy_shift(-1 * Elo1)  # If Elo1 = -3.846, shifts by +3.846
```

This shifts H → H' where H' has all positive eigenvalues. **H' becomes the physical system** we analyze.

### OperatorRepresentation's Perspective
The class was designed for a different use case where:
- You have a physical Hamiltonian H (the "unshifted" reference)
- You might work with a shifted version H_shifted = H - E·I (e.g., for better numerical conditioning)
- You want to convert between them

This assumes there's a canonical "unshifted" H that's the physical system.

---

## The Correct Interpretation

For our error analysis use case:
- Both H_exact and U_approx are derived from H' (shifted Hamiltonian)
- H' IS the physical system - there's no separate "original" H
- The shift was internal bookkeeping, not a transformation we need to undo

Therefore:
- Both operators are on the **same scale** (H')
- No relative adjustment is needed
- Set `energy_shifted=False` and `energy_shift=0.0` for both

---

## Why Negative Sign Doesn't Help

One might think: "If conventions are opposite, use negative energy_shift!"

Tested this: `energy_shift=-3.846` with `energy_shifted=True`

**Result**: Still wrong! The issue isn't just sign - it's that the conceptual model is different:
- Driver: shift creates the physical Hamiltonian
- OperatorRepresentation: shift is a deviation from the physical Hamiltonian

---

## Solutions

### Current Workaround (Implemented) ✅
```python
exact_op = OperatorRepresentation(
    data=exact_matrix,  # H' from file
    operator_type='hamiltonian',
    energy_shifted=False,  # Treat as baseline
    energy_shift=0.0  # No adjustment needed
)

approx_op = OperatorRepresentation(
    data=unitary_matrix,  # U(H') from Trotter
    operator_type='time_evolution',
    energy_shifted=False,  # Treat as baseline
    energy_shift=0.0  # No adjustment needed
)
```

**Result**: Frobenius error = 0.038 ✅

This works because we're saying "both operators are on the same baseline scale, don't adjust them."

### Alternative Fix (Not Implemented)
Change OperatorRepresentation convention to match driver:
```python
# In operators.py line 359, change comment and implementation:
# For Hamiltonians: H_shifted = H + E*I, so λ_shifted = λ + E

def _apply_energy_shift(self, eigenvalues, operator_type):
    if operator_type == 'hamiltonian':
        return eigenvalues + self.energy_shift  # Changed from -
```

**Pros**: Would match driver convention  
**Cons**: Breaking change, affects all users of OperatorRepresentation

### Proper Fix (Future Work)
Clarify the semantics and add explicit parameters:
```python
class OperatorRepresentation:
    def __init__(self, ..., shift_convention='lower'):
        """
        shift_convention: 'lower' or 'higher'
            'lower': shifted version has lower eigenvalues (H_shifted = H - E·I)
            'higher': shifted version has higher eigenvalues (H_shifted = H + E·I)
        """
```

---

## Testing Results

| Configuration | Frobenius Error | Interpretation |
|--------------|----------------|----------------|
| energy_shift=0.0, energy_shifted=False | 0.038 | ✅ Both on same baseline |
| energy_shift=+3.846, energy_shifted=True | 1.406 | ❌ Wrong: unshifts in wrong direction |
| energy_shift=-3.846, energy_shifted=True | 1.406 | ❌ Wrong: still incorrect |

Only the first configuration gives correct results because it avoids the sign convention mismatch entirely.

---

## Recommendation

**Document the convention mismatch** but keep the current workaround:
1. Add comments in error_analysis() explaining why energy_shift=0.0
2. Update OperatorRepresentation docstring to clarify its convention
3. Consider adding a sign_convention parameter in future versions

**Do NOT** change OperatorRepresentation's convention without careful review of all uses.

---

## Key Takeaway

The sign convention mismatch exists, but the **deeper issue** is conceptual:
- Driver treats shifted Hamiltonian as the physical system
- OperatorRepresentation treats unshifted Hamiltonian as the physical system

For error analysis, both input operators describe the **same physical system**, so they need **no relative adjustment** regardless of which convention we use. Setting energy_shift=0.0 sidesteps the convention mismatch entirely.
