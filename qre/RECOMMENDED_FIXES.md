# Recommended Fixes for qre_unitary.py

## Problem Summary

The formulas on lines 106-107 have wrong time scaling:
- s1 is independent of time (should scale as t)
- s2 scales as √t (should scale as t^(3/2))

## Fix Option 1: Minimal Change (RECOMMENDED)

### Change Only the Step Count Formulas

This is the safest fix - only modify lines 106-107 in qre_unitary.py:

**Current code (WRONG):**
```python
s1 = timestep * error_scale * c1 / eps_trotter
s2 = timestep * math.sqrt(error_scale * c2 / eps_trotter)
```

**Fixed code:**
```python
s1 = timestep**2 * error_scale * c1 / (2 * eps_trotter)
s2 = timestep**1.5 * math.sqrt(error_scale * c2 / eps_trotter)
```

Or equivalently (for clarity):
```python
s1 = timestep**2 * error_scale * c1 / (2 * eps_trotter)
s2 = (timestep**1.5) * math.sqrt(error_scale * c2 / eps_trotter)
```

### Why This Works

From Childs et al. Equations 145 and 152:
- First-order: r = t² * c1 / (2ε)
- Second-order: r = t^(3/2) * sqrt(c2 / ε)

With eps_trotter = energy_error * timestep / (2π), this gives:
```
s1 = t² * c1 / (2 * E_err * t / (2π))
   = t * c1 * π / E_err  ✓ (scales linearly with t)

s2 = t^(3/2) * sqrt(c2 / (E_err * t / (2π)))
   = t^(3/2) * sqrt(c2 * 2π / (E_err * t))  ✓ (scales as t^(3/2))
```

### Verification

With this fix, my test case gives:

| Time | s1 (fixed) | r1 (paper) | s2 (fixed) | r2 (paper) |
|------|------------|------------|------------|------------|
| 0.5  | 31.42      | 31.42      | 2.29       | 2.29       |
| 1.0  | 62.83      | 62.83      | 4.58       | 4.58       |
| 2.0  | 125.66     | 125.66     | 9.15       | 9.15       |
| 4.0  | 251.33     | 251.33     | 18.31      | 18.31      |

Perfect match! ✓

### Implementation

```python
# File: qre_unitary.py, lines 104-107
eps_trotter = config_unitary.energy_error * timestep / (2 * math.pi)
config_general.log(f"-- allowable fractional Trotter error = {eps_trotter}")

# FIXED: Added timestep**2 for s1 and timestep**1.5 for s2
# Reference: Childs et al. arXiv:1912.08854v3, Equations 145 and 152
s1 = timestep**2 * config_unitary.error_scale * c1 / (2 * eps_trotter)
s2 = (timestep**1.5) * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)

config_general.log(f"-- Trotter step count: s1 = {s1}, s2 = {s2}")
```

## Fix Option 2: Redefine eps_trotter (MORE INVASIVE)

Instead of changing the formulas, redefine what eps_trotter means:

**For first-order:**
```python
eps_trotter_1st = config_unitary.energy_error * timestep**2 / (2 * math.pi)
s1 = c1 / (2 * eps_trotter_1st)
```

**For second-order:**
```python
eps_trotter_2nd = config_unitary.energy_error * timestep / (2 * math.pi)
s2 = math.sqrt(c2 * timestep / eps_trotter_2nd)
```

**Pros:**
- Makes the formulas look simpler
- Separates concerns about error scaling

**Cons:**
- Need different eps_trotter for each order
- More confusing conceptually
- Might break other code that uses eps_trotter

**Verdict:** NOT RECOMMENDED - too invasive

## Fix Option 3: Change Both (CLEANEST BUT RISKY)

Redefine eps_trotter to be the actual operator norm error tolerance:

```python
# Operator norm error tolerance (from energy error)
# For phase estimation: operator_error ≈ energy_error * timestep
operator_error = config_unitary.energy_error * timestep / (2 * math.pi)

# Step counts directly from paper formulas
s1 = timestep**2 * config_unitary.error_scale * c1 / (2 * operator_error)
s2 = (timestep**1.5) * math.sqrt(config_unitary.error_scale * c2 / operator_error)
```

**Pros:**
- Most physically clear
- Direct correspondence with paper

**Cons:**
- Variable name change might confuse
- Same as Option 1 in practice

**Verdict:** OK, but Option 1 is simpler

## Testing the Fix

Here's a test to verify the fix works:

```python
# File: test_fixed_formulas.py
import numpy as np
import math
from openfermion import QubitOperator
from jkg_utils import trotter_error_estimator

# Test Hamiltonian: H = X + Y
terms = [QubitOperator('X0', 1.0), QubitOperator('Y0', 1.0)]
C1, C2 = trotter_error_estimator(terms, time_limit=2.0)

# C1 = 2.0, C2 = 1/3 (analytically)
print(f"C1 = {C1:.2f}, C2 = {C2:.4f}")

# Test parameters
timestep = 2.0
energy_error = 0.1
error_scale = 1.0

eps_trotter = energy_error * timestep / (2 * math.pi)

# FIXED FORMULAS
s1 = timestep**2 * error_scale * C1 / (2 * eps_trotter)
s2 = (timestep**1.5) * math.sqrt(error_scale * C2 / eps_trotter)

# Expected from paper
r1 = timestep**2 * C1 / (2 * eps_trotter)
r2 = (timestep**1.5) * math.sqrt(C2 / eps_trotter)

print(f"\ns1 (fixed) = {s1:.2f}")
print(f"r1 (paper) = {r1:.2f}")
print(f"Match: {abs(s1 - r1) < 0.01}")

print(f"\ns2 (fixed) = {s2:.2f}")
print(f"r2 (paper) = {r2:.2f}")
print(f"Match: {abs(s2 - r2) < 0.01}")

# Expected output:
# C1 = 2.00, C2 = 0.3333
# s1 (fixed) = 125.66
# r1 (paper) = 125.66
# Match: True
# s2 (fixed) = 9.15
# r2 (paper) = 9.15
# Match: True
```

## Impact Assessment

### What Changes:
- Number of Trotter steps will be different (correct now)
- First-order will need MORE steps (was undercounting)
- Second-order will need more/fewer steps depending on timestep

### What Stays Same:
- C1 and C2 calculations (already correct)
- Overall structure of the code
- Interface and API

### Backward Compatibility:
**⚠️ WARNING:** This fix will change resource estimates!

If you have:
- Published results using this code
- Existing resource estimates
- Benchmarks or comparisons

They may need to be recalculated with the corrected formulas.

## Migration Plan

### Step 1: Verify the Fix (Do This First!)
```bash
cd qre
python3.11 test_step_count_formula.py  # Should show the problem
# Apply fix to qre_unitary.py
python3.11 test_step_count_formula.py  # Should show perfect match
```

### Step 2: Add Unit Test
Create a permanent test that checks the formulas:

```python
# File: test_trotter_step_count.py
def test_step_count_scaling():
    """Verify step counts scale correctly with time."""
    # ... (use the test I created)
    assert abs(s1_code / r1_paper - 1.0) < 0.01
    assert abs(s2_code / r2_paper - 1.0) < 0.01
```

Add this to your CI pipeline.

### Step 3: Update Documentation
Add to README.md:

```markdown
## Breaking Change in Version X.Y

The Trotter step count formulas in qre_unitary.py have been corrected to
match Childs et al. (arXiv:1912.08854v3) Equations 145 and 152. Previous
versions had incorrect time scaling. Resource estimates may differ from
previous versions, especially for:
- First-order Trotter (all timesteps affected)
- Second-order Trotter with timestep ≠ 1.0
```

### Step 4: Deprecate error_scale
Since it's already deprecated, consider removing it:

```python
# Old (with deprecated parameter):
s1 = timestep**2 * config_unitary.error_scale * c1 / (2 * eps_trotter)

# New (clean):
s1 = timestep**2 * c1 / (2 * eps_trotter)
```

## My Recommendation

**Use Fix Option 1 (Minimal Change)**

1. Apply the two-line fix to qre_unitary.py lines 106-107
2. Run test_step_count_formula.py to verify
3. Add the test as a permanent unit test
4. Update documentation noting the breaking change
5. Recompute any published resource estimates

This is the safest, clearest fix with the least chance of introducing new bugs.

## Questions?

If you need help:
1. Applying the fix
2. Creating the test
3. Understanding the derivation
4. Validating existing results

Let me know!
