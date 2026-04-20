# CRITICAL FINDING: Formula Issues in qre_unitary.py

## Summary

✅ **jkg_utils.py is correct** - C1 and C2 match the paper exactly

❌ **qre_unitary.py has wrong formulas** - The step count calculations don't match the paper

## Test Results

Test case: H = X₀ + Y₀ with t=1.0, E_err=0.1

### Second-Order: ✅ CORRECT at t=1.0
- Code formula: s2 = 4.58 steps
- Paper formula: r = 4.58 steps  
- **Match!** ✓

### First-Order: ❌ WRONG
- Code formula: s1 = 125.66 steps
- Paper formula: r = 62.83 steps
- **Off by factor of 2!** ❌

## The Smoking Gun: Time Dependence

Testing different evolution times reveals the dimensional problem:

| Time (t) | s1 (code) | r1 (paper) | Ratio | s2 (code) | r2 (paper) | Ratio |
|----------|-----------|------------|-------|-----------|------------|-------|
| 0.5      | 125.66    | 31.42      | 4.00  | 3.24      | 2.29       | 1.41  |
| 1.0      | 125.66    | 62.83      | 2.00  | 4.58      | 4.58       | 1.00  |
| 2.0      | 125.66    | 125.66     | 1.00  | 6.47      | 9.15       | 0.71  |
| 4.0      | 125.66    | 251.33     | 0.50  | 9.15      | 18.31      | 0.50  |

**Key observations:**

1. **s1 is CONSTANT** (always 125.66) - doesn't depend on time at all! ❌
2. **s2 scales as √t** but should scale as t^(3/2) ❌
3. The ratios change with time, proving dimensional inconsistency

## Root Cause Analysis

### First-Order Formula in Code
```python
s1 = timestep * error_scale * c1 / eps_trotter
   = t * c1 / (E_err * t / 2π)
   = c1 * 2π / E_err
```
**Result: Independent of t** ❌

**Expected from paper:** r = t² * c1 / (2ε)

### Second-Order Formula in Code
```python
s2 = timestep * sqrt(error_scale * c2 / eps_trotter)
   = t * sqrt(c2 / (E_err * t / 2π))
   = sqrt(t) * sqrt(c2 * 2π / E_err)
```
**Result: Scales as √t** ❌

**Expected from paper:** r = t^(3/2) * sqrt(c2 / ε)

## What This Means

### If You're Using First-Order Trotter:
- The code will give the **same** number of steps regardless of evolution time
- This is **physically wrong** - longer times need more steps
- The formula needs a factor of **t²** 

### If You're Using Second-Order Trotter:
- The code gives the **right answer at t=1.0** (by accident?)
- For other times, it's wrong by a factor of **√t**
- The formula needs an extra factor of **t**

## Possible Corrections

### Option 1: Fix the formulas in qre_unitary.py

```python
# Current (WRONG):
s1 = timestep * error_scale * c1 / eps_trotter
s2 = timestep * math.sqrt(error_scale * c2 / eps_trotter)

# Should be:
s1 = timestep**2 * error_scale * c1 / (2 * eps_trotter)
s2 = timestep**2 * math.sqrt(error_scale * c2 * timestep / eps_trotter)
```

### Option 2: Fix the definition of eps_trotter

Perhaps eps_trotter should be defined differently to make the formulas work:

```python
# Current:
eps_trotter = energy_error * timestep / (2 * math.pi)

# Alternative 1 (for first-order):
eps_trotter = energy_error * timestep**2 / (2 * math.pi)

# Alternative 2 (adjust both):
# Would need different eps values for first and second order
```

## Immediate Actions Required

1. **DO NOT USE FIRST-ORDER TROTTER** with this code - the formula is wrong

2. **Second-order is OK ONLY if timestep ≈ 1.0** in your units

3. **Contact the code author URGENTLY** to understand:
   - Was this bug known?
   - Are there compensating factors elsewhere?
   - What were the intended formulas?

4. **Check all previous results** that used first-order or non-unit timesteps

## Questions for Code Author

1. Why is `timestep` in the numerator instead of squared?
2. What is the physical meaning of `eps_trotter`?
3. Why does second-order happen to work at t=1.0?
4. Are there any published results using this code that need correction?
5. What is the relationship between `energy_error` and operator norm error?

## Verification

This test used:
- Analytical calculation of C1=2.0, C2=1/3 (verified by Monte Carlo)
- Direct application of Childs et al. Equations 145 and 152
- Multiple evolution times to test scaling

The discrepancy is **not** due to:
- Incorrect C1 or C2 values (verified ✓)
- Different error definitions (tested multiple interpretations)
- Numerical precision (differences are large and systematic)

## Recommendation

**STOP using qre_unitary.py until the formulas are corrected or explained.**

The second-order formula happens to give reasonable results near t=1.0, but this appears accidental. The first-order formula is clearly wrong.
