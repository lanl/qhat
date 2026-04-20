# Usage Verification Status

## Summary

✅ **The C1 and C2 coefficients returned by jkg_utils.py are CORRECT** according to the paper.

⚠️ **The usage of these coefficients in qre_unitary.py needs verification** - there are dimensional issues.

## What's Correct

### jkg_utils.py and jkg_utils_fast.py
- ✅ C1 = Σᵢ<ⱼ ||[Hᵢ, Hⱼ]|| matches Equation 145
- ✅ C2 = C21/12 + C22/24 matches Equation 152
- ✅ Coefficients 1/12 and 1/24 are correct
- ✅ All tests pass

## What Needs Review

### qre_unitary.py lines 104-107

```python
eps_trotter = config_unitary.energy_error * timestep / (2 * math.pi)
s1 = timestep * config_unitary.error_scale * c1 / eps_trotter
s2 = timestep * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)
```

### Issue 1: Dimensional Analysis Problem

**First-order formula (s1):**
- Expected from paper: `r = t² * c1 / (2ε)`
- Code computes: `s1 = t * c1 / eps_trotter`
- **Analysis:** Dimensions work out IF `eps_trotter ∝ t`, which it is (`E_err * t / 2π`)
- **Conclusion:** Dimensionally consistent ✓

**Second-order formula (s2):**
- Expected from paper: `r = sqrt(t³ * c2 / ε)`
- Code computes: `s2 = t * sqrt(c2 / eps_trotter)`
- **Analysis:** 
  - s2 = t * sqrt(c2 / eps_trotter)
  - s2 = t * sqrt(c2 / (E_err * t / 2π))
  - s2 = t * sqrt(c2 * 2π / (E_err * t))
  - s2 = sqrt(t² * c2 * 2π / E_err)
  - s2 = sqrt(t) * sqrt(t * c2 * 2π / E_err)
- Expected: sqrt(t³ * c2 / ε) = sqrt(t) * sqrt(t² * c2 / ε)
- **Problem:** The power of t doesn't match! ❌

### Issue 2: Relationship Between Operator Error and Energy Error

The paper gives bounds on **operator norm** ||S(t) - e^(-itH)||.

The code uses **energy_error** (in Hartrees) and converts it to `eps_trotter`.

**Question:** What is the relationship between these two error measures?

**Hypothesis:** For phase estimation:
- Operator error δ causes phase error Δφ ~ δ
- Phase φ = E*t, so energy error ΔE ~ δ/t
- If we want ΔE ≤ E_err, then δ ≤ E_err * t
- The factor of 2π might relate to the phase wrap-around range

**However:** This relationship is not documented in the code, and the dimensional analysis for s2 doesn't work out.

### Issue 3: error_scale is Deprecated

From README.md: "error_scale: This option is deprecated."

Currently defaults to 1.0, so it doesn't affect anything. But its presence suggests the formulas might have been adjusted in the past.

## What Paper Says

From Childs et al. Equations 145 and 152:

**First-order with r steps:**
```
Total error ≤ r * (t/r)²/2 * c1 = t² * c1 / (2r)
Setting error = ε: r = t² * c1 / (2ε)
```

**Second-order with r steps:**
```
Total error ≤ r * (t/r)³ * c2 = t³ * c2 / r²
Setting error = ε: r² = t³ * c2 / ε
                  r = t^(3/2) * sqrt(c2 / ε)
```

## Comparison with Code

**First-order:**
```python
s1 = t * error_scale * c1 / eps_trotter
   = t * c1 / (E_err * t / 2π)
   = c1 * 2π / E_err
```
This is proportional to `c1 / E_err` but independent of t! ❌

Expected: `r ∝ t² * c1 / E_err` (assuming ε ∝ E_err)

**Second-order:**
```python
s2 = t * sqrt(error_scale * c2 / eps_trotter)
   = t * sqrt(c2 * 2π / (E_err * t))
   = sqrt(t) * sqrt(c2 * 2π / E_err)
```
This is proportional to `sqrt(t) * sqrt(c2 / E_err)` ✓ (assuming ε ∝ E_err)

Expected: `r ∝ t^(3/2) * sqrt(c2 / E_err)`

The power of t is off by a factor of t!

## Possible Explanations

1. **Different error model:** The code might be using a different relationship between operator error and energy error than I assumed.

2. **Effective time:** Perhaps `timestep` represents effective time per Trotter step, not total evolution time?

3. **Normalization convention:** There might be a normalization factor I'm missing.

4. **Error in the code:** The formulas might be incorrect (though the code has been used, so this seems less likely).

5. **Missing factor:** The relationship between eps_trotter and the operator norm error might include an extra factor of t.

## Recommendations

### Immediate Actions

1. ✅ Keep the corrected jkg_utils.py (C1 is correct, not C1/2)
2. ⚠️ **DO NOT modify qre_unitary.py yet** - need more information

### To Verify the Usage

You should:

1. **Check with code author** - They can explain the derivation of s1 and s2
2. **Look for derivation notes** - "Brendan's notes" mentioned in the old TODO
3. **Check DARPA Report #1** mentioned in the old TODO
4. **Test with known cases** - If you have examples where correct # of steps is known
5. **Check related papers** - Look for papers that cite Childs et al. and discuss practical implementation

### Questions to Answer

1. What is the exact relationship between operator norm error ||S(t) - e^(-itH)|| and the energy_error parameter?
2. Does the factor of 2π in eps_trotter have a specific physical meaning?
3. Why does s1 not depend on t (after simplification)?
4. Is there missing documentation about the formulas used?

## Test Case Suggestion

Create a simple test where:
- The exact answer is known analytically
- You can verify the number of Trotter steps needed for a given error
- Compare with what the code computes

This would immediately reveal if the formulas are correct or if there's a bug.
