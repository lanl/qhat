# Verification Results: Code vs. Childs et al. Paper

## Equations Found in Paper

### Equation 145 (First-Order Lie-Trotter Formula)

For Hamiltonian H = Σᵢ₌₁ᴸ Hᵢ and first-order formula S₁(t) = Πᵢ₌₁ᴸ e^{-itHᵢ}:

```
||S₁(t) - e^{-itH}|| ≤ (t²/2) * Σ₁₁₌₁ᴸ (Σ₁₂₌₁₁₊₁ᴸ ||[H₁₂, H₁₁]||)
```

Expanding the double sum:
```
= (t²/2) * Σᵢ<ⱼ ||[Hⱼ, Hᵢ]||
= (t²/2) * Σᵢ<ⱼ ||[Hᵢ, Hⱼ]||  (since ||[Hⱼ, Hᵢ]|| = ||[Hᵢ, Hⱼ]||)
```

### Equation 152 (Second-Order Suzuki Formula)

For second-order formula S₂(t) = Π₁ᴸ e^{-it/2 Hᵢ} Πᵢ₌₁ᴸ e^{-it/2 Hᵢ}:

```
||S₂(t) - e^{-itH}|| ≤ (t³/12) * Σᵢ₌₁ᴸ ||[Σₖ₌ᵢ₊₁ᴸ Hₖ, [Σⱼ₌ᵢ₊₁ᴸ Hⱼ, Hᵢ]]||
                      + (t³/24) * Σᵢ₌₁ᴸ ||[Hᵢ, [Hᵢ, Σⱼ₌ᵢ₊₁ᴸ Hⱼ]]||
```

## What the Code Computes

From `jkg_utils.py` lines 119-209:

```python
C1 = Σᵢ<ⱼ ||[Hᵢ, Hⱼ]||

C21 = Σₖ<ᵢ,ₖ<ⱼ ||[Hᵢ, [Hⱼ, Hₖ]]||

C22 = Σₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]||

Returns: (C1, C21/12 + C22/24)
```

## Analysis

### ✅ First-Order Formula (C1)

**Paper:** Error ≤ (t²/2) * Σᵢ<ⱼ ||[Hᵢ, Hⱼ]||

**Code:** Computes C1 = Σᵢ<ⱼ ||[Hᵢ, Hⱼ]|| and returns it directly

**Conclusion:** ✅ **CORRECT!** The code correctly computes the coefficient for first-order error.

**Resolution of TODO:** The TODO comment on line 207 questioning "Should we return C1_est/2?" is **INCORRECT**. The code should return C1 as it currently does. The factor of 1/2 is already in the error formula (t²/2), not in the definition of C1.

### ⚠️ Second-Order Formula (C2)

**Paper equation 152 uses nested commutators with SUMS inside:**
- First term: Σᵢ ||[Σₖ>ᵢ Hₖ, [Σⱼ>ᵢ Hⱼ, Hᵢ]]||
- Second term: Σᵢ ||[Hᵢ, [Hᵢ, Σⱼ>ᵢ Hⱼ]]||

**Applying triangle inequality to expand these:**

First term:
```
Σᵢ ||[Σₖ>ᵢ Hₖ, [Σⱼ>ᵢ Hⱼ, Hᵢ]]|| 
  ≤ Σᵢ Σₖ>ᵢ Σⱼ>ᵢ ||[Hₖ, [Hⱼ, Hᵢ]]||
```

Second term:
```
Σᵢ ||[Hᵢ, [Hᵢ, Σⱼ>ᵢ Hⱼ]]|| 
  ≤ Σᵢ Σⱼ>ᵢ ||[Hᵢ, [Hᵢ, Hⱼ]]||
```

**Code computes:**
- C21 = Σₖ<ᵢ,ₖ<ⱼ ||[Hᵢ, [Hⱼ, Hₖ]]|| (all triples where k < i AND k < j)
- C22 = Σₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]|| (all pairs where k < j)

### 🔍 Index Ordering Issue

**Problem:** The paper's formula after triangle inequality gives:
```
Σᵢ₌₁ᴸ Σⱼ₌ᵢ₊₁ᴸ Σₖ₌ᵢ₊₁ᴸ ||[Hₖ, [Hⱼ, Hᵢ]]||
```
This sums over: **i < j, i < k** (where j and k are independent, both > i)

**Code computes C21 as:**
```
Σₖ<ᵢ,ₖ<ⱼ ||[Hᵢ, [Hⱼ, Hₖ]]||
```
This sums over: **k < i, k < j** (where i and j are independent, both > k)

**These are related by relabeling!** If we relabel indices in the code's formula:
- Let k → i (innermost operator)
- Let i → k (outer operator in nested commutator)
- Let j → j (middle operator)

Then the code's C21 becomes:
```
Σᵢ<ₖ,ᵢ<ⱼ ||[Hₖ, [Hⱼ, Hᵢ]]|| = Σᵢ Σⱼ>ᵢ Σₖ>ᵢ ||[Hₖ, [Hⱼ, Hᵢ]]||
```

✅ **This matches the paper!**

### ✅ Second Term (C22)

**Paper (second term):**
```
Σᵢ₌₁ᴸ Σⱼ₌ᵢ₊₁ᴸ ||[Hᵢ, [Hᵢ, Hⱼ]]||
```

**Code computes:**
```
C22 = Σₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]||
```

Relabeling k → i:
```
C22 = Σᵢ<ⱼ ||[Hᵢ, [Hᵢ, Hⱼ]]||
```

✅ **This matches the paper!**

### ✅ Coefficients 1/12 and 1/24

**Paper equation 152:**
```
Error ≤ (t³/12) * [first term] + (t³/24) * [second term]
```

**Code returns:**
```
c2 = C21/12 + C22/24
```

✅ **The coefficients are correct!**

## Final Verification

### First-Order (Equation 145)
- ✅ Code correctly computes C1 = Σᵢ<ⱼ ||[Hᵢ, Hⱼ]||
- ✅ Error bound: ||S₁(t) - e^{-itH}|| ≤ (t²/2) * C1
- ❌ TODO comment suggesting C1/2 should be **REMOVED** - it's incorrect

### Second-Order (Equation 152)
- ✅ Code correctly computes C21 (matches paper after index relabeling)
- ✅ Code correctly computes C22 (matches paper after index relabeling)  
- ✅ Coefficients 1/12 and 1/24 are correct
- ✅ Error bound: ||S₂(t) - e^{-itH}|| ≤ t³ * (C21/12 + C22/24)

## How the Error Bounds Are Used

From `qre_unitary.py` lines 106-107:
```python
s1 = timestep * config_unitary.error_scale * c1 / eps_trotter
s2 = timestep * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)
```

For **first-order Trotter** with r steps:
```
Error ≤ r * (t/r)² * (c1/2) = t² * c1 / (2r)
Setting Error = ε: r ≥ t² * c1 / (2ε)
```
So: `r ~ t * c1 / ε` ✅ This matches `s1 ~ timestep * c1 / eps_trotter`

For **second-order Trotter** with r steps:
```
Error ≤ r * (t/r)³ * c2 = t³ * c2 / r²
Setting Error = ε: r ≥ t * sqrt(c2 / ε)
```
So: `r ~ t * sqrt(c2 / ε)` ✅ This matches `s2 ~ timestep * sqrt(c2 / eps_trotter)`

## Recommendations

1. **REMOVE the TODO comment** on line 207 of `jkg_utils.py` - it's incorrect
2. **Add equation references** to the code comments pointing to equations 145 and 152
3. **Add a comment** explaining the index relabeling between code and paper notation
4. ✅ **No changes needed to the actual computations** - they are correct!

## Updated Comment for Code

Suggested comment to add to `jkg_utils.py`:
```python
# This function implements the Monte Carlo estimation of Trotter error coefficients
# from Childs et al., "Theory of Trotter Error" (arXiv:1912.08854v3)
#
# For first-order formulas (Equation 145):
#   Error ≤ (t²/2) * C1, where C1 = Σᵢ<ⱼ ||[Hᵢ, Hⱼ]||
#
# For second-order formulas (Equation 152):
#   Error ≤ (t³/12) * C21 + (t³/24) * C22, where
#   C21 = Σₖ<ᵢ,ₖ<ⱼ ||[Hᵢ, [Hⱼ, Hₖ]]||  (all triples with k < i and k < j)
#   C22 = Σₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]||       (all pairs with k < j)
#
# Note: The index ordering differs from the paper by relabeling, but the sums are equivalent.
```
