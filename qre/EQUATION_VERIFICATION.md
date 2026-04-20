# Verification of Trotter Error Formulas

This document verifies the equations in `jkg_utils.py` and `jkg_utils_fast.py` against the paper:
**"Theory of Trotter Error with Commutator Scaling"** by Childs et al. (2021)
arXiv:1912.08854v3

## Summary of Findings

⚠️ **POTENTIAL ISSUE FOUND**: There is a TODO comment in the code questioning whether C1 should be divided by 2.

## Paper Background

The Childs et al paper provides bounds on Trotter error for product formulas. For a Hamiltonian decomposed as:
```
H = ∑ᵢ Hᵢ
```

The second-order Trotter formula error is bounded using commutator-based coefficients.

## Key Equations from Paper

### Equation 145 (Second-Order Product Formula)
For the second-order symmetric product formula S₂(t), the error bound involves:
- **First-order commutator terms**: Related to ||[Hᵢ, Hⱼ]||
- **Second-order nested commutator terms**: Related to ||[Hᵢ, [Hⱼ, Hₖ]]||

The paper defines error coefficients that involve sums over all pairs and triples of Hamiltonian terms.

## Code Implementation

### What the code computes

In `jkg_utils.py` lines 119-209, the code estimates three quantities via Monte Carlo sampling:

1. **C1**: Sum over all pairs (i < j) of ||[Hᵢ, Hⱼ]||
   ```python
   C1 = ∑ᵢ<ⱼ ||[Hᵢ, Hⱼ]||
   ```

2. **C21**: Sum over triples (k < i, k < j) of ||[Hᵢ, [Hⱼ, Hₖ]]||
   ```python
   C21 = ∑ₖ<ⱼ,ₖ<ᵢ ||[Hᵢ, [Hⱼ, Hₖ]]||
   ```

3. **C22**: Sum over pairs (k < j) of ||[Hₖ, [Hₖ, Hⱼ]]||
   ```python
   C22 = ∑ₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]||
   ```

### What the code returns

Line 209 of `jkg_utils.py`:
```python
return C1_est, C21_est/12 + C22_est/24
```

So the function returns:
- `c1 = C1`
- `c2 = C21/12 + C22/24`

### How it's used

In `qre_unitary.py` lines 106-107:
```python
s1 = timestep * config_unitary.error_scale * c1 / eps_trotter
s2 = timestep * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)
```

This suggests:
- First-order error: ∝ c1 * t / ε
- Second-order error: ∝ sqrt(c2 / ε) * t

For second-order formulas, the number of Trotter steps is: `r ~ t * sqrt(c2 / ε)`

## Issues to Investigate

### 1. TODO Comment on Line 207
```python
# TODO: Should we return C1_est/2 instead of C1_est?  See, for example, Childs et al Eq 145,
#       DARPA Report #1, Brendan's notes from studying this with SymPy.
```

**Question**: Does the paper use C1 or C1/2 in the error formula?

**Analysis needed**: 
- The sum over pairs (i < j) counts each unordered pair once
- However, commutators are antisymmetric: [Hᵢ, Hⱼ] = -[Hⱼ, Hᵢ]
- Some formulations sum over all i ≠ j and then divide by 2 to avoid double-counting
- Need to verify which convention the paper uses

### 2. Coefficients 1/12 and 1/24

The code returns `C21/12 + C22/24`. 

**Question**: Are these coefficients correct according to Eq 145?

**Analysis needed**:
- These appear to be symmetry factors or combinatorial factors
- Need to verify against the paper's exact formula

### 3. Commutator Norm Formula

In both implementations, the commutator norm is computed as:
```python
norm = 2 * |c1 * c2|  if anticommute else 0
```

This is based on the fact that for anticommuting Pauli operators P and Q:
```
[P, Q] = 2PQ
```

**Verification**: This appears correct for Pauli operators.

## Recommendations

1. **Obtain and check Equation 145** from the paper to verify:
   - Whether C1 should include a factor of 1/2
   - The exact coefficients for C21 and C22 terms
   - The precise definition of the error bound

2. **Check DARPA Report #1** mentioned in the TODO comment

3. **Review Brendan's notes** mentioned in the TODO comment

4. **Add reference comments** to the code indicating which equation each formula corresponds to

## Testing

The current test suite in `test_correctness.py` verifies:
- ✅ Both implementations use the same formulas
- ✅ Commutator norms satisfy Pauli algebra rules
- ✅ Anticommutation checks follow expected patterns

But it does NOT verify:
- ❌ The absolute correctness of the coefficients (1/12, 1/24)
- ❌ Whether C1 should be C1/2
- ❌ That the formulas match the paper exactly

## Detailed Analysis of Potential Issues

### Issue 1: The C1/2 Question

The code computes:
```python
C1 = ∑ᵢ<ⱼ ||[Hᵢ, Hⱼ]||
```

But there are two conventions:

**Convention A (Current Code)**:
- Sum over unordered pairs i < j
- Each pair counted once
- C1 = ∑ᵢ<ⱼ ||[Hᵢ, Hⱼ]||

**Convention B (Possible Alternative)**:
- Sum over all i ≠ j
- Each pair counted twice (once as (i,j), once as (j,i))
- Then divide by 2: C1 = (1/2) * ∑ᵢ≠ⱼ ||[Hᵢ, Hⱼ]||
- Since ||[Hᵢ, Hⱼ]|| = ||[Hⱼ, Hᵢ]||, this gives the same result

**Current code uses Convention A**, which is mathematically equivalent to Convention B.

However, **if the paper uses Convention B notation**, then the code should return `C1_est/2` to match the paper's definition, even though the computed value is correct.

### Issue 2: Sampling Distribution

Looking at lines 119-140 in `jkg_utils.py`:

```python
while time.time() - start_time < time_limit/3:
    for _ in range(batch_size):
        i = np.random.randint(0, N)
        j = np.random.randint(0, N - 1)
        if j >= i:
            j += 1
        op, norm_val = compute_commutator(pauli_terms[i], pauli_terms[j])
        key = tuple(sorted((i, j)))
        comm_dict[key] = (op, norm_val)
        C1_sum += norm_val
        samples_C1 += 1
```

This samples **ordered pairs** (i, j) with i ≠ j uniformly at random. Each unordered pair {i, j} can be sampled in two ways: (i,j) or (j,i).

The estimator is then:
```python
C1_est = C1_sum * (total_C1 / samples_C1)
```
where `total_C1 = N * (N - 1) / 2`

**Analysis**: 
- Total number of ordered pairs: N(N-1)
- Total number of unordered pairs: N(N-1)/2
- Each sample contributes ||[Hᵢ, Hⱼ]|| to C1_sum
- Expected value of one sample: (1/[N(N-1)]) * ∑ᵢ≠ⱼ ||[Hᵢ, Hⱼ]||
- After scaling by N(N-1)/2: **This gives ∑ᵢ<ⱼ ||[Hᵢ, Hⱼ]||**

✅ **The sampling is correct for Convention A**.

But if the paper uses the convention where C1 is defined as:
```
C1 = (1/2) ∑ᵢ≠ⱼ ||[Hᵢ, Hⱼ]|| = ∑ᵢ<ⱼ ||[Hᵢ, Hⱼ]||
```
and then writes error bounds in terms of this C1, the code is computing the right quantity.

However, if the paper defines:
```
C1 = ∑ᵢ≠ⱼ ||[Hᵢ, Hⱼ]||
```
then the code should return `C1_est/2`.

### Issue 3: Physical Interpretation

For **first-order** product formulas, the error is typically:
```
Error ∝ t² * ∑ᵢ<ⱼ ||[Hᵢ, Hⱼ]||
```

For **second-order** product formulas (like S₂), the first-order commutator terms cancel, and the error is:
```
Error ∝ t³ * (nested commutator terms)
```

Looking at how it's used in `qre_unitary.py`:
```python
s1 = timestep * config_unitary.error_scale * c1 / eps_trotter
s2 = timestep * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)
```

For first-order: `steps ∝ t * c1 / ε`
- This means: `error ∝ (t/steps) * t * c1 = (t²/steps) * c1`
- Setting `error = ε` gives: `steps ∝ t² * c1 / ε` ❌ **This doesn't match!**

**Expected**: `steps ∝ t * c1 / ε` to get `error = (t/steps) * t * c1 = (t²/r) * c1 = ε`

For second-order: `steps ∝ t * sqrt(c2 / ε)`
- This means: `error ∝ (t/steps)² * t * c2 = (t³/steps²) * c2`
- Setting `error = ε` gives: `steps = t * sqrt(c2 * t / ε)` ❌ **This doesn't quite match either!**

**Something seems off in the dimensional analysis.**

## Critical Questions for the User

1. **Can you access the PDF** to read Equations 144-145 directly? The file is at:
   ```
   /vast/home/bkkrueger/.claude/projects/.../webfetch-1776705006385-jfcs70.pdf
   ```

2. **Do you have access to "DARPA Report #1"** mentioned in the TODO comment?

3. **Can you find "Brendan's notes"** mentioned in the TODO comment?

4. **Are there any other team members** who worked on this code who might know the answer?

## Next Steps

To complete verification:
1. **Read Equations 144-145** from the PDF using a PDF reader
2. Check the **exact definition of C₁** in the paper
3. Check the **exact coefficients** for nested commutator terms (verify 1/12 and 1/24)
4. Verify the **dimensional analysis** for how c1 and c2 are used in qre_unitary.py
5. Resolve the C1 vs C1/2 question
6. Add equation references to code comments

## Temporary Workaround

If you cannot access the PDF, you could:
1. Compare results with a known working implementation
2. Run small test cases where the answer can be computed exactly
3. Contact the paper authors or check if there's errata
