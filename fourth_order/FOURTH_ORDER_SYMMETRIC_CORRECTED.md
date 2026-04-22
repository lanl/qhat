# Fourth-Order Trotter: Symmetric 5-Step Suzuki Formula

## Important Note About Negative Time Steps

The symmetric 5-step Suzuki formula is:

```
S4(t) = S2(p·t) · S2(p·t) · S2((1-4p)·t) · S2(p·t) · S2(p·t)
```

where **p = 1/(4 - 4^(1/3)) ≈ 0.4145**

This gives **(1-4p) ≈ -0.658**, which is **NEGATIVE**.

### Why Is There a Negative Time Step?

The Suzuki construction **requires** negative time steps to achieve fourth-order accuracy. The middle step goes "backward in time" and this is essential for canceling the third-order error terms.

**This is not a bug—it's a feature!**

### Alternatives Without Negative Time Steps

If you want to avoid negative time steps, you have two options:

1. **Use a different fourth-order formula** (e.g., Forest-Ruth FSAL methods)
2. **Stay with second-order** (which never needs negative time)

However, most practical fourth-order methods DO use negative substeps.

## Clarification Needed

The phrase "time remains within the time slice" is ambiguous:

**Interpretation A:** Cumulative time stays in [0, t]
- This IS satisfied: even though middle step is negative, cumulative progress stays within bounds

**Interpretation B:** All fractional times are positive  
- This is NOT satisfied by standard Suzuki formula
- Would need a different method (like forward-difference formulas)

## The Formula (Assuming Negative Steps Are OK)

### Error Bound

```
Error ≤ |t|⁵ · C₄
```

Note the absolute value on t, since we're using negative substeps.

### Coefficient Formula

Since (1-4p) < 0, the weight factor has a sign issue. The correct approach is to use absolute values in the error analysis:

```
C₄ = |C41/α₁ + C42/α₂ + C43/α₃ + C44/α₄|
```

where the denominators αᵢ account for the weight:

```
weight = |4p⁵ + (1-4p)⁵| = |0.0489 - 0.1233| = 0.0744
```

This gives effective denominators:
```
α₁ ≈ 120/0.0744 ≈ 1613
α₂ ≈ 240/0.0744 ≈ 3227  
α₃ ≈ 480/0.0744 ≈ 6454
α₄ ≈ 240/0.0744 ≈ 3227
```

### Practical Formula

```
C₄ = |C41/1600 + C42/3200 + C43/6400 + C44/3200|
```

(Rounded to nice numbers)

### Step Count

```
r = ⌈|t|^(5/4) · (C₄/ε)^(1/4)⌉
```

Each fourth-order step requires **5 second-order evaluations**.

## Alternative: Positive-Time-Only Fourth-Order

If you specifically need all time steps to be positive, you'd need a different composition. One option:

### 7-Step Positive-Time Formula (Yoshida)

```
S4(t) = S2(w₁t) · S2(w₂t) · S2(w₃t) · S2(w₄t) · S2(w₃t) · S2(w₂t) · S2(w₁t)
```

where:
```
w₁ = w₂ = w₃ = 1/(2 - 2^(1/3)) ≈ 1.3512
w₄ = 1 - 4w₁ ≈ -4.4048
```

**Problem:** w₄ is still negative! In fact, **all known fourth-order methods require negative time steps**.

### Forward-Difference Methods

Some fourth-order methods use only forward steps but require:
- Extra function evaluations
- Storage of intermediate states
- More complex bookkeeping

These are typically used in ODE solvers, not Hamiltonian simulation.

## Recommendation

**I need clarification on what you want:**

1. **Standard Suzuki with negative middle step?**
   - Most common in quantum simulation
   - Formula above applies
   - 5 S2 evaluations per S4 step

2. **Alternative method without negative steps?**
   - Would need different formula entirely
   - Likely more expensive
   - May need higher-order methods

3. **Just use second-order?**
   - Always positive time steps
   - Well-understood error bounds
   - Current implementation already works

## My Recommendation

For **Hamiltonian simulation**, the negative time step in Suzuki's formula is standard practice and not problematic:

- The Hamiltonian is Hermitian → time evolution is unitary in both directions
- exp(-iH·(-τ)) = exp(iH·τ) is just the adjoint
- Negative time just means running the exponential "in reverse"

**This is completely valid and widely used.**

If this is acceptable, I can provide the complete error formula with proper treatment of the negative step. If you specifically need positive-time-only methods, we should discuss alternative approaches.

## What Do You Want?

Please clarify:
1. Are negative time steps acceptable?
2. Do you want the standard 5-step Suzuki formula (with negative middle step)?
3. Or do you need a different fourth-order method entirely?

Based on your answer, I can provide the exact formula and implementation guidance.
