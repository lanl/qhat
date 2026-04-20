"""
Analyze how trotter_error_estimator results are used in qre_unitary.py
and verify dimensional consistency with the paper.
"""

print("="*70)
print("DIMENSIONAL ANALYSIS OF QRE_UNITARY.PY")
print("="*70)

print("""
From qre_unitary.py lines 104-107:

    eps_trotter = config_unitary.energy_error * timestep / (2 * math.pi)
    s1 = timestep * config_unitary.error_scale * c1 / eps_trotter
    s2 = timestep * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)

Let's analyze the dimensions:
""")

print("\nVariables:")
print("  timestep (t):         [time]")
print("  energy_error (E_err): [energy]")
print("  c1:                   [energy]")
print("  c2:                   [energy]")
print("  error_scale:          [dimensionless]")
print()

print("Derived quantities:")
print("  eps_trotter = E_err * t / (2π)")
print("              = [energy] * [time] / [1]")
print("              = [dimensionless] (in units where ℏ=1)")
print()

print("First-order step estimate:")
print("  s1 = t * error_scale * c1 / eps_trotter")
print("     = [time] * [1] * [energy] / [1]")
print("     = [energy * time] = [dimensionless]")
print("     But this is the NUMBER OF STEPS, which should be dimensionless ✓")
print()

print("Second-order step estimate:")
print("  s2 = t * sqrt(error_scale * c2 / eps_trotter)")
print("     = [time] * sqrt([1] * [energy] / [1])")
print("     = [time] * sqrt([energy])")
print("     = [time] * [energy^(1/2)]")
print("     In units where ℏ=1: [time] ~ [energy^(-1)]")
print("     So: [energy^(-1)] * [energy^(1/2)] = [energy^(-1/2)] ≠ [1] ❌")
print("     This is NOT dimensionless!")
print()

print("="*70)
print("CHECKING AGAINST PAPER FORMULAS")
print("="*70)

print("""
From Childs et al. (arXiv:1912.08854v3):

First-order error with r steps of size Δt = t/r:
  ||S₁(t) - e^(-itH)|| ≤ r * (Δt)²/2 * c1 = t² * c1 / (2r)

Second-order error with r steps:
  ||S₂(t) - e^(-itH)|| ≤ r * (Δt)³ * c2 = t³ * c2 / r²

Setting error equal to ε (operator norm error):
  First-order:  r = t² * c1 / (2ε)
  Second-order: r = t * sqrt(c2/ε) * sqrt(t) = sqrt(t³ * c2 / ε)

But the code uses eps_trotter = E_err * t / (2π), not just ε.
""")

print("\nKey question: What is the relationship between operator norm error")
print("(||S(t) - e^(-itH)||) and energy error (E_err)?")
print()
print("Hypothesis: For phase estimation, the phase error Δφ relates to")
print("operator error δ approximately as Δφ ~ δ. Since φ = E*t, we have")
print("ΔE ~ δ/t. If we want ΔE ≤ E_err, then δ ≤ E_err * t.")
print()
print("But eps_trotter = E_err * t / (2π), suggesting a factor of 2π somewhere.")
print("This might relate to the energy range being 2π/t (line 103).")

print("\n" + "="*70)
print("VERIFICATION NEEDED")
print("="*70)
print("""
To fully verify the usage in qre_unitary.py, we need to:

1. Check if there's documentation explaining the relationship between
   operator norm error and energy error for the specific application

2. Verify the dimensional analysis for s2 (currently looks wrong)

3. Look for any papers or references that explain the formulas used

4. Check if error_scale is meant to compensate for dimensional issues

5. Test with known examples where the correct number of steps is known

RECOMMENDATION: This requires input from the code author or domain expert
to clarify the intended formulas and their derivation.
""")
