# Two-body BCH cancellation case study

## Verdict

The two-body contribution hypothesis is **supported as a direct participating component of the BCH mechanism**.

The strongest evidence is not that the 2B-bearing vector always cancels more strongly by itself. It is that signed fermionic ordering makes the 2B-bearing vector combine more destructively with the remaining 1B-only BCH vector. The static amount of two-body BCH mass has no performance association.

The primary order-sensitive quantity is `G_2B = log10(R_2B,JW / R_2B,F)`, where `R_2B` is the HF-state BCH2 cancellation ratio of the exact BCH amplitude containing at least one two-body fermionic coefficient contribution. Positive `G_2B` means that signed fermionic ordering cancels that component more strongly than JW magnitude-descending ordering.

## Exact decomposition

For each final JW coefficient, `c_i = c_i^(1B) + c_i^(2B)`. Every noncommuting-pair coefficient product is therefore decomposed as

`c_i c_j = c_i^(1B)c_j^(1B) + [c_i^(1B)c_j^(2B)+c_i^(2B)c_j^(1B)] + c_i^(2B)c_j^(2B)`.

The three complex HF-state BCH vectors sum back to the full measured BCH vector. This is coefficient-provenance decomposition, not a graph support proxy. Ordering changes only the orientation sign of each fixed Pauli-pair amplitude.

## Coverage and numerical checks

- Cases: 20 (10 favorable, 10 matched controls)
- Matched pairs: 10
- Active-space range: 12-20 qubits
- Maximum final-Pauli coefficient reconstruction error: 8.882e-16
- Maximum floating-point reconciliation applied before exact reconstruction: 1.244e-08
- Maximum BCH pair-component reconstruction error: 2.220e-16
- Maximum discrepancy from independently saved full BCH ratios: 3.500e-16

## Primary result: two-body cancellation vs Trotter advantage

- Raw case-level: r=0.74, p=0.000193, q=0.000677, n=20
- Matched-pair adjusted: r=0.737, p=0.000207, q=0.000723
- Controlling 1B-only cancellation gain: r=0.68, p=0.000962, q=0.00289
- Controlling both matched pair and 1B-only gain: r=0.641, p=0.00232, q=0.00696
- Pair-difference correlation: r=0.751, p=0.0124, n=10 pairs
- Favorable-minus-control mean `G_2B`: 0.0466, bootstrap 95% CI [-0.0292, 0.121], exact paired permutation p=0.271
- Favorable case has larger `G_2B` in 5/10 matched pairs.

Median `G_2B` is 0.057 in favorable cases and 0.0052 in controls.
The significant continuous correlation therefore does not amount to clean favorable/control classification by 2B-internal cancellation alone.

## Does the two-body component track the measured full BCH effect?

- `G_2B` vs full BCH cancellation advantage: r=0.644, p=0.00217, q=0.0076
- Matched-pair adjusted: r=0.632, p=0.00282

## Cross-component interference checks

The destructive cosine is positive when the aggregated 1B-only and 2B-bearing vectors point against each other. The leave-in statistic is positive when adding the 2B-bearing vector lowers the full norm relative to the 1B-only counterfactual.

- Signed-vs-JW gain in destructive cosine vs Trotter advantage: r=0.601, p=0.00509, q=0.00594; controlling 1B-only gain, r=0.525, p=0.0174, q=0.0347; favorable direction in 9/10 pairs, exact paired p=0.00391
- Signed-vs-JW gain in two-body leave-in norm reduction vs advantage: r=0.903, p=4.87e-08, q=3.41e-07; controlling 1B-only gain, r=0.871, p=5.91e-07, q=3.54e-06; favorable direction in 9/10 pairs, exact paired p=0.00391

## Illustrative matched cases

These examples are descriptive selections from the completed panel, not additional inferential tests.

- Largest performance contrast (20q_polyatomic): `H2O_s-1.00_hgbs-5_as-008-012` vs `B-H_1.29_hgbs-5_as-006-014` has a log-error-advantage difference of 2.62. Their `G_2B` values are 0.187 and -0.00126; their 2B leave-in gains are 0.107 and -0.036.
- Internal-cancellation counterexample (12q_pair_a): `Li-Li_2.66_hgbs-5_as-006-006` performs better even though its `G_2B` is -0.000101, below the control value 0.142. Its 2B leave-in gain is nevertheless less unfavorable (-0.00188 vs -0.0109). This is why the cross-vector mechanism is more complete than a 2B-internal-cancellation-only story.

## Specificity controls

- Static 2B-bearing BCH mass fraction vs advantage: r=0.0718, p=0.764, q=0.764
- 1B-only cancellation gain also correlates with advantage: r=0.615, p=0.00391, q=0.00547
- The 2B internal, cross-interference, and leave-in signals remain after controlling the 1B-only gain, so they are not merely copies of that negative-control component. The result nevertheless identifies 2B as a participant, not as the sole source of cancellation.

## Interpretation

The combined evidence identifies the two-body-bearing portion of the leading BCH vector as part of the ordering mechanism. It does not mean that a larger static two-body weighted fraction is beneficial. The relevant quantity is the change in vector cancellation caused by the order. Conversely, a null primary association means the existing full BCH effect cannot be localized consistently to two-body-bearing terms with this decomposition.

This remains an HF-state leading-BCH mechanism study. It is not a causal term-deletion experiment, and higher BCH orders or other initial states can redistribute contributions.

## Files

- `two_body_bch_case_metrics.csv`: exact per-case component metrics.
- `two_body_bch_matched_pairs.csv`: favorable-control differences.
- `two_body_bch_statistics.csv`: raw, matched, paired, FDR, and bootstrap tests.
- `two_body_bch_case_study.png`: diagnostic figure.
