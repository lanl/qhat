# BCH cancellation hypothesis validation

## Result

This validation panel **supports both the within-Hamiltonian and cross-case claims**. The primary analysis changes only the Pauli schedule within each fixed Hamiltonian and compares every alternative with the signed fermionic reference.

## Coverage

- Manifest cases: 20
- Successful benchmark rows: 60
- Valid raw alternative/reference deltas: 40
- Case/schedule/condition medians: 40
- Primary case/schedule inference units: 40
- Historical performance labels reproduced from current tensors: 19/20
- Historical-label drift: F-F_1.28_hgbs-5_as-014-002
- Trotter step counts: [100]
- Evolution times: [1.0]

## Primary within-case test

The variables are `delta log10 cancellation = log10(R_alt/R_signed)` and `delta log10 error = log10(epsilon_alt/epsilon_signed)`. A positive association is the predicted direction.

- Case-centered Pearson r: 0.835
- p-value: 2.2e-11
- Case-block bootstrap 95% interval: [0.0562, 0.962]
- Direction-concordant deltas: 90.0%
- Signed stronger-cancellation and lower-error deltas: 50.0%
- Deterministic-control sensitivity: r=0.835, p=2.2e-11, n=40
- Random-parent-block sensitivity: not rerun in this deterministic full panel; see the pilot analysis.
- Step-specific Pearson r values: r100_t1=0.835

## Held-out cross-case check

The panel cases were not among the eight BCH cases used to formulate the hypothesis. Favorable/control labels were used only to construct the matched panel; the reported advantage below is freshly recomputed from this run's JW-signed, JW-magnitude, and signed-reference rows.

- Best-JW/signed cancellation ratio vs best-JW/signed error ratio: r=0.807, p=1.73e-05, rho=0.806, n=20
- Case bootstrap 95% interval: [0.349, 0.967]
- Leave-one-case-out Pearson r range: [0.694, 0.907]
- Relative-cancellation direction predicts the winning ordering in 17/20 cases (one-sided binomial p=0.00129)
- Matched-pair delta correlation: r=0.713, p=0.0205, n=10 pairs
- Matched-pair bootstrap 95% interval: [-0.179, 0.984]
- Matched-pair direction concordance: 8/10 (one-sided binomial p=0.0547)
- Absolute signed cancellation strength (supplemental): r=0.242, p=0.304, n=20

## Schedule medians relative to signed fermionic

| Schedule | BCH cancellation ratio | Error ratio | Concordance |
|---|---:|---:|---:|
| JW magnitude | 1.038 | 1.017 | 85.0% |
| JW signed | 1.179 | 0.9837 | 95.0% |

## Interpretation limits

- This is a BCH-held-out test, not a fully prospective performance test: favorable/control labels came from the existing Trotter sweep.
- Random seeds are summarized by their median before the primary correlation; repeated step counts are also collapsed before primary inference, preventing either from inflating sample size.
- `R_BCH` is evaluated on the HF initial state. A later campaign should repeat the analysis with propagated local-error vectors over time.
- This full matched panel establishes the current cross-case result; external confirmation still requires prospectively selected tensors.

## Files

- `within_case_deltas.csv`: every raw alternative/reference comparison.
- `primary_aggregated_deltas.csv`: one median per case/schedule/condition.
- `schedule_summary.csv`: schedule-level ratios and concordance.
- `case_summary.csv`: case-level outcome and cancellation results.
- `matched_pair_summary.csv`: within-pair cancellation and error deltas.
- `validation_statistics.csv`: primary and held-out correlations.
- `cancellation_validation.{png,pdf}`: diagnostic figure.
