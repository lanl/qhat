# BCH cancellation hypothesis validation

## Result

This validation panel **supports the within-Hamiltonian cancellation link, but does not yet establish the cross-case claim**. The primary analysis changes only the Pauli schedule within each fixed Hamiltonian and compares every alternative with the signed fermionic reference.

## Coverage

- Manifest cases: 4
- Successful benchmark rows: 180
- Valid raw alternative/reference deltas: 168
- Case/schedule/condition medians: 72
- Primary case/schedule inference units: 24
- Historical performance labels reproduced from current tensors: 4/4
- Historical-label drift: none
- Trotter step counts: [50, 100, 200]
- Evolution times: [1.0]

## Primary within-case test

The variables are `delta log10 cancellation = log10(R_alt/R_signed)` and `delta log10 error = log10(epsilon_alt/epsilon_signed)`. A positive association is the predicted direction.

- Case-centered Pearson r: 0.891
- p-value: 5.47e-09
- Case-block bootstrap 95% interval: [0.835, 0.988]
- Direction-concordant deltas: 70.8%
- Signed stronger-cancellation and lower-error deltas: 45.8%
- Deterministic-control sensitivity: r=0.937, p=9.16e-08, n=16
- Random-parent-block sensitivity (seeds and Trotter conditions collapsed, centered within case): r=0.988, p=5e-16, n=20
- Step-specific Pearson r values: r50_t1=0.891, r100_t1=0.891, r200_t1=0.891

## Held-out cross-case check

The panel cases were not among the eight BCH cases used to formulate the hypothesis. Favorable/control labels were used only to construct the matched panel; the reported advantage below is freshly recomputed from this run's JW-signed, JW-magnitude, and signed-reference rows.

- Best-JW/signed cancellation ratio vs best-JW/signed error ratio: r=0.943, p=0.057, rho=1, n=4
- Case bootstrap 95% interval: [0.785, 1]
- Leave-one-case-out Pearson r range: [0.851, 0.996]
- Relative-cancellation direction predicts the winning ordering in 3/4 cases (one-sided binomial p=0.312)
- Matched-pair delta correlation: r=NA, p=NA, n=2 pairs
- Matched-pair bootstrap 95% interval: [NA, NA]
- Matched-pair direction concordance: 2/2 (one-sided binomial p=0.25)
- Absolute signed cancellation strength (supplemental): r=0.568, p=0.432, n=4

## Schedule medians relative to signed fermionic

| Schedule | BCH cancellation ratio | Error ratio | Concordance |
|---|---:|---:|---:|
| Fermionic magnitude | 0.8903 | 0.9367 | 75.0% |
| JW magnitude | 0.8973 | 1.036 | 75.0% |
| JW signed | 1.179 | 1.309 | 100.0% |
| Random parent blocks | 5.947 | 99.12 | 100.0% |
| Round robin | 2.217 | 5.237 | 75.0% |
| Within-parent shuffle | 1 | 1 | 0.0% |

## Interpretation limits

- This is a BCH-held-out test, not a fully prospective performance test: favorable/control labels came from the existing Trotter sweep.
- Random seeds are summarized by their median before the primary correlation; repeated step counts are also collapsed before primary inference, preventing either from inflating sample size.
- `R_BCH` is evaluated on the HF initial state. A later campaign should repeat the analysis with propagated local-error vectors over time.
- A pilot panel can establish workflow correctness and effect direction, but the full 20-case panel is needed for a report-level mechanism claim.

## Files

- `within_case_deltas.csv`: every raw alternative/reference comparison.
- `primary_aggregated_deltas.csv`: one median per case/schedule/condition.
- `schedule_summary.csv`: schedule-level ratios and concordance.
- `case_summary.csv`: case-level outcome and cancellation results.
- `matched_pair_summary.csv`: within-pair cancellation and error deltas.
- `validation_statistics.csv`: primary and held-out correlations.
- `cancellation_validation.{png,pdf}`: diagnostic figure.
