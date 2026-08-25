# Two-body weighted-fraction validation

## Verdict

The broad hypothesis is **not supported by the total two-body coefficient-mass fraction alone**. Its raw association with the fermionic error disappears after molecule and active-space-size adjustment. However, the internal composition of the two-body mass does contain reproducible cross-case signal. That narrower signal does not survive as an independent explanation after controlling the measured BCH cancellation ratio in the 20-case matched panel.

## Definitions

The primary fraction is

`f_2b = sum_{two-body parents}|c_p| / sum_{one- and two-body parents}|c_p|`.

Every `c_p` is the canonical real coefficient of one complete Hermitian fermionic parent, matching the signed-ascending production ordering. The constant term is excluded. Two-body composition is also split by HF excitation rank (0/1/2) and distinct spin-orbital support size (2/3/4), with each split normalized by total two-body |coefficient| mass.

Outcomes are `log10(fermionic one-minus-overlap)`, the matching fixed JW magnitude-descending error, and `log10(JW-magnitude error / fermionic error)`.

## Coverage

- Basis: hgbs-5
- Valid performance cases considered: 109
- Fermionic parent structures recovered: 85
- Molecules represented: 11
- Active-space sizes represented: 9
- Unavailable/error cases: 24

## Primary total two-body fraction test

### Absolute fermionic error

- Raw: r=0.418, p=6.98e-05, q=7.97e-05, n=85
- Active-size adjusted: r=0.353, p=0.00162, q=0.00325, n=85
- Molecule + active-size adjusted: r=-0.006, p=0.96, q=0.986, n=85, molecule-bootstrap 95% CI [-0.183, 0.283], LOMO r range [-0.051, 0.086]

The positive raw coefficient means that more total two-body weight is associated with *larger*, not smaller, fermionic error. The complete loss of association after the full adjustment shows that this is a between-molecule/active-space composition effect rather than an ordering-specific predictor.

### Relative fermionic advantage over JW magnitude descending

- Molecule + active-size adjusted: r=-0.179, p=0.146, q=0.293, n=85, molecule-bootstrap 95% CI [-0.462, 0.070], LOMO r range [-0.252, -0.099]

Thus total two-body coefficient mass does not verify a general fermionic-ordering advantage.

## Exploratory two-body composition results

All rows below use molecule + active-space-size adjustment and the reported q-values correct across the eight tested weighted/count features for the same outcome and specification.

### HF excitation-rank composition vs absolute fermionic error

- Rank 0: r=-0.618, p=2.48e-08, q=1.98e-07, n=85, molecule-bootstrap 95% CI [-0.790, -0.374], LOMO r range [-0.670, -0.546]
- Rank 1: r=0.579, p=2.92e-07, q=1.17e-06, n=85, molecule-bootstrap 95% CI [0.323, 0.776], LOMO r range [0.504, 0.639]
- Rank 2: r=0.523, p=5.53e-06, q=1.47e-05, n=85, molecule-bootstrap 95% CI [0.299, 0.762], LOMO r range [0.462, 0.604]

More rank-0 two-body mass accompanies lower fermionic error, whereas more single- and double-excitation mass accompanies higher error. These are compositional associations: the three fractions sum to one.

### Four-orbital two-body mass vs relative advantage

- Relative advantage: r=0.527, p=4.51e-06, q=3.61e-05, n=85, molecule-bootstrap 95% CI [0.183, 0.690], LOMO r range [0.391, 0.623]
- Absolute fermionic error: r=0.002, p=0.986, q=0.986, n=85
- Absolute JW-magnitude error: r=0.435, p=0.000231, q=0.00185, n=85

The four-orbital fraction tracks a larger relative fermionic advantage, but it does not lower the absolute fermionic error. The relative effect arises because the JW magnitude-descending error increases with this fraction while the fermionic error is nearly unchanged.

## Is two-body composition independent of BCH cancellation?

The available 20-case matched cancellation panel was joined without loss. Partial correlations control the independently measured `log10(JW-magnitude/signed BCH cancellation ratio)`.

- Total two-body mass fraction, BCH-adjusted: partial r=0.034, p=0.889, n=20
- Four-orbital two-body fraction, BCH-adjusted: partial r=-0.040, p=0.869, n=20
- Four-orbital fraction, matched-pair + BCH adjusted: partial r=-0.320, p=0.368, n=20

The current data therefore do **not** verify two-body composition as an independent second mechanism beyond BCH cancellation. The full-sweep composition signal is useful as a structural marker, but it may be mediated by the same order-induced cancellation mechanism or by another shared Hamiltonian property.

## Interpretation limits

- This establishes cross-case association, not a causal mechanism independent of BCH cancellation.
- The total two-body fraction was the primary test. The excitation- and support-composition findings are exploratory and need a prospectively selected panel.
- The fixed effects treat molecule and active-space size as categorical controls. Molecule-cluster bootstrap intervals and leave-one-molecule-out ranges quantify robustness to repeated cases.
- Source tensors/pickles were unavailable for some otherwise valid performance rows; no structural values were imputed.
- Only HGBS-5 could be reconstructed from the current repository sources; the STO-6G panel remains an external replication target.

## Files

- `two_body_case_metrics.csv`: case-level parent composition and performance outcomes.
- `two_body_correlations.csv`: raw and adjusted correlation tests, FDR q-values, cluster-bootstrap intervals, and LOMO ranges.
- `two_body_bch_independence.csv`: 20-case partial correlations controlling independently measured BCH cancellation.
- `two_body_weighted_fraction_validation.png`: raw/adjusted diagnostic figure.
