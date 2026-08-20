# Coefficient-weighted fermionic graph centralization

## Definition

Vertices are active spin orbitals, partitioned at the Hartree--Fock occupied/virtual boundary. For each complete Hermitian fermionic parent, the absolute coefficient of its lexicographically canonical component is added to every occupied--virtual orbital pair jointly touched by that parent. Parents confined to one partition do not create edges. This matches the occupied-to-virtual incidence network behind the qualitative hub-dominance hypothesis.

For one partition `P` with `m` orbitals, strengths `s_i = sum_j w_ij`, maximum strength `s_max`, and total bipartite edge weight `W`, its weighted Freeman centralization is

`C_P^w = sum_{i in P}(s_max - s_i) / ((m - 1) W)`.

The graph-level statistic is `max(C_occupied^w, C_virtual^w)`: a graph is hub dominated if either partition concentrates its incident coefficient weight on one orbital. Partition-wise normalization is coefficient-scale invariant, equals zero for equal strength within each partition, equals one for a partition star, and does not mistake an unequal occupied/virtual node count for hub structure.

The response is `E_JW-magnitude / E_fermionic`; values above one favor fermionic-aware ordering. Correlations use its base-10 logarithm, because the ratios span orders of magnitude.

## Coverage

- Graphs built: 93/120 HGBS-5 cases.
- Numerically valid graph/ratio pairs: 85.
- Direct active tensors: 72.
- Reconstructed exactly in memory from retained full molecular pickles: 21.
- Unavailable/error cases: 27.

## Correlation results

- Matched B-H/Be-H/Li-H, unadjusted: n=26, Pearson r=0.217 (p=0.286), Spearman ρ=0.262 (p=0.196).
- Matched B-H/Be-H/Li-H, after removing active-space-size fixed effects: n=26, Pearson r=0.238 (p=0.342), Spearman ρ=0.378 (p=0.122).
- Full available sweep, unadjusted: n=85, Pearson r=-0.065 (p=0.554), Spearman ρ=0.104 (p=0.344).
- Full available sweep, after removing active-space-size fixed effects: n=85, Pearson r=-0.065 (p=0.577), Spearman ρ=-0.068 (p=0.557).
- Full available sweep, after removing molecule and active-space-size fixed effects: n=85, Pearson r=0.165 (p=0.182), Spearman ρ=0.171 (p=0.168).
- Direct active tensors only, unadjusted: n=70, Pearson r=-0.111 (p=0.362), Spearman ρ=0.056 (p=0.646).
- Direct active tensors only, after removing active-space-size fixed effects: n=70, Pearson r=-0.091 (p=0.48), Spearman ρ=-0.036 (p=0.781).

## Interpretation

The hypothesized negative relationship is not supported. The full-sweep estimate points weakly in the proposed direction but is statistically unresolved, and the matched/fixed-effect tests do not recover a negative association.

As a confounding diagnostic, centralization changes with active-space size (Pearson r=0.069, p=0.532), while log ordering advantage changes with active-space size (r=0.446, p=1.86e-05).

The virtual-side centralization—the most direct version of the original virtual-hub claim—also gives no resolved full-sweep association: r=-0.120 (p=0.275). The case table also reports occupied-side, partition-mean, and non-partition-normalized variants for sensitivity analysis.

The fixed-effect rows are the more relevant tests of the mechanistic claim: they compare cases at matched active-space size, and the two-way version also removes stable molecule-family differences. These are observational associations, not a causal identification of the ordering mechanism.

## Missing-source limitation

Coefficient tensors and full molecular pickles were not retained for: C-C, N-N, O-O. Those rows remain in the case and coverage tables with an explicit unavailable reason; they are not silently dropped.

See `fermionic_graph_case_metrics.csv` for the auditable case data, `matched_heteronuclear_cases.csv` for the initial matched comparison, `correlation_summary.csv` for all statistical tests, and `coverage_summary.csv` for source coverage.
