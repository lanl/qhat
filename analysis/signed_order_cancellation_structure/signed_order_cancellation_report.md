# Signed-order cancellation-structure test

## Result

The revised cancellation-structure hypothesis is **not verified as stated**. The independent BCH result supports its core cancellation premise, but the proposed order-sensitive shared-support metrics do not consistently identify that cancellation structure across the sweep.

## Definition

Nodes are complete Hermitian fermionic parents. For parents `i,j`, the support-overlap interaction weight is `w_ij = |c_i c_j| * |S_i intersect S_j|`. This coefficient-weighted shared-support graph is a pre-BCH proxy, not an exact commutator graph. The signed order is ascending canonical signed coefficient; the control order is descending coefficient magnitude.

The primary metrics include the opposite-sign share of interaction weight, coefficient-weighted mean edge span, weight within 5% and 10% of the parent sequence, local opposite-sign weight/enrichment, and signed-minus-magnitude order improvements. Span is normalized by `N-1`.

## Coverage

- Input HGBS-5 cases: 120
- Parent structures built: 93
- Valid performance comparisons: 85
- Independent BCH cancellation cases: 8
- Unavailable source Hamiltonians: 27

## Performance correlations

The response is `log10(epsilon_JW-magnitude / epsilon_fermionic)`, so positive values favor the fermionic-aware ordering. Reported `r` values below retain each metric's natural sign; the ranking uses the predeclared favorable direction (lower span, higher other metrics).

### Full sweep, raw

- Opposite-sign share of interaction weight: r=0.350, p=0.00102, q=0.00365, n=85
- Signed-order interaction weight within 5%: r=0.182, p=0.0954, q=0.172, n=85
- Signed-order interaction weight within 10%: r=0.141, p=0.198, q=0.297, n=85
- Opposite-sign enrichment in 10% window: r=-0.078, p=0.477, q=0.477, n=85
- Opposite-sign local interaction weight (10%): r=-0.084, p=0.447, q=0.477, n=85

### Full sweep, molecule + active-size adjusted

- Opposite-sign share of interaction weight: r=0.344, p=0.00431, q=0.0194, n=85
- Opposite-sign enrichment in 10% window: r=0.112, p=0.366, q=0.389, n=85
- Opposite-sign local interaction weight (10%): r=0.107, p=0.389, q=0.389, n=85
- Signed-vs-magnitude locality improvement (10%): r=-0.147, p=0.236, q=0.304, n=85
- Signed-vs-magnitude opposite-sign locality improvement: r=-0.187, p=0.129, q=0.194, n=85

### Matched B-H / Be-H / Li-H, active-size adjusted

- Opposite-sign local interaction weight (10%): r=0.401, p=0.0992, q=0.463, n=26
- Opposite-sign enrichment in 10% window: r=0.397, p=0.103, q=0.463, n=26
- Signed-order interaction weight within 10%: r=0.101, p=0.691, q=0.793, n=26
- Signed-order interaction weight within 5%: r=-0.036, p=0.886, q=0.886, n=26
- Signed-vs-magnitude span improvement: r=-0.096, p=0.705, q=0.793, n=26

### Absolute fermionic infidelity, molecule + active-size adjusted

- Signed-order interaction weight within 10%: r=-0.387, p=0.00121, q=0.00368, n=85
- Signed-order interaction weight within 5%: r=-0.384, p=0.00132, q=0.00368, n=85
- Opposite-sign enrichment in 10% window: r=-0.377, p=0.00164, q=0.00368, n=85
- Signed-order weighted mean span: r=0.346, p=0.0041, q=0.00738, n=85
- Signed-vs-magnitude opposite-sign locality improvement: r=-0.332, p=0.0061, q=0.00915, n=85

## Independent cancellation check

For the eight ablation cases, `R_BCH = ||sum_k v_k|| / sum_k ||v_k||` is the independently measured HF-state BCH2 cancellation ratio. Smaller `R_BCH` means more destructive cancellation; the analysis uses `-log10(R_BCH)` as cancellation strength.

- Measured cancellation strength vs performance advantage: r=0.852, p=0.00723, q=0.00723, n=8
- Measured cancellation strength vs absolute fermionic infidelity: r=-0.870, p=0.005, q=0.005, n=8

None of the nine structural proxies significantly predicts the measured cancellation outcomes after FDR correction. The strongest suggestive relation is opposite-sign local interaction weight (10%) vs absolute cancellation strength (r=0.699, p=0.0806, q=0.433, n=7), but it is not conclusive.

Several full-sweep locality results point opposite to the proposed direction after molecule and active-size adjustment: larger signed edge span and smaller local interaction weight accompany more relative advantage. Those same locality metrics do track lower absolute fermionic infidelity, but not the relative advantage or measured cancellation consistently. Therefore the successful BCH cancellation cannot be reduced to these parent-support span/locality statistics.

The detailed CSV also tests every structural feature against both absolute signed-order cancellation strength and the signed-order cancellation advantage over fermionic magnitude ordering.

## Interpretation guardrails

- Shared orbital support is necessary for many interactions but is not an exact commutator magnitude or sign.
- The independent cancellation subset contains eight deliberately selected cases, so its correlations are mechanism checks rather than population estimates.
- Nine related structural metrics are tested; FDR-adjusted q-values are reported within each scope/model.
- Cases missing both active tensors and reconstructable molecular pickles remain unavailable and are not imputed.

## Files

- `signed_order_case_metrics.csv`: case-level structure, outcomes, and cancellation joins.
- `signed_order_performance_correlations.csv`: matched and full-sweep performance tests.
- `signed_order_cancellation_correlations.csv`: independent BCH tests.
- `matched_heteronuclear_by_size.csv`: B-H/Be-H/Li-H cases grouped by active-space size.
- `favorable_structure_group_comparison.csv`: favorable-vs-other feature medians and one-sided rank tests.
- `signed_order_cancellation_diagnostics.{png,pdf}`: diagnostic plot.
