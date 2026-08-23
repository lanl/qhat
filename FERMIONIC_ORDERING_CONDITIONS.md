# When Fermionic-Aware Ordering Is Likely to Help

## Scope and comparison

This report evaluates `fermionic_signed_coefficient_lexicographic`, the
fermionic-aware ordering in which complete Hermitian fermionic terms are
ordered first and then induce the final JW Pauli-string order.

The fair baseline is the strongest deterministic JW result available for the
same Hamiltonian and numerical settings:

```text
E_best-JW = min(E_JW-signed, E_JW-magnitude-descending)
fermionic advantage = E_best-JW / E_fermionic
```

An advantage larger than 1 favors fermionic-aware ordering. Comparisons use
`one_minus_overlap`, `T = 1`, `r = 100`, and coefficient tolerance `1e-12`.
Cases for which either compared error is at or below `1e-12` are not classified
as wins or losses.

## Main result

- 230 exact JW case/settings pairs were found.
- Fermionic-aware results exist for 215 pairs.
- 201 comparisons remain above the numerical floor.
- Fermionic-aware ordering wins 59, ties 14, and loses 128 comparisons against
  the strongest JW baseline.
- Its global median advantage is 0.952. It is therefore not a universally
  better ordering; its benefit is strongly conditional.

## Conditions associated with good performance

### 1. Large HGBS-5 active spaces

This is the clearest cross-molecule condition in the current data.

- HGBS-5, 18--20 qubits: 17/28 wins, median advantage 1.14, geometric-mean
  advantage 3.59.
- All HGBS-5 cases: 39/109 wins, compared with 20/92 for STO-6G.
- STO-6G, 4--8 qubits: 0/30 wins, with median advantage 0.487.

The effect is not caused by qubit count alone. Large active space and the
orbital/basis structure must occur together.

### 2. Be2 and Li2 are the most consistently favorable molecular families

- Be2: 9/12 wins and median advantage 1.126. All six HGBS-5 Be2 cases win.
- Li2: 9/14 wins and median advantage 1.092. All seven HGBS-5 Li2 cases win.
- F2: 8/16 wins with three ties. Its median is approximately 1, but the
  occupied-heavy HGBS-5 cases at 16, 18, and 20 qubits improve by roughly
  1.8e3, 3.1e3, and 4.4e3, respectively.

These families are the best candidates for demonstrating repeatable rather
than isolated fermionic-ordering gains.

### 3. Some polyatomics become favorable only in larger HGBS-5 spaces

- NH3 HGBS-5: 5/5 valid cases win; median advantage 4.54. The combined NH3
  median falls below 1 because every valid STO-6G case loses.
- H2O HGBS-5 at 16, 18, and 20 qubits improves by approximately 109, 111, and
  141. Small H2O spaces are often much worse, so the overall H2O median is
  only 0.423.
- CH4 HGBS-5 has isolated 18- and 20-qubit gains of approximately 51 and 23,
  while most smaller cases lose.

For polyatomics, the useful prediction is therefore not "this molecule is
good," but "this molecule may become good after the active space contains
enough of the fermionic interaction structure."

## Conditions associated with poor performance

- Small active spaces: only 2/56 cases at 4--8 qubits are fermionic wins.
- Small STO-6G spaces: no wins among the 30 valid 4--8-qubit comparisons.
- BH, LiH, and most BeH cases favor the strongest JW ordering.
- H2O, NH3, and CH4 in small or balanced STO-6G active spaces can be worse by
  one to two orders of magnitude.
- N2, O2, C2, and BeH2 mostly cluster near or below advantage 1; isolated wins
  do not establish a reliable family-wide benefit.

## Mechanistic hypothesis to test

The data support the hypothesis that fermionic-aware ordering helps when a
large active space contains many Pauli descendants whose useful cancellation
or commutator structure is preserved by keeping descendants of related
fermionic parents close together. Small spaces have too little such structure,
and direct Pauli-coefficient ordering is usually sufficient or better.

This is a hypothesis, not yet a proof. A convincing mechanism claim should
also require:

1. the same case to remain favorable over several Trotter step counts;
2. an aligned reduction in the leading BCH-error proxy;
3. robustness to signed versus magnitude fermionic parent ordering; and
4. confirmation that the gain is not a numerical-floor or single-time effect.

### Weighted interaction-graph follow-up

The proposed "distributed versus hub-dominated" explanation was tested on a
coefficient-weighted occupied-to-virtual spin-orbital graph. Each complete
Hermitian fermionic parent contributes its absolute canonical coefficient to
every occupied--virtual support pair it touches. Strength-weighted Freeman
centralization is normalized separately within the occupied and virtual
partitions, and the graph statistic is the larger of the two partition values.
It was then compared with

```text
E_JW-magnitude / E_fermionic
```

for the HGBS-5 sweep.

The result does not confirm the hypothesis. The 85 available, numerically valid
cases have no resolved aggregate correlation between graph centralization and
log advantage (Pearson `r = -0.065`, `p = 0.554`). The matched B-H/Be-H/Li-H
estimate points in the opposite direction (`r = 0.217`, `p = 0.286`), as does
the molecule-and-active-space-size fixed-effect estimate (`r = 0.165`,
`p = 0.182`). The virtual-side centralization—the closest direct test of the
original virtual-hub claim—is also null (`r = -0.120`, `p = 0.275`). In this
sweep, occupied-to-virtual hub dominance is therefore not a systematic
predictor of fermionic-ordering benefit.

The definition, coverage limitations, case table, correlation table, and
figure are in `analysis/fermionic_graph_centralization/`.

### Signed-order cancellation-structure follow-up

The revised mechanism was tested on a parent interaction network whose nodes
are complete Hermitian fermionic parents and whose edge weights are

```text
w_ij = |c_i c_j| * |support_i intersect support_j|.
```

For the signed ascending-coefficient order, the analysis measured the weighted
opposite-sign interaction share, normalized mean edge span, interaction weight
within 5% and 10% of the parent sequence, local opposite-sign enrichment, and
the corresponding improvements over fermionic magnitude ordering. It compared
these quantities against `E_JW-magnitude / E_fermionic`, absolute fermionic
infidelity, and the independently measured HF-state BCH2 cancellation ratio.

The strict structural hypothesis is not verified. There is strong direct
evidence for the cancellation premise itself: across the eight independent BCH
cases, cancellation strength `-log10(R_BCH)` correlates with greater
JW-magnitude/fermionic advantage (`r = 0.852`, `p = 0.0072`) and lower absolute
fermionic infidelity (`r = -0.870`, `p = 0.0050`). However, none of the nine
shared-support structural proxies predicts measured cancellation after FDR
correction. The strongest is local opposite-sign interaction weight
(`r = 0.699`, `p = 0.081`, `q = 0.433`, `n = 7`).

Across the 85 usable full-sweep cases, the global opposite-sign interaction
share correlates with relative advantage after molecule and active-size
adjustment (`r = 0.344`, `q = 0.019`), but that feature is order invariant.
The genuinely order-sensitive locality results are inconsistent with relative
advantage: larger edge span and smaller local weight often accompany more
advantage. They do track lower absolute fermionic infidelity, while no metric
survives FDR correction in the active-size-matched B-H/Be-H/Li-H comparison.
Thus the useful cancellation appears to depend on Pauli-commutator orientations
or BCH-vector alignment that shared parent support and coefficient signs alone
do not encode.

The full case table, matched-size table, correlation/FDR tables, favorable-case
group comparison, report, and figure are in
`analysis/signed_order_cancellation_structure/`.

### Within-Hamiltonian cancellation-ablation validation

The stronger test keeps each final JW Hamiltonian fixed and changes only its
Pauli schedule. A BCH-held-out 12-qubit pilot paired two favorable cases
(Li2 6+6 and BeH 5+7) with two negative controls (H2O 4+8 and LiH 4+8).
Each case used the signed fermionic reference, fermionic magnitude, JW signed,
JW magnitude, parent-descendant round robin, five randomized parent-block
orders, and five within-parent shuffles.

The pilot supports the order--cancellation--error link. Results were repeated
at 50, 100, and 200 Trotter steps. After collapsing step-count repeats to one
median per case and alternative schedule, the case-centered association between
`log10(R_alt/R_signed)` and `log10(E_alt/E_signed)` is Pearson `r = 0.891`
(`p = 5.47e-9`; case-block bootstrap 95% interval `[0.835, 0.988]`). The
step-specific estimates are all `r = 0.891` to the reported precision. The four
deterministic controls give `r = 0.937`, while the randomized parent-block
orders, also collapsed across step counts and centered within case, give
`r = 0.988`. Randomizing intact parent blocks increases the median BCH
cancellation ratio by `5.95x` and the median Trotter error by `99.1x`; shuffling
descendants within each parent changes neither quantity materially. This
localizes the important structure to the relative ordering of fermionic-parent
blocks, not arbitrary descendant order inside a block.

The decisive cross-case test extends the deterministic signed-reference,
JW-signed, and JW-magnitude comparison to 20 BCH-held-out HGBS-5 cases arranged
as ten active-size-matched favorable/control pairs from 12 through 20 qubits.
The best-JW/signed BCH cancellation ratio strongly tracks the freshly
recomputed best-JW/signed error ratio: Pearson `r = 0.807` (`p = 1.73e-5`) and
Spearman `rho = 0.806` (`p = 1.78e-5`). A 10,000-sample case bootstrap gives a
95% interval `[0.349, 0.967]`, and every leave-one-case-out estimate remains
strongly positive (`r = 0.694` to `0.907`). The matched-pair difference test is
also positive (`r = 0.713`, `p = 0.0205`, ten pairs), with the predicted
direction in eight of ten pairs. The sign of the relative cancellation ratio
alone identifies the winning ordering in 17 of 20 cases (one-sided binomial
`p = 0.00129`). The ten-pair bootstrap interval remains wide and crosses zero
(`[-0.179, 0.984]`), so the full-case result is stronger than the pair-only
uncertainty estimate.

Absolute signed cancellation strength remains a poor predictor (`r = 0.242`,
`p = 0.304`). The predictive quantity is therefore the cancellation produced
by one ordering relative to its competing ordering on the same Hamiltonian,
not an absolute molecular cancellation score. Current tensors reproduce 19 of
20 historical favorable/control labels. The exception is F2 14+2, whose old
large favorable label is not reproduced; all reported correlations use fresh
errors and cancellation ratios, so that drift is not silently inherited.

Together, the pilot intervention and full matched panel support the mechanism
claim: fermionic-aware ordering succeeds when its parent-block sequence
produces stronger destructive BCH cancellation than the best direct JW
ordering. The remaining limitation is external validation on prospectively
selected Hamiltonians and, later, cancellation measured along propagated
states rather than only on the HF initial state.

The runner, analysis, raw results, statistical tables, report, and figure are
in `analysis/cancellation_hypothesis_validation/`.

## Practical case-selection rule

Prioritize new overnight studies in this order:

1. HGBS-5 active spaces with 18--20 qubits;
2. Be2 and Li2 HGBS-5 sweeps;
3. occupied-heavy F2 spaces;
4. NH3, H2O, and CH4 HGBS-5 spaces at 16--20 qubits;
5. matched step-scaling and BCH runs for the strongest cases above.

Deprioritize 4--8-qubit STO-6G cases unless they are needed as negative
controls.

## Reproducible outputs

Run from the repository root:

```bash
python analysis/analyze_jw_magnitude_baseline.py
python analysis/analyze_fermionic_aware_performance.py
python analysis/analyze_fermionic_graph_centralization.py
python analysis/analyze_signed_order_cancellation_structure.py
python analysis/run_cancellation_hypothesis_validation.py --panel pilot --execute
python analysis/analyze_cancellation_hypothesis_validation.py
python analysis/run_cancellation_hypothesis_validation.py --panel full \
  --samples 0 --schedules fermionic_signed_reference jw_signed_baseline \
  jw_magnitude_baseline --execute
python analysis/analyze_cancellation_hypothesis_validation.py \
  --input analysis/cancellation_hypothesis_validation/full_ablation_results.csv \
  --manifest analysis/cancellation_hypothesis_validation/full_panel_manifest.csv \
  --outdir analysis/cancellation_hypothesis_validation/full_analysis
```

The case table, molecule summaries, condition summaries, and figures are in
`analysis/fermionic_aware_performance/`.
