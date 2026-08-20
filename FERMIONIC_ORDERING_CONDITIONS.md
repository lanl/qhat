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
```

The case table, molecule summaries, condition summaries, and figures are in
`analysis/fermionic_aware_performance/`.
