# Pre-BCH Structural Features for Fermionic-Aware Ordering

## Scope

All predictor features below are computed **before** BCH evaluation and use only fermionic parent coefficients, orbital indices, occupied/virtual excitation structure, and parent/orbital graph structure. No BCH norm, commutator-vector norm, or BCH cancellation quantity is used as a feature.

When an outcome CSV is available, the external target is `one_minus_overlap` at T=1 and r=100, with the fair JW baseline `min(JW signed, JW magnitude-descending)`.

## Four-case outcome check

| Molecule | active occ+virt | parents | best-JW / fermionic | label |
|---|---:|---:|---:|---|
| H2O | 8+10 | 1367 | 111 | win |
| NH3 | 8+10 | 3907 | 42.66 | win |
| BeH2 | 6+12 | 1131 | 0.9112 | loss |
| O2 | 8+10 | 1293 | 0.9989 | loss |

## Best pre-BCH separators

The ranking favors features for which H2O and NH3 are mutually consistent but O2 is well separated; when all four external advantages are available, a small Spearman correlation term is added. With only four cases these are descriptive diagnostics, not a fitted predictive model.

| feature | H2O | NH3 | BeH2 | O2 | direction | score |
|---|---:|---:|---:|---:|---|---:|
| `excitation_rank2_coefficient_mass_fraction` | 0.0014 | 7.411e-04 | 0.0030 | 0.0060 | O2 higher | 0.758 |
| `virtual_orbital_coefficient_mass_gini` | 0.1668 | 0.1250 | 0.1157 | 0.3895 | O2 higher | 0.665 |
| `virtual_orbital_coefficient_mass_top1_fraction` | 0.1449 | 0.1277 | 0.1003 | 0.2883 | O2 higher | 0.651 |
| `ov_coupling_top_singular_fraction` | 0.5389 | 0.6333 | 0.7766 | 0.8529 | O2 higher | 0.646 |
| `virtual_orbital_coefficient_mass_entropy` | 0.9809 | 0.9880 | 0.9888 | 0.8513 | H2O/NH3 higher | 0.624 |
| `parent_orbital_incidence_top_singular_fraction` | 0.0952 | 0.0963 | 0.0979 | 0.1044 | O2 higher | 0.570 |
| `mixed_occ_virtual_coefficient_mass_fraction` | 0.1782 | 0.1623 | 0.2249 | 0.2132 | O2 higher | 0.564 |
| `ov_coupling_effective_rank_fraction` | 0.3760 | 0.2891 | 0.2703 | 0.1698 | H2O/NH3 higher | 0.542 |
| `same_sign_coeffproduct_edge_fraction` | 0.8264 | 0.8287 | 0.8472 | 0.8542 | O2 higher | 0.435 |
| `excitation_rank1_coefficient_mass_fraction` | 0.0508 | 0.0419 | 0.0250 | 0.0371 | H2O/NH3 higher | 0.414 |
| `coefficient_sign_mass_imbalance` | 0.9002 | 0.9107 | 0.9146 | 0.9392 | O2 higher | 0.382 |
| `coefficient_top1_mass_fraction` | 0.0630 | 0.0648 | 0.1119 | 0.0513 | H2O/NH3 higher | 0.241 |

## Structural interpretation

1. **O2 is much more concentrated on a small set of virtual-orbital hubs.** Its virtual coefficient-mass entropy is 0.851, versus 0.981/0.988 for H2O/NH3; its virtual Gini is 0.390, versus 0.167/0.125; and one virtual orbital carries 28.8% of the weighted incidence, versus 14.5%/12.8%.

2. **O2's occupied-to-virtual parent network is lower-rank/more hub dominated.** The top singular channel carries 85.3% for O2 versus 53.9%/63.3% for H2O/NH3. The normalized effective rank is 0.170, versus 0.376/0.289.

3. **O2 also puts more coefficient mass into rank-2 occupied/virtual transfer parents.** The fraction is 0.601% for O2, compared with 0.144%/0.074% for H2O/NH3.

**Working pre-BCH hypothesis:** signed fermionic parent ordering is more likely to help when significant parent coefficient weight is distributed across many virtual orbitals and many occupied-to-virtual structural channels, rather than being concentrated in a small, low-rank set of hubs. A distributed network gives the parent sequence more structural degrees of freedom to change downstream interference; an O2-like hub-dominated network provides fewer such degrees of freedom.

BeH2 is useful as a check because it is mixed: its virtual-orbital coefficient mass is fairly distributed like H2O/NH3, while its occupied-to-virtual coupling is more concentrated. That is exactly the kind of intermediate structural signature expected for a case whose fermionic advantage is weaker or less robust.

## Caveats / next test

This four-case analysis is hypothesis generation, not a proof. The next validation should apply these same pre-BCH features to the larger HGBS-5 sweep (especially Be2, Li2, F2, H2O, NH3, CH4, N2, and O2) and test whether virtual-hub concentration and O-V effective rank predict the externally measured best-JW/fermionic advantage out of sample.

## Warnings

- NH3_s-1.00_hgbs-5_as-008-010: reconstructed parent count 3907 differs from the reference count 3899 by +8. This can indicate a tensor revision; structural features use the supplied tensor exactly.
