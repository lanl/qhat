# First reproducibility experiment

Independent-repeat checks: PASS
Completed worker jobs: 6/6
Historical metric comparisons flagged: 6

Primary physical error: 1 - abs(overlap). Raw vector error is NOT used:
the existing QHAT exact/Trotter routines use different scalar energy offsets.

Repeat PASS establishes only reproducibility with the CURRENT frozen input/code.
Historical agreement does not establish identical historical tensors or code.
Historical disagreement does not identify its cause. Never replace old results.

## Historical differences

- f2 / analysis/fermionic_body_rank_ablation_20case.csv / fermionic_signed_coefficient_lexicographic / one_minus_overlap: historical=5.8586358e-11, current=1.0581537e-07 (different)
- f2 / analysis/fermionic_body_rank_ablation_20case.csv / fermionic_signed_coefficient_lexicographic / bch2_hf_state_norm: historical=0.003177539, current=0.099855912 (different)
- nh3 / analysis/fermionic_body_rank_ablation_20case.csv / fermionic_signed_coefficient_lexicographic / one_minus_overlap: historical=6.5371342e-10, current=9.6330777e-10 (different)
- nh3 / analysis/fermionic_body_rank_ablation_20case.csv / fermionic_signed_coefficient_lexicographic / bch2_hf_state_norm: historical=0.0077932561, current=0.0091384476 (different)
- nh3 / analysis/fermionic_body_rank_ablation_20case.csv / fermionic_signed_coefficient_lexicographic / number_of_pauli_terms: historical=5064, current=5098 (different)
- nh3 / analysis/fermionic_body_rank_ablation_20case.csv / fermionic_signed_coefficient_lexicographic / number_of_fermionic_terms: historical=2620, current=2628 (different)

## Execution / repeat failures


## Next decision

If repeats fail: investigate environment, numerical tolerances, and code first.
If repeats pass but history differs: trace tensor/orbital/code revisions before pooling data.
Do not start a large new sweep until input lineage and scalar-phase conventions are resolved.
