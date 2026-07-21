# Focused Fermionic-vs-JW Trotter Benchmark

## Benchmark definition

The benchmark uses 100 seeded random permutations of each decomposition's nonidentity vertices. The colored baseline flattens deterministic largest-degree-first color groups. The same sampled permutations are reused across first/second order and all step counts.

- Evolution time: `1.0`
- Steps: `[1, 2, 5, 10]`
- Implemented orders: first and second only
- Random seed: `20260721`
- ‘Best observed’ means best among the colored baseline and sampled orders; it is not a proven optimum.

Fermionic tensor monomials are first put into the same descending-index normal-order convention used by `h2_fermionic.ipynb`, and algebraically equivalent permutations are combined before Hermitian pairing. This keeps the comparison independent of OpenFermion without changing its term convention.

For color block `B_c = sum(H_i for i in color c)`, the report records pairwise `||[B_c,B_d]||₂`, the norm of the ordered commutator sum as a first-order diagnostic, and the sum of left/right nested-commutator norms as a second-order diagnostic. These are coefficient-weighted explanatory proxies, not rigorous rankings of the full finite-time error.

## Selected cases

| case_id | molecule | bond_length | basis | active_occupied | active_vacant | qubits | selection_reason |
| --- | --- | --- | --- | --- | --- | --- | --- |
| similar_4q | Be-Be | 2.244 | sto-6g | 2 | 2 | 4 | 4-qubit case where both colored methods have similar error |
| fermionic_first_order_6q | B-B | 2.38 | sto-6g | 2 | 4 | 6 | 6-qubit case where fermionic coloring is better at first order |
| jw_second_order_6q | B-B | 2.38 | hgbs-5 | 2 | 4 | 6 | 6-qubit case where direct JW coloring is much better at second order |

## similar_4q: Be-Be at L=2.244

4-qubit case where both colored methods have similar error

### Graph and coefficient-weighted diagnostics

| scheme | vertices | noncommuting_edges | colors | block_commutator_norm_sum | ordered_commutator_sum_norm | nested_commutator_norm_sum |
| --- | --- | --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 13 | 12 | 2 | 0.0183315 | 0.0183315 | 0.00702131 |
| Direct JW Pauli coloring | 15 | 16 | 2 | 0.0183315 | 0.0183315 | 0.00702131 |

The color order and vertex membership used by the baseline are:

| scheme | color_order_json | color_groups_json |
| --- | --- | --- |
| Fermionic coloring → JW | [0,1] | {"0": [0, 6, 7, 8, 11], "1": [1, 2, 3, 4, 5, 9, 10, 12]} |
| Direct JW Pauli coloring | [0,1] | {"0": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], "1": [11, 12, 13, 14]} |

Pairwise color-block diagnostics:

| scheme | left_color | right_color | commutator_spectral_norm | normalized_commutator_norm | left_nested_commutator_norm | right_nested_commutator_norm |
| --- | --- | --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 0 | 1 | 0.0183315 | 0.000370942 | 0.00357641 | 0.00344489 |
| Direct JW Pauli coloring | 0 | 1 | 0.0183315 | 0.00323914 | 0.00344489 | 0.00357641 |

### Colored baseline errors

| scheme | formula_order | steps | colored_error |
| --- | --- | --- | --- |
| Fermionic coloring → JW | 1 | 1 | 0.00914708 |
| Fermionic coloring → JW | 1 | 2 | 0.00457004 |
| Fermionic coloring → JW | 1 | 5 | 0.00182762 |
| Fermionic coloring → JW | 1 | 10 | 0.000913784 |
| Fermionic coloring → JW | 2 | 1 | 0.000322947 |
| Fermionic coloring → JW | 2 | 2 | 8.06642e-05 |
| Fermionic coloring → JW | 2 | 5 | 1.2903e-05 |
| Fermionic coloring → JW | 2 | 10 | 3.22564e-06 |
| Direct JW Pauli coloring | 1 | 1 | 0.00914708 |
| Direct JW Pauli coloring | 1 | 2 | 0.00457004 |
| Direct JW Pauli coloring | 1 | 5 | 0.00182762 |
| Direct JW Pauli coloring | 1 | 10 | 0.000913784 |
| Direct JW Pauli coloring | 2 | 1 | 0.000330279 |
| Direct JW Pauli coloring | 2 | 2 | 8.24962e-05 |
| Direct JW Pauli coloring | 2 | 5 | 1.31961e-05 |
| Direct JW Pauli coloring | 2 | 10 | 3.29891e-06 |

### Ordering robustness

| scheme | formula_order | colored_error | random_best_observed_error | random_median_error | random_p90_error | colored_over_best_observed | scheme_best_over_global_best | colored_over_global_best | fraction_random_better_than_colored |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 1 | 0.000913784 | 0.000136987 | 0.00383981 | 0.00774683 | 6.67061 | 1 | 6.67061 | 0.2 |
| Fermionic coloring → JW | 2 | 3.22564e-06 | 4.09577e-07 | 2.11419e-05 | 8.02321e-05 | 7.87552 | 1.07365 | 8.45553 | 0.1 |
| Direct JW Pauli coloring | 1 | 0.000913784 | 0.000169673 | 0.000458109 | 0.000643954 | 5.38555 | 1.23861 | 6.67061 | 0.96 |
| Direct JW Pauli coloring | 2 | 3.29891e-06 | 3.81483e-07 | 1.2198e-06 | 2.4163e-06 | 8.64759 | 1 | 8.64759 | 1 |

### Particle-number leakage

| scheme | formula_order | colored_particle_number_leakage | random_median_particle_number_leakage | random_worst_particle_number_leakage |
| --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 1 | 2.44249e-15 | 2.88658e-15 | 5.55112e-15 |
| Fermionic coloring → JW | 2 | 0 | 0 | 5.55112e-16 |
| Direct JW Pauli coloring | 1 | 0 | 2.22045e-15 | 9.65894e-15 |
| Direct JW Pauli coloring | 2 | 0 | 0 | 4.996e-15 |

## fermionic_first_order_6q: B-B at L=2.38

6-qubit case where fermionic coloring is better at first order

### Graph and coefficient-weighted diagnostics

| scheme | vertices | noncommuting_edges | colors | block_commutator_norm_sum | ordered_commutator_sum_norm | nested_commutator_norm_sum |
| --- | --- | --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 28 | 90 | 4 | 0.0939502 | 0.0471468 | 0.166109 |
| Direct JW Pauli coloring | 34 | 168 | 3 | 0.118061 | 0.0579512 | 0.265016 |

The color order and vertex membership used by the baseline are:

| scheme | color_order_json | color_groups_json |
| --- | --- | --- |
| Fermionic coloring → JW | [0,1,2,3] | {"0": [0, 5, 6, 8, 10, 11, 17, 27], "1": [3, 4, 9, 13, 14, 19, 20], "2": [1, 2, 7, 21, 22, 23, 26], "3": [12, 15, 16, 18, 24, 25]} |
| Direct JW Pauli coloring | [0,1,2] | {"0": [0, 5, 6, 7, 8, 9, 12, 13, 16, 21, 22, 23, 24, 25], "1": [3, 4, 10, 11, 14, 15, 26, 27, 28, 29], "2": [1, 2, 17, 18, 19, 20, 30, 31, 32, 33]} |

Pairwise color-block diagnostics:

| scheme | left_color | right_color | commutator_spectral_norm | normalized_commutator_norm | left_nested_commutator_norm | right_nested_commutator_norm |
| --- | --- | --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 0 | 1 | 0.0315619 | 0.00047337 | 0.0269385 | 0.0269385 |
| Fermionic coloring → JW | 0 | 2 | 0.0311463 | 0.000461159 | 0.0225262 | 0.0335881 |
| Fermionic coloring → JW | 0 | 3 | 4.79014e-05 | 2.26753e-07 | 1.80065e-06 | 1.2208e-07 |
| Fermionic coloring → JW | 1 | 2 | 0.0311463 | 0.0327626 | 0.0225262 | 0.0335881 |
| Fermionic coloring → JW | 1 | 3 | 4.79014e-05 | 1.61094e-05 | 1.80065e-06 | 1.2208e-07 |
| Fermionic coloring → JW | 2 | 3 | 5.20417e-18 | 1.72779e-18 | 1.62691e-19 | 1.54074e-33 |
| Direct JW Pauli coloring | 0 | 1 | 0.0429253 | 0.000532645 | 0.0496067 | 0.0496067 |
| Direct JW Pauli coloring | 0 | 2 | 0.0375681 | 0.000509606 | 0.0413035 | 0.0415978 |
| Direct JW Pauli coloring | 1 | 2 | 0.0375681 | 0.0300806 | 0.0413035 | 0.0415978 |

### Colored baseline errors

| scheme | formula_order | steps | colored_error |
| --- | --- | --- | --- |
| Fermionic coloring → JW | 1 | 1 | 0.0220589 |
| Fermionic coloring → JW | 1 | 2 | 0.0115865 |
| Fermionic coloring → JW | 1 | 5 | 0.00469809 |
| Fermionic coloring → JW | 1 | 10 | 0.00235359 |
| Fermionic coloring → JW | 2 | 1 | 0.00400425 |
| Fermionic coloring → JW | 2 | 2 | 0.00102894 |
| Fermionic coloring → JW | 2 | 5 | 0.000165896 |
| Fermionic coloring → JW | 2 | 10 | 4.15193e-05 |
| Direct JW Pauli coloring | 1 | 1 | 0.0261914 |
| Direct JW Pauli coloring | 1 | 2 | 0.0141233 |
| Direct JW Pauli coloring | 1 | 5 | 0.00576771 |
| Direct JW Pauli coloring | 1 | 10 | 0.00289236 |
| Direct JW Pauli coloring | 2 | 1 | 0.00491341 |
| Direct JW Pauli coloring | 2 | 2 | 0.00127083 |
| Direct JW Pauli coloring | 2 | 5 | 0.000205269 |
| Direct JW Pauli coloring | 2 | 10 | 5.13866e-05 |

### Ordering robustness

| scheme | formula_order | colored_error | random_best_observed_error | random_median_error | random_p90_error | colored_over_best_observed | scheme_best_over_global_best | colored_over_global_best | fraction_random_better_than_colored |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 1 | 0.00235359 | 0.000711566 | 0.00207881 | 0.00309533 | 3.30762 | 1 | 3.30762 | 0.63 |
| Fermionic coloring → JW | 2 | 4.15193e-05 | 3.16406e-06 | 2.7681e-05 | 5.49229e-05 | 13.1221 | 1 | 13.1221 | 0.77 |
| Direct JW Pauli coloring | 1 | 0.00289236 | 0.00106728 | 0.00163185 | 0.00212629 | 2.71002 | 1.49991 | 4.06479 | 1 |
| Direct JW Pauli coloring | 2 | 5.13866e-05 | 1.80394e-05 | 3.83541e-05 | 5.86027e-05 | 2.84858 | 5.70133 | 16.2407 | 0.79 |

### Particle-number leakage

| scheme | formula_order | colored_particle_number_leakage | random_median_particle_number_leakage | random_worst_particle_number_leakage |
| --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 1 | 0 | 0 | 0 |
| Fermionic coloring → JW | 2 | 0 | 0 | 6.66134e-16 |
| Direct JW Pauli coloring | 1 | 1.52101e-14 | 1.13669e-07 | 1.36352e-06 |
| Direct JW Pauli coloring | 2 | 8.32667e-15 | 1.78715e-10 | 2.68754e-09 |

## jw_second_order_6q: B-B at L=2.38

6-qubit case where direct JW coloring is much better at second order

### Graph and coefficient-weighted diagnostics

| scheme | vertices | noncommuting_edges | colors | block_commutator_norm_sum | ordered_commutator_sum_norm | nested_commutator_norm_sum |
| --- | --- | --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 28 | 90 | 4 | 0.0382388 | 0.0175823 | 0.0595734 |
| Direct JW Pauli coloring | 34 | 168 | 3 | 0.0156569 | 0.00727631 | 0.0110153 |

The color order and vertex membership used by the baseline are:

| scheme | color_order_json | color_groups_json |
| --- | --- | --- |
| Fermionic coloring → JW | [0,1,2,3] | {"0": [0, 5, 6, 8, 10, 11, 17, 27], "1": [3, 4, 9, 13, 14, 19, 20], "2": [1, 2, 7, 21, 22, 23, 26], "3": [12, 15, 16, 18, 24, 25]} |
| Direct JW Pauli coloring | [0,1,2] | {"0": [0, 5, 6, 7, 8, 9, 12, 13, 16, 21, 22, 23, 24, 25], "1": [3, 4, 10, 11, 14, 15, 26, 27, 28, 29], "2": [1, 2, 17, 18, 19, 20, 30, 31, 32, 33]} |

Pairwise color-block diagnostics:

| scheme | left_color | right_color | commutator_spectral_norm | normalized_commutator_norm | left_nested_commutator_norm | right_nested_commutator_norm |
| --- | --- | --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 0 | 1 | 0.0123686 | 0.000185322 | 0.00966731 | 0.00966731 |
| Fermionic coloring → JW | 0 | 2 | 0.0123027 | 0.000146815 | 0.00862086 | 0.0114383 |
| Fermionic coloring → JW | 0 | 3 | 0.000632314 | 5.23201e-06 | 1.00805e-05 | 5.01588e-05 |
| Fermionic coloring → JW | 1 | 2 | 0.0123027 | 0.010514 | 0.00862086 | 0.0114383 |
| Fermionic coloring → JW | 1 | 3 | 0.000632314 | 0.000374683 | 1.00805e-05 | 5.01588e-05 |
| Fermionic coloring → JW | 2 | 3 | 3.72966e-17 | 1.76021e-17 | 5.29369e-19 | 1.94134e-31 |
| Direct JW Pauli coloring | 0 | 1 | 0.00606966 | 0.000223194 | 0.00234657 | 0.00234657 |
| Direct JW Pauli coloring | 0 | 2 | 0.00479361 | 0.000259002 | 0.0020817 | 0.00107941 |
| Direct JW Pauli coloring | 1 | 2 | 0.00479361 | 0.0458411 | 0.0020817 | 0.00107941 |

### Colored baseline errors

| scheme | formula_order | steps | colored_error |
| --- | --- | --- | --- |
| Fermionic coloring → JW | 1 | 1 | 0.00835222 |
| Fermionic coloring → JW | 1 | 2 | 0.00433998 |
| Fermionic coloring → JW | 1 | 5 | 0.00175461 |
| Fermionic coloring → JW | 1 | 10 | 0.000878642 |
| Fermionic coloring → JW | 2 | 1 | 0.00120826 |
| Fermionic coloring → JW | 2 | 2 | 0.000308021 |
| Fermionic coloring → JW | 2 | 5 | 4.95527e-05 |
| Fermionic coloring → JW | 2 | 10 | 1.23978e-05 |
| Direct JW Pauli coloring | 1 | 1 | 0.003607 |
| Direct JW Pauli coloring | 1 | 2 | 0.00181352 |
| Direct JW Pauli coloring | 1 | 5 | 0.00072652 |
| Direct JW Pauli coloring | 1 | 10 | 0.000363337 |
| Direct JW Pauli coloring | 2 | 1 | 0.00012377 |
| Direct JW Pauli coloring | 2 | 2 | 3.09796e-05 |
| Direct JW Pauli coloring | 2 | 5 | 4.9584e-06 |
| Direct JW Pauli coloring | 2 | 10 | 1.23966e-06 |

### Ordering robustness

| scheme | formula_order | colored_error | random_best_observed_error | random_median_error | random_p90_error | colored_over_best_observed | scheme_best_over_global_best | colored_over_global_best | fraction_random_better_than_colored |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 1 | 0.000878642 | 0.000243566 | 0.000729663 | 0.00108606 | 3.60741 | 1.54775 | 5.58338 | 0.71 |
| Fermionic coloring → JW | 2 | 1.23978e-05 | 7.48332e-07 | 6.82414e-06 | 1.53783e-05 | 16.5673 | 1.09614 | 18.16 | 0.8 |
| Direct JW Pauli coloring | 1 | 0.000363337 | 0.000157367 | 0.000235524 | 0.000299506 | 2.30884 | 1 | 2.30884 | 0.99 |
| Direct JW Pauli coloring | 2 | 1.23966e-06 | 6.827e-07 | 2.03102e-06 | 3.42349e-06 | 1.81582 | 1 | 1.81582 | 0.13 |

### Particle-number leakage

| scheme | formula_order | colored_particle_number_leakage | random_median_particle_number_leakage | random_worst_particle_number_leakage |
| --- | --- | --- | --- | --- |
| Fermionic coloring → JW | 1 | 0 | 0 | 3.66374e-15 |
| Fermionic coloring → JW | 2 | 1.66533e-15 | 0 | 5.77316e-15 |
| Direct JW Pauli coloring | 1 | 1.9984e-15 | 5.47944e-09 | 8.11316e-08 |
| Direct JW Pauli coloring | 2 | 5.9952e-15 | 6.32439e-13 | 1.7396e-11 |

## Interpretation

The B–B STO-6G and HGBS-5 cases have the same molecule, geometry, active space, graph counts, and deterministic color memberships. Nevertheless, the preferred decomposition reverses when the basis changes. This is evidence that coefficient-weighted block commutators—and for second order, nested commutators—capture behavior that graph topology alone misses. The three-case study supports this explanation but does not establish a universal predictor.

Fermionic vertices preserve particle number individually. Random JW Pauli-string orders can temporarily separate cancellation partners and leak out of the target particle-number sector. The colored JW baseline may suppress that leakage when the relevant strings remain grouped. In these cases the leakage is much smaller than the spectral error, so it is a physical robustness diagnostic rather than the primary explanation for the operator-norm error.

The random-ordering ensemble also shows whether a colored result is robust or merely one favorable schedule. Compare the colored error with the random median and best observed error rather than treating one deterministic ordering as definitive.

Canonical normal-order aggregation is itself important: compared with the earlier raw tensor-slot convention, it reduces the selected four-qubit fermionic graph from 23 vertices/32 edges to 13/12 and each selected six-qubit fermionic graph from 55/360 to 28/90. The colored errors remain unchanged in these three cases because the combined permutations occupied the same commuting blocks, while the random-order statistics change because the sampling unit is now a genuine canonical fermionic term. This is another concrete reason that raw graph size is representation-dependent.

## Scope limits

- No fourth-order product formula is included.
- No dense all-case 8-qubit sweep is included.
- The random ensemble samples vertex orders, not additional coloring algorithms.
- Best observed errors are sampled results, not mathematical optima.
- Nested-commutator sums are pairwise heuristics, not the complete schedule-weighted BCH error generator.
- Normalized block-commutator values include full block norms and are sensitive to scalar energy shifts; the absolute norms are the primary diagnostics.
- No gate-count or circuit-depth claim is made.
