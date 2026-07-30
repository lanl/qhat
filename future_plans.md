# Future Plans

## Project Goal

This project studies whether graph-coloring information from fermionic Hamiltonian terms can be used to construct effective Pauli-string orderings for Trotterized quantum simulation.

The current implementation compares several ordering strategies for the Jordan–Wigner Hamiltonian, including:

* Raw Jordan–Wigner insertion order
* Jordan–Wigner Pauli-string graph coloring
* Fermionic-term graph coloring mapped to a Jordan–Wigner Pauli-string order
* Reversed color-group order
* Random color-group order
* Random within-group order
* Random group and within-group order

The fermionic method should be described as a **fermionic-coloring-induced Jordan–Wigner Pauli ordering**. The current code does not exponentiate each complete fermionic operator as a single block.

## Current Status

The benchmark is currently running across a large collection of molecular Hamiltonians and active-space configurations.

For each Hamiltonian, the code records:

* Molecular and active-space metadata
* Number of qubits
* Number of fermionic and Pauli terms
* Graph size, density, and number of colors
* Color-group sizes and execution order
* Random seeds and ordering hashes
* Trotter step count and time step
* BCH-based error indicators
* Exact symmetry-sector evolution
* State overlap and infidelity
* Phase-aligned state-vector error
* Particle-number and spin-sector leakage
* Runtime information

The remaining cases should be completed before drawing final conclusions.

## Immediate Next Steps

* [ ] Complete the remaining benchmark cases.
* [ ] Confirm that every expected Hamiltonian case completed successfully.
* [ ] Record and investigate failed, skipped, or incomplete cases.
* [ ] Verify that no output rows contain unexpected `NaN`, `inf`, or missing values.
* [ ] Preserve the raw benchmark output before performing aggregation or filtering.
* [ ] Combine all per-case outputs into one master dataset.

## Implementation Validation

Before using the complete dataset for scientific conclusions, perform the following validation checks.

### 1. Hamiltonian Reconstruction

Verify that the Hermitian fermionic terms reconstruct the original fermionic Hamiltonian within the selected numerical tolerance.

### 2. Complete Pauli Permutation

For every ordering strategy, verify that the generated Pauli order is a complete permutation of the final combined Jordan–Wigner Hamiltonian:

* No missing Pauli strings
* No duplicated Pauli strings
* No extra Pauli strings
* Identical coefficients across ordering methods

### 3. Dense and Matrix-Free Agreement

For representative small cases, compare the matrix-free implementation against an independently constructed dense-matrix calculation using the same:

* Hamiltonian
* Pauli-string order
* Evolution time
* Trotter formula
* Trotter step count

The state infidelity and phase-aligned state-vector error should agree to numerical precision.

### 4. Trotter Convergence

For representative cases, repeat the calculation at multiple Trotter step counts.

For first-order Trotterization, verify that the leading state-vector error decreases approximately as

[
O!\left(\frac{1}{r}\right),
]

while the state infidelity decreases approximately as

[
O!\left(\frac{1}{r^2}\right),
]

before reaching the numerical-precision floor.

### 5. Symmetry Preservation

Check that the following remain near zero:

* Particle-number leakage
* Spin-sector leakage

Cases with unexpectedly large leakage should be inspected separately.

### 6. Commuting-Group Invariance

For a valid Pauli coloring, changing the order of Pauli strings inside one commuting group should not materially change the Trotter result.

This should be tested as an implementation sanity check. Small differences at machine precision are acceptable.

## Post-Processing Plan

After all cases finish, aggregate the results by:

* Molecule
* Bond length
* Basis
* Active occupied orbitals
* Active vacant orbitals
* Number of qubits
* Graph level
* Ordering schedule
* Trotter step count

For randomized schedules, calculate:

* Mean
* Median
* Standard deviation
* Minimum
* Maximum
* Interquartile range

The primary performance metric should be `state_infidelity`.

Supporting metrics should include:

* `bch2_hf_state_norm`
* `leading_state_error_norm_estimate`
* `phase_aligned_state_2norm_error`
* Particle-number leakage
* Spin-sector leakage
* Runtime

The unaligned state-vector norm should not be interpreted by itself because it can be dominated by an irrelevant global phase.

## Main Scientific Comparisons

The main comparison will be between:

1. Raw Jordan–Wigner ordering
2. Jordan–Wigner Pauli coloring
3. Fermionic-coloring-induced Jordan–Wigner ordering

For each Hamiltonian, determine whether the fermionic-induced ordering:

* Outperforms the raw Jordan–Wigner order
* Outperforms the Jordan–Wigner coloring order
* Performs similarly to the Pauli-based methods
* Produces highly variable results across randomized color-group orders
* Performs worse than the Pauli-based methods

## Baselines

The final study should include the following baselines when available:

* Raw insertion order
* Magnitude ordering
* Lexicographical ordering
* Fully random Pauli ordering
* Jordan–Wigner graph coloring
* Fermionic-induced graph coloring
* Tranter ordering or another chemistry-informed ordering method

All methods must use the same Hamiltonian coefficients, initial state, evolution time, Trotter formula, and number of Trotter steps.

## Questions to Investigate

The final analysis should address the following questions:

1. Does fermionic coloring reduce Trotter error more often than Jordan–Wigner Pauli coloring?
2. Is performance determined mainly by the number of colors, or by the order of the color groups?
3. Does a smaller fermionic graph lead to a better Pauli ordering?
4. How strongly does the result depend on randomized group order?
5. Does the BCH Hartree–Fock state norm predict the measured state infidelity?
6. When does fermionic coloring produce cancellation among coefficient-weighted Pauli commutators?
7. Which molecular, active-space, or graph features are associated with favorable fermionic-induced orderings?
8. Are there cases where the BCH metric and measured infidelity rank the orderings differently?

## Planned Figures

The final report should include:

* Infidelity versus Trotter steps
* BCH Hartree–Fock state norm versus Trotter steps
* Fermionic-induced infidelity versus Jordan–Wigner-coloring infidelity
* Distribution of randomized color-group results
* Win-rate summary across Hamiltonians
* Performance grouped by number of qubits
* Performance grouped by graph density and number of colors
* Best, median, and worst representative cases

For randomized experiments, plots should show the median or mean together with a variability interval.

## Interpretation Guidelines

The following claims should be avoided unless supported by the complete dataset:

* Fewer colors always imply lower Trotter error.
* A smaller graph always produces a better Trotter ordering.
* Fermionic coloring is universally better than Pauli coloring.
* A lower BCH norm guarantees lower finite-time infidelity.

A more appropriate working hypothesis is:

> Fermionic-term coloring can produce useful Jordan–Wigner Pauli orderings when the induced sequence improves cancellation among signed, coefficient-weighted commutators, especially those that act strongly on the selected initial state.

## Possible Method Extensions

After validating the current first-order benchmark, future extensions may include:

* Second-order Suzuki–Trotter formulas
* Multiple evolution times
* Additional initial states
* Optimized color-group ordering
* Simulated annealing over group order
* Coefficient-weighted graph objectives
* Commutator-weighted graph objectives
* Direct exponentiation of complete mapped fermionic blocks
* Jordan–Wigner and Bravyi–Kitaev comparisons under matched induced orderings
* Larger active spaces and additional molecular families

Direct fermionic-block exponentiation should be treated as a separate method from the current fermionic-induced Pauli ordering.

## Reproducibility

Before sharing final results:

* [ ] Store the command used for each benchmark run.
* [ ] Record the Git commit hash.
* [ ] Record the Python and dependency versions.
* [ ] Preserve all random seeds.
* [ ] Preserve ordering hashes.
* [ ] Document the coefficient tolerance.
* [ ] Keep raw outputs separate from processed summaries.
* [ ] Add a small reproducible example to the repository.
* [ ] Add automated tests for the ordering and evolution kernels.

## Expected Final Deliverables

The completed project should produce:

* A validated benchmark dataset
* A summary CSV with aggregated statistics
* A set of reproducible plotting scripts
* Representative case studies
* A comparison of fermionic-induced and Pauli-based orderings
* A documented validation procedure
* A concise statement describing when fermionic coloring appears useful
