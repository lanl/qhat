# QHAT Fermionic vs. Pauli Noncommutation Graph Insights

## Purpose

This report summarizes what the generated QHAT bond-length sweep reveals about:

- noncommutation graphs built from Hermitian fermionic terms;
- noncommutation graphs built directly from Jordan–Wigner (JW) Pauli strings; and
- the Trotter error obtained when each graph coloring determines the product-formula order.

The detailed, executed analysis is in [`qhat_L_sweep_trotter_demo.ipynb`](qhat_L_sweep_trotter_demo.ipynb). The numerical results are in [`qhat_full_library_coloring_trotter_errors.csv`](qhat_full_library_coloring_trotter_errors.csv).

## Dataset and scope

The sweep includes all nine molecule families supported by the selected `build_config_L_sweep.py` library:

- H–H
- He–H
- He–He
- Li–H
- Li–Li
- Be–H
- Be–Be
- B–H
- B–B

For each family, the sweep attempted three bond lengths, the STO-6G and HGBS-5 bases, multiple active spaces, and both JW and Bravyi–Kitaev mappings.

Overall generation produced:

- 312 attempted configurations;
- 306 successful `.dat` Hamiltonians;
- 153 unique molecule/bond/basis/active-space cases when JW and BK are treated as mappings of the same case;
- 153 JW files and 153 BK files; and
- 153 matching active-space tensor files.

Six He–He/STO-6G configurations failed because the closed-shell system had no vacant orbital for the requested active space. They are not treated as successful results.

The complete Pauli graph inventory covers all 306 generated Hamiltonians. Dense exact-evolution comparisons cover 108 JW cases at 2, 4, and 6 qubits. The generated 8-qubit cases are included in the graph inventory but not in the dense spectral-norm Trotter comparison.

## What the two graphs represent

### Fermionic graph

Each vertex is a complete Hermitian fermionic term. An edge connects two vertices when their commutator is nonzero:

\[
(T_i,T_j)\in E_f \quad\Longleftrightarrow\quad [T_i,T_j]\ne 0.
\]

The fermionic terms are reconstructed from QHAT's active-space one- and two-body tensors. They are mapped to JW matrices without using OpenFermion in the analysis notebook.

### JW Pauli graph

Each vertex is an individual Pauli string in the final QHAT JW `.dat` Hamiltonian. Two Pauli strings are adjacent when they anticommute:

\[
(P_i,P_j)\in E_P \quad\Longleftrightarrow\quad [P_i,P_j]\ne 0.
\]

For either graph, a color class is an internally commuting group. The notebook uses greedy coloring, so the reported color counts are not guaranteed chromatic numbers.

## Validation

Before comparing Trotter errors, the notebook reconstructs the Hamiltonian from QHAT's fermionic tensors, applies JW, and compares it against the corresponding `.dat` matrix.

The maximum reconstruction discrepancy across the exact-evolution cases was:

\[
4.26\times 10^{-14}.
\]

This confirms that the fermionic and Pauli decompositions represent the same Hamiltonian to numerical precision.

## Main graph result

Direct JW Pauli decomposition generally gives a smaller and less connected noncommutation graph than retaining the original fermionic grouping.

| Active space | Median fermionic vertices | Median Pauli vertices | Median fermionic edges | Median Pauli edges |
|---:|---:|---:|---:|---:|
| 2 qubits | 5 | 4 | 0 | 0 |
| 4 qubits | 23 | 15 | 32 | 16 |
| 6 qubits | 55 | 34 | 360 | 168 |

Across all 108 exact-comparison cases:

- the Pauli graph had fewer vertices in every case;
- the Pauli graph had fewer noncommuting edges in 96 cases;
- the edge counts were equal in 12 cases;
- Pauli coloring used fewer colors in 45 cases;
- both approaches used the same number of colors in 63 cases; and
- fermionic coloring never used fewer colors.

At four qubits, the fermionic graph had approximately twice as many edges as the Pauli graph. At six qubits, it had approximately 2.14 times as many.

The median greedy-color counts were:

| Active space | Fermionic colors | Pauli colors |
|---:|---:|---:|
| 2 qubits | 1 | 1 |
| 4 qubits | 2 | 2 |
| 6 qubits | 4 | 3 |

## Why the Pauli graph can be smaller

JW preserves operator algebra: if two complete fermionic operators commute, their complete JW images commute. However, the two graphs use different vertex definitions.

A Hermitian fermionic term may map to a sum of Pauli strings. When all mapped terms are collected into the full Hamiltonian:

- identical Pauli strings originating from different fermionic terms can combine;
- some contributions can cancel;
- individual Pauli components can reveal finer commuting relationships; and
- treating an entire mapped fermionic term as one Trotter factor hides this finer structure.

Therefore, algebraic commutation is preserved, but graph size and coloring are not preserved when the operator is repartitioned into different vertices.

## Trotter-error result

The analysis compares two product-formula decompositions:

1. color Hermitian fermionic terms, map each complete term to JW, and preserve the fermionic color order;
2. color the final JW Pauli strings directly and use the Pauli color order.

Both are compared with exact time evolution using the spectral operator norm.

### Second order, 10 steps

Among the 96 nontrivial 4- and 6-qubit cases:

- direct JW Pauli coloring had lower error in 88 cases;
- fermionic coloring had lower error in 8 cases; and
- the median fermionic-color error was approximately 1.92 times the Pauli-color error.

| Qubits | Median fermionic-color error | Median Pauli-color error | Median fermionic/Pauli ratio |
|---:|---:|---:|---:|
| 4 | \(3.19\times10^{-6}\) | \(1.69\times10^{-6}\) | 1.92 |
| 6 | \(4.20\times10^{-5}\) | \(3.83\times10^{-5}\) | 1.96 |

At two qubits, the relevant terms commute, and measured errors are numerical roundoff near \(10^{-15}\).

### First order

The first-order median errors were essentially equal. Individual cases favored different orderings, but there was no comparably strong systematic winner. The advantage of direct Pauli coloring was much clearer for the symmetric second-order formula.

## Trends by molecule

For second-order Trotterization at 10 steps, the median ratio

\[
\frac{\text{fermionic-color error}}{\text{Pauli-color error}}
\]

was:

| Molecule | Median ratio |
|---|---:|
| He–H | 1.42 |
| B–B | 1.49 |
| He–He | 1.58 |
| Be–Be | 1.59 |
| Be–H | 1.76 |
| B–H | 1.90 |
| Li–Li | 1.94 |
| H–H | 1.99 |
| Li–H | 2.08 |

A value above one means the direct Pauli-color ordering had the lower median error.

## Active-space scaling

Graph complexity rises rapidly with active-space size. The median fermionic edge count increased from 32 at four qubits to 360 at six qubits. The median Pauli edge count increased from 16 to 168.

Larger active spaces introduce:

- more Hamiltonian terms;
- more possible noncommuting pairs;
- more commuting groups/colors; and
- generally larger Trotter errors at a fixed number of steps.

In this dataset, graph size, edge count, and color count are positively correlated with Trotter error. This correlation should not be interpreted as proof that graph density alone determines error, because all of these quantities also grow with system size.

## Main interpretation

The strongest empirical conclusion is:

> Commutation at the fermionic-term level is preserved by JW, but retaining complete fermionic terms as indivisible Trotter factors hides finer Pauli-level commutation, combination, and cancellation structure.

Consequently:

- fermionic coloring gives a chemically meaningful operator decomposition;
- direct Pauli coloring usually produces fewer graph edges and sometimes fewer colors;
- direct Pauli coloring generally gives lower second-order Trotter error in this dataset; and
- neither fewer colors nor fewer edges alone guarantees lower error in every individual case.

The numerical coefficients and nested commutators between color blocks still matter. Graph coloring is therefore a useful ordering heuristic, not a complete error model.

## Limitations and cautions

- The coloring algorithm is greedy and does not prove the minimum possible number of colors.
- Dense exact-evolution comparisons stop at six qubits; eight-qubit cases currently have graph results only.
- Only JW is used for the fermionic-versus-Pauli Trotter comparison. BK files are included in the overall Pauli graph inventory.
- Results depend on the selected bases, active spaces, bond-length range, coefficient threshold, term grouping, total evolution time, and Trotter formula.
- A different ordering of the same color blocks can change product-formula error even when the graph coloring is unchanged.
- Lower graph density does not by itself prove a better product formula; commutator magnitudes are also important.

## Suggested follow-up studies

Useful next investigations include:

1. compare several greedy-coloring strategies and exact coloring for small cases;
2. reorder color blocks using weighted commutator norms rather than color identifiers;
3. compare JW and BK Pauli colorings directly;
4. compute sparse or matrix-free errors for the generated eight-qubit cases;
5. track trends separately by bond length and basis;
6. compare unweighted graph metrics with coefficient-weighted noncommutation graphs; and
7. test whether nested-commutator metrics predict the observed second-order error more accurately than edge or color counts.
