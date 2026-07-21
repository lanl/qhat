# Reproducing the QHAT L-Sweep and Trotter Analysis

## Overview

This guide explains how to reproduce two related datasets:

1. QHAT molecular Hamiltonian files, including fermionic tensors and mapped Pauli-string `.dat` files.
2. The analyzer CSV comparing Trotter errors from fermionic-term coloring and direct JW Pauli-string coloring.

The executed reference notebook is [`qhat_L_sweep_trotter_demo.ipynb`](qhat_L_sweep_trotter_demo.ipynb). Its saved analyzer output is [`qhat_full_library_coloring_trotter_errors.csv`](qhat_full_library_coloring_trotter_errors.csv), and the findings are summarized in [`QHAT_FERMIONIC_VS_PAULI_NONCOMMUTATION_INSIGHTS.md`](QHAT_FERMIONIC_VS_PAULI_NONCOMMUTATION_INSIGHTS.md).

## Requirements

Use a Python environment containing the QHAT Hamiltonian-generation dependencies, including:

- `basis-set-exchange`
- `mendeleev`
- `numpy`
- `openfermionpyscf`
- `pyscf`
- `pandas`
- `networkx`
- `scipy`
- `matplotlib`
- Jupyter and `nbconvert`

For example, if the repository's Conda environment is available:

```bash
conda activate qhat
```

## 1. Generate the QHAT molecular library

From the repository root, enter the Hamiltonian generator:

```bash
cd hamiltonian_generator
```

Run the complete nine-family sweep:

```bash
python build_config_L_sweep.py \
  H-H,He-H,He-He,Li-H,Li-Li,Be-H,Be-Be,B-H,B-B \
  --L-steps 3 \
  --L-min-frac 0.8 \
  --L-max-frac 1.4 \
  --max-active 8 \
  --library ../demo_L_sweep_library \
  --hamgen-dir . \
  --run
```

This command creates configurations and immediately runs `hamgen.py` on each one.

### Molecule families

The selected sweep contains:

- H–H
- He–H
- He–He
- Li–H
- Li–Li
- Be–H
- Be–Be
- B–H
- B–B

The builder supports the `X-H` and `X-X` families. It does not generate arbitrary polyatomic geometries.

### Bond lengths

For atoms with covalent radii \(r_1\) and \(r_2\), the builder estimates

\[
L_{\mathrm{eq}}=r_1+r_2.
\]

It then generates evenly spaced bond lengths between

\[
L_{\min}=0.8L_{\mathrm{eq}}
\]

and

\[
L_{\max}=1.4L_{\mathrm{eq}}.
\]

With `--L-steps 3`, each molecule receives three bond lengths. Increase this value for a denser bond scan.

### Basis and active-space variants

The builder automatically considers:

- STO-6G;
- HGBS-5;
- valid active occupied/vacant spaces up to `--max-active`; and
- JW and BK mappings.

For example,

```text
as-002-004
```

means two active occupied spin orbitals and four active vacant spin orbitals, giving six active spin orbitals/qubits.

Not every requested combination is physically or numerically valid. In the reference run, six He–He/STO-6G configurations failed because the closed-shell system did not have a vacant orbital for the requested active space.

## 2. Understand the generated files

A typical output directory is:

```text
demo_L_sweep_library/
└── Li-Li/
    └── 2.93/
        └── hgbs-5/
            ├── Li-Li_2.93_hgbs-5_as-002-004_JW.config
            ├── Li-Li_2.93_hgbs-5_as-002-004_jw.dat
            ├── Li-Li_2.93_hgbs-5_as-002-004_bk.dat
            ├── Li-Li_2.93_hgbs-5_as-002-004.tensors.npz
            └── ...
```

Important file types are:

- `.config`: input to `hamgen.py`;
- `.tensors.npz`: the active-space constant, one-body tensor, and two-body tensor;
- `_jw.dat`: the JW Pauli Hamiltonian;
- `_bk.dat`: the BK Pauli Hamiltonian;
- `.npy`: the generated initial state;
- `.pickle` and `.hdf5`: reusable but potentially very large intermediate calculations; and
- `.log`: generation diagnostics.

The generation pipeline is:

```text
molecular geometry
    → Hartree–Fock calculation
    → active-space Hamiltonian
    → fermionic coefficient tensors
    → JW/BK mapping
    → Pauli-string .dat file
```

A `.dat` file contains coefficients and Pauli strings:

```text
-1.47260413107609658e+01 IIII
 8.58532520960764434e-02 ZIII
-3.88057932231605741e-03 XXYY
```

Each row represents \(c_iP_i\), where \(c_i\) is the coefficient and \(P_i\) is a Pauli string.

## 3. Build the fermionic noncommutation graph

The notebook reads each `.tensors.npz` file and reconstructs terms of the form

\[
h_{pq}a_p^\dagger a_q
\]

and

\[
h_{pqrs}a_p^\dagger a_q^\dagger a_r a_s.
\]

Before graph construction, tensor permutations are canonicalized into the same
descending-index normal-order convention used by OpenFermion in
`h2_fermionic.ipynb`. Fermionic antisymmetry signs are applied, repeated
equal-action operators are discarded, and equivalent monomials are combined.
The implementation performs this limited tensor canonicalization directly and
does not import OpenFermion.

A monomial and its Hermitian conjugate are grouped into one Hermitian vertex:

\[
T_i=c_iO_i+c_i^*O_i^\dagger.
\]

Two vertices are adjacent when

\[
[T_i,T_j]\ne0.
\]

The notebook constructs JW matrices for these complete fermionic vertices without importing OpenFermion. Since JW preserves commutators, testing the mapped matrix commutator is equivalent to testing the original fermionic commutator.

## 4. Build the Pauli-string noncommutation graph

For the direct JW graph, every individual Pauli string in `_jw.dat` is a vertex.

Two Pauli strings anticommute when the number of positions at which both are nonidentity and different is odd. Therefore:

```text
XI and ZI → anticommute
XI and IZ → commute
XX and ZZ → commute
```

An edge is added for every anticommuting pair. The notebook colors both graphs with NetworkX's greedy `largest_first` strategy. Because edges represent noncommutation, every color class is an internally commuting group.

The greedy color count is not guaranteed to equal the graph's chromatic number.

## 5. Validate the tensor reconstruction

Before calculating Trotter error, the notebook verifies that the tensor reconstruction and QHAT `.dat` file represent the same Hamiltonian:

\[
H_{\mathrm{fermionic}\rightarrow\mathrm{JW}}
\approx
H_{\mathrm{QHAT\ .dat}}.
\]

The maximum discrepancy in the reference results was approximately

\[
4.26\times10^{-14},
\]

which is consistent with floating-point roundoff.

## 6. Construct the product formulas

The analysis compares two decompositions:

```text
Hermitian fermionic terms
    → fermionic noncommutation coloring
    → map complete colored terms to JW
    → preserve the fermionic color order
```

and

```text
QHAT JW Pauli strings
    → Pauli noncommutation coloring
    → preserve the Pauli color order
```

### First order

For color blocks \(H_1,\ldots,H_m\), one first-order step is

\[
S_1(\Delta t)
=e^{-iH_m\Delta t}\cdots e^{-iH_1\Delta t}.
\]

For \(r\) steps,

\[
U_1(t,r)=S_1(t/r)^r.
\]

### Second order

The symmetric second-order step uses a forward and reverse sequence:

\[
S_2(\Delta t)
=e^{-iH_1\Delta t/2}\cdots e^{-iH_m\Delta t/2}
e^{-iH_m\Delta t/2}\cdots e^{-iH_1\Delta t/2}.
\]

The current analyzer implements first and second order only. It does not calculate fourth-order error.

## 7. Calculate Trotter error

The exact evolution is

\[
U_{\mathrm{exact}}=e^{-iHt}.
\]

The reported error is the spectral operator norm

\[
\epsilon
=\left\|U_{\mathrm{Trotter}}-U_{\mathrm{exact}}\right\|_2.
\]

The reference settings are:

```python
CASE_MAX_QUBITS = 6
CASE_ORDERS = (1, 2)
CASE_STEPS = (1, 2, 5, 10)
CASE_TOTAL_TIME = 1.0
```

The six-qubit limit applies only to dense exact evolution. The notebook can still inventory and color larger generated Hamiltonians.

## 8. Run the analyzer notebook

If the molecular library has already been generated, keep generation disabled and enable the case-study runner:

```python
RUN_GENERATION = False
RUN_ALL_CASE_STUDIES = True
CASE_MAX_QUBITS = 6
SAVE_CASE_RESULTS = True
```

Return to the repository root and execute the notebook:

```bash
cd ..

jupyter nbconvert \
  --to notebook \
  --execute \
  --inplace \
  qhat_L_sweep_trotter_demo.ipynb \
  --ExecutePreprocessor.timeout=3600
```

This creates or updates:

```text
qhat_full_library_coloring_trotter_errors.csv
```

## 9. Understand the CSV size

The reference CSV contains 1,728 rows:

\[
108\ \text{cases}
\times2\ \text{coloring schemes}
\times2\ \text{Trotter orders}
\times4\ \text{step counts}
=1728.
\]

The 108 exact-evolution cases consist of all successful JW cases at two, four, and six qubits. Each row records:

- molecule;
- bond length;
- basis;
- active occupied and vacant orbitals;
- qubit count;
- coloring scheme;
- graph vertices, edges, and colors;
- Trotter order and step count;
- spectral error; and
- tensor-to-JW reconstruction error.

## 10. Run the focused three-case robustness benchmark

After generating the L-sweep library, run the representative benchmark from the
repository root:

```bash
python run_trotter_benchmark.py --random-orderings 100
python summarize_trotter_results.py
```

The default selection contains one four-qubit near-tie, one six-qubit case where
fermionic coloring is better at first order, and one six-qubit case where direct
JW coloring is much better at second order. The runner evaluates first and
second order at 1, 2, 5, and 10 steps. It saves the graph/color schedule,
color-block and nested-commutator diagnostics, 100 seeded random vertex orders,
state error, infidelity, and particle-number leakage.

The portable outputs are written to `trotter_benchmark_results/`:

- `ordering_trials.csv`: colored and random-order measurements;
- `graph_diagnostics.csv`: vertices, edges, colors, color membership/order, and
  coefficient-weighted commutator summaries;
- `block_commutators.csv`: pairwise color-block diagnostics;
- `robustness_summary.csv`: best, median, percentile, and relative-error metrics;
- `benchmark_summary.md`: the compact numeric report; and
- `benchmark_config.json`: exact cases, seed, time, and step settings.

The notebook's **Focused benchmark** section reads these outputs and plots the
colored schedule against the random median and 10th–90th percentile band.

The generated `demo_L_sweep_library/` directory is intentionally ignored by
Git. Consequently, a clean clone must generate the L-sweep inputs before
rerunning the benchmark. The regression tests run the data-independent checks
everywhere and automatically add the three numeric case checks when those local
inputs are present.

## 11. Disk usage and Git guidance

The complete generated library can require several gigabytes because of the Hartree–Fock HDF5 and pickle intermediates. In the reference run, the molecular library occupied approximately 6.3 GB, while the final CSV, notebook, and Markdown reports were small enough to review and publish normally.

Do not commit the full generated library unless the repository explicitly uses an appropriate large-file storage policy. The portable analysis artifacts are:

- `qhat_L_sweep_trotter_demo.ipynb`;
- `qhat_full_library_coloring_trotter_errors.csv`;
- `QHAT_FERMIONIC_VS_PAULI_NONCOMMUTATION_INSIGHTS.md`; and
- `trotter_benchmark_results/*.csv` and `benchmark_summary.md`;
- `run_trotter_benchmark.py` and `summarize_trotter_results.py`; and
- this reproduction guide.

## Troubleshooting

### Missing dependencies

If `basis_set_exchange`, `mendeleev`, `pyscf`, or `openfermionpyscf` cannot be imported, run the sweep from the environment in which QHAT's Hamiltonian-generation dependencies are installed.

### OpenFermion writes into a read-only installation directory

Some OpenFermion environments choose their installed package data directory when `MolecularData` has no explicit filename. If generation fails with a permission error inside `site-packages/openfermion`, provide a writable filename when constructing `MolecularData` in `hamgen.py`, such as the configured QHAT file stub:

```python
molecule_data = MolecularData(
    geometry,
    state.config_hamiltonian.basis,
    1 + electron_count % 2,
    0,
    filename=state.config_general.file_stub,
)
```

This keeps the intermediate HDF5 file beside the configured QHAT outputs.

### He–He/STO-6G active-space failure

The closed-shell He–He/STO-6G case can have no vacant orbital for the requested active space. Treat those configurations as skipped rather than as missing successful results.

### Dense exact evolution is slow

Matrix dimensions grow as \(2^n\). Keep `CASE_MAX_QUBITS` conservative for spectral-norm calculations, or replace dense evolution with sparse or matrix-free methods for larger active spaces.
