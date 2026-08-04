# Running the L-sweep Fermionic-Term Trotter Benchmark

This runbook shows how to execute the four-method L-sweep benchmark, resume a
run safely, select a subset of methods, and regenerate plots. Run every command
from the QHAT repository root.

## 1. Prepare the environment

Activate an environment containing QHAT's analysis and Hamiltonian-generation
dependencies, including NumPy, SciPy, NetworkX, pandas, Matplotlib, and
OpenFermion. For the repository's Conda environment:

```bash
conda activate qhat
```

Confirm that the command-line programs load:

```bash
python analysis/benchmark_L_sweep_trotter.py --help
python analysis/visualize_L_sweep_trotter.py --help
python analysis/visualize_all_L_sweep_cases.py --help
```

## 2. Confirm that tensor inputs are available

The benchmark reads QHAT `*.tensors.npz` files containing `constant`,
`one_body`, and `two_body` arrays. The default input root is
`hamiltonian_generator/library`.

```bash
find hamiltonian_generator/library -name '*.tensors.npz' | head
```

If that command prints nothing, generate or restore the L-sweep tensor
artifacts before attempting the full run. The benchmark does not reconstruct
missing molecular tensors from historical CSV files.

A normal path has this layout:

```text
hamiltonian_generator/library/
└── H-H/
    └── 0.70/
        └── sto-6g/
            └── H-H_0.70_sto-6g_as-002-002.tensors.npz
```

## 3. Understand the four methods

The canonical method names accepted by `--orderings` are:

- `jw_raw`: final combined non-identity JW Pauli strings in insertion order.
- `jw_coloring`: the same Pauli factors ordered by greedy coloring of the
  Pauli-string noncommutation graph.
- `fermionic_coloring`: legacy fermionic-induced JW ordering. It colors
  complete Hermitian fermionic terms but still exponentiates individual final
  JW Pauli strings.
- `fermionic_term_coloring`: genuine fermionic-term factorization. It colors
  complete Hermitian fermionic terms and exponentiates each complete mapped JW
  matrix as one dense factor.

`jw_raw` is automatically added when omitted because historical ratio columns
use it as a baseline.

## 4. Run a small smoke test

The following command selects the smallest discovered tensor case, evaluates
all four methods, uses first- and second-order formulas, and uses two small step
counts:

```bash
python analysis/benchmark_L_sweep_trotter.py \
  --library hamiltonian_generator/library \
  --limit 1 \
  --orderings \
    jw_raw \
    jw_coloring \
    fermionic_coloring \
    fermionic_term_coloring \
  --formula-orders 1 2 \
  --steps 1 2 \
  --time 1 \
  --output /tmp/qhat_l_sweep_smoke.csv
```

To include the fourth-order Suzuki formula in the same smoke matrix:

```bash
python analysis/benchmark_L_sweep_trotter.py \
  --library hamiltonian_generator/library \
  --limit 1 \
  --orderings \
    jw_raw \
    jw_coloring \
    fermionic_coloring \
    fermionic_term_coloring \
  --formula-orders 1 2 4 \
  --steps 1 2 \
  --time 1 \
  --output /tmp/qhat_l_sweep_smoke_all_orders.csv
```

Add `--skip-operator-norm` when only state-error metrics are needed. This
avoids the spectral-norm calculation but does not avoid dense matrix
exponentials for `fermionic_term_coloring`.

## 5. Resume without duplicating completed rows

Repeat the exact command with `--resume`:

```bash
python analysis/benchmark_L_sweep_trotter.py \
  --library hamiltonian_generator/library \
  --limit 1 \
  --orderings \
    jw_raw \
    jw_coloring \
    fermionic_coloring \
    fermionic_term_coloring \
  --formula-orders 1 2 4 \
  --steps 1 2 \
  --time 1 \
  --output /tmp/qhat_l_sweep_smoke_all_orders.csv \
  --resume
```

Resume keys include case, method (`ordering`), formula order, step count, and
evolution time. Completed rows are skipped. Known legacy CSV headers are
migrated to the appended schema before new rows are written; incompatible or
reordered headers are rejected.

## 6. Run only the genuine fermionic-term method

Use an explicit subset when developing or profiling the dense-factor path:

```bash
python analysis/benchmark_L_sweep_trotter.py \
  --library hamiltonian_generator/library \
  --limit 1 \
  --orderings fermionic_term_coloring \
  --formula-orders 1 2 \
  --steps 1 2 \
  --time 1 \
  --output /tmp/qhat_l_sweep_fermionic_term_only.csv
```

The output will also contain `jw_raw` rows for the required baseline ratios.
No `jw_coloring` or legacy `fermionic_coloring` rows are added implicitly.

## 7. Filter cases

`--case-pattern` matches a substring anywhere in the tensor path. For example:

```bash
python analysis/benchmark_L_sweep_trotter.py \
  --library hamiltonian_generator/library \
  --case-pattern 'H-H_0.70_sto-6g_as-002-002' \
  --formula-orders 1 2 4 \
  --steps 10 100 \
  --time 1 \
  --output /tmp/qhat_h2_selected_case.csv
```

Combine `--case-pattern` with `--limit` to stage larger families gradually.

## 8. Run the historical full matrix

The historical L-sweep result set contains 153 cases, evolution times 1/2/5,
formula orders 1/2/4, and step counts 10/100/1000. With four methods, a
complete rerun contains:

```text
153 cases * 3 times * 3 formula orders * 3 step counts * 4 methods
= 16,524 successful rows
```

Use separate resume-safe files for each evolution time:

```bash
for qhat_time in 1 2 5; do
  python analysis/benchmark_L_sweep_trotter.py \
    --library hamiltonian_generator/library \
    --orderings \
      jw_raw \
      jw_coloring \
      fermionic_coloring \
      fermionic_term_coloring \
    --formula-orders 1 2 4 \
    --steps 10 100 1000 \
    --time "$qhat_time" \
    --output "analysis/l_sweep_trotter_state_t${qhat_time}_v2.csv" \
    --resume
done
```

Do not reuse the historical three-method filenames unless an in-place schema
migration and mixed old/new output are intentional.

## 9. Generate per-case and summary plots

Plot one case and create the ranking, win-rate, and graph-complexity summaries:

```bash
python analysis/visualize_L_sweep_trotter.py \
  --csv /tmp/qhat_l_sweep_smoke_all_orders.csv \
  --case-id H-H_0.70_sto-6g_as-002-002 \
  --metric state_infidelity \
  --output-dir /tmp/qhat_l_sweep_smoke_figures
```

Plot the full three-time result matrix:

```bash
python analysis/visualize_L_sweep_trotter.py \
  --csv \
    analysis/l_sweep_trotter_state_t1_v2.csv \
    analysis/l_sweep_trotter_state_t2_v2.csv \
    analysis/l_sweep_trotter_state_t5_v2.csv \
  --metric state_infidelity \
  --output-dir analysis/l_sweep_figures_v2
```

Create one figure for every discovered case and write a manifest:

```bash
python analysis/visualize_all_L_sweep_cases.py \
  --csv \
    analysis/l_sweep_trotter_state_t1_v2.csv \
    analysis/l_sweep_trotter_state_t2_v2.csv \
    analysis/l_sweep_trotter_state_t5_v2.csv \
  --expected-cases 153 \
  --metric state_infidelity \
  --output-dir analysis/l_sweep_all_case_figures_v2
```

Both visualization programs also accept historical three-method CSV files and
method subsets.

## 10. Inspect the output schema

Every new row has `schema_version=2`. The fields that make the factorization
semantics explicit are:

```text
ordering
factorization_level
ordering_method
mapping
number_of_factors
factor_exponential_type
factor_reconstruction_error
graph_level
schema_version
```

For a successful `fermionic_term_coloring` row, expect:

```text
factorization_level=fermionic_term
ordering_method=greedy_coloring
mapping=jordan_wigner
factor_exponential_type=dense_hermitian_expm
graph_level=fermionic_term
```

`factor_reconstruction_error` should be close to floating-point roundoff. A
nonzero identity coefficient is reported separately but omitted from the exact
and product-formula Hamiltonians because it only contributes a global phase.

## 11. Run the relevant tests

From the repository root, make the parent directory importable as `qhat`:

```bash
PYTHONPATH=.. pytest -q analysis/tests/test_L_sweep_fermionic_terms.py
```

Run the full local analysis/common suite:

```bash
PYTHONPATH=.. pytest -q analysis/tests common/tests
```

If Numba cannot write its cache in the active environment, redirect the cache
to a writable temporary directory:

```bash
MPLCONFIGDIR=/tmp/qhat-matplotlib \
NUMBA_CACHE_DIR=/tmp/qhat-numba \
PYTHONPATH=.. \
pytest -q analysis/tests common/tests
```

## 12. Computational scaling

All methods construct dense exact evolution matrices. The analytic Pauli
factor path is comparatively cheap, but each complete fermionic-term factor
uses a general dense matrix exponential. For `n` qubits, one dense factor
requires `O(2^(2n))` matrix storage and approximately `O(2^(3n))` exponential
time. Start with `--limit 1`, then increase the case count while monitoring
memory and runtime. The benchmark imposes no hidden qubit limit and does not
silently skip expensive cases.
