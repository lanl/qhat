# SOSSA Resource Estimates

This repository contains prototype Python code for estimating quantum resources for spectrum-amplified sum-of-squares electronic-structure simulation, abbreviated here as **SOSSA**.

The current implementation is centered on a single Python class, `ResourceEstimates`, in the `SOSSA_resource_estimates.py` module. The class now handles both parts of the workflow:

1. Build and solve the spin-free or spin-explicit SOS semidefinite program.
2. Factor the SDP Gram matrices into SOS generator coefficients.
3. Estimate the sparse support size `d_sparse` from the nonzero entries of the lowercase generator coefficients `g`, `gbar`, `d`, and `q`.
4. Configure and compute sparse-qubitization-style Toffoli and logical-qubit resource estimates.

The implementation is intended for research prototyping and model comparison, not as a finalized fault-tolerant resource estimator.

---

## Repository contents

```text
.
├── SOSSA_resource_estimates.py
├── test_run_SOSSA.ipynb
└── README.md
```

### `SOSSA_resource_estimates.py`

Contains the single main class:

```python
ResourceEstimates
```

`ResourceEstimates` combines the functionality that was previously split across separate SDP and resource-counting classes.

It builds and solves the SOS semidefinite program based on the spin-free / spin-explicit algebra in Appendix G of the SOSSA paper.

The SDP variables are uppercase Gram/PSD matrices:

```text
S_G = [[G, Gp],
       [Gpp, Gppp]]
```

and either

```text
D_up, D_dn, Q_up, Q_dn
```

for spin-explicit calculations, or

```text
D, Q
```

for spin-free calculations.

After solving, the code factors the PSD matrices to obtain lowercase SOS generator coefficients:

```text
g, gbar, d, q
```

These are the coefficients appearing directly in the SOS generators `O_alpha`. The code estimates

```text
d_sparse = nnz(g) + nnz(gbar) + nnz(d) + nnz(q)
```

after applying a user-specified threshold.

The same class then converts the SOSSA quantities into rough Toffoli and logical-qubit estimates using a sparse-QROM/QROAM-inspired block-encoding model.

Resource estimates can be computed in two ways:

1. **Preferred sparse mode**

   After `solve()` succeeds, call `estimate_resources(...)`. By default, this uses `Lambda`, `E_SOS`, total SOS rank `R`, and `d_sparse` from the most recent successful solve.

2. **Dense fallback mode**

   If `d_sparse` is not supplied or available, the table size is estimated from a dense formula depending on the number of spin orbitals and the SOS rank.

### `test_run_SOSSA.ipynb`

Example notebook for building molecular integrals, running the SOS SDP, extracting `d_sparse`, and computing SOSSA resource estimates with the unified `ResourceEstimates` class.

---

## Installation

Create a fresh Python environment:

```bash
conda create -n sossa python=3.10
conda activate sossa
```

Install core dependencies:

```bash
pip install numpy scipy cvxpy pyscf jupyter matplotlib
```

Optional but recommended:

```bash
pip install mosek
```

`MOSEK` generally performs better for semidefinite programs, but it requires a valid MOSEK license. The code can also use `SCS`, which is open source but may be slower or less accurate for larger SDPs.

---

## Minimal usage

```python
from pyscf import gto, scf, ao2mo

from SOSSA_resource_estimates import ResourceEstimates
```

Build a molecule and molecular-orbital integrals:

```python
mol = gto.M(
    atom="H 0 0 0; H 0 0 0.74",
    unit="Angstrom",
    basis="sto-3g",
    charge=0,
    spin=0,
    verbose=0,
)

nao = mol.nao_nr()

mf = scf.RHF(mol)
mf.conv_tol = 1e-10
E_hf = mf.kernel()

C = mf.mo_coeff

hcore_ao = mol.intor("int1e_kin") + mol.intor("int1e_nuc")
h1_mo = C.T @ hcore_ao @ C

eri_mo = ao2mo.kernel(mol, C)
h2_mo = ao2mo.restore(1, eri_mo, nao)
```

Create the unified SOSSA estimator and solve the SOS SDP:

```python
solver = ResourceEstimates(
    norb=nao,
    h1=h1_mo,
    h2=h2_mo,
    spin_explicit=False,
    enforce_spin_free=True,
    use_triu_pairs=True,
    h2_abs_threshold=1e-6,
    support_threshold=1e-8,
    solver="SCS",
    max_iterations=200000,
    verbose=False,
)

out = solver.solve()
```

Inspect the SDP result:

```python
print(out["diagnostics"])
print(out["ESOS"])
print(out["Lambda"])
print(out["total_rank"])
print(out["sos_support"])
```

Extract the sparse support if desired:

```python
d_sparse = out["sos_support"]["d_sparse"]
```

Run the SOSSA resource estimator from the same object:

```python
epsilon = 1.0e-3
M_spin_orbitals = 2 * nao

resource_summary = solver.estimate_resources(
    E_0=E_hf,
    eps_qpe=epsilon,
    M_spin_orbitals=M_spin_orbitals,
    bits_precision=20,
    optimum=True,
    reduced_diagonal_doublets=True,
    real_coefficients=True,
    include_a16_overheads=True,
    include_qpe_qubits=True,
)

print(resource_summary)
```

`estimate_resources(...)` automatically uses `Lambda`, `E_SOS`, `R`, and `d_sparse` from the most recent successful `solve()` call unless you explicitly override them.

You can also configure the resource-estimation parameters first and call `summary()` later:

```python
solver.configure_resource_estimates(
    E_0=E_hf,
    eps_qpe=epsilon,
    M_spin_orbitals=M_spin_orbitals,
    bits_precision=20,
)

print(solver.summary())
```

---

## Example output fields

The SDP solver returns a dictionary with fields like:

```python
out["ESOS"]
out["lambda_sqrt2"]
out["Lambda"]
out["coeffs"]
out["diagnostics"]
out["total_rank"]
out["sos_support"]
```

The support dictionary contains information such as:

```python
out["sos_support"]["support_model"]
out["sos_support"]["threshold"]
out["sos_support"]["d_sparse"]
out["sos_support"]["nnz_g"]
out["sos_support"]["nnz_gbar"]
out["sos_support"]["nnz_d"]
out["sos_support"]["nnz_q"]
```

The resource-estimate summary contains fields like:

```python
resource_summary["lambda_eff"]
resource_summary["num_walk_steps"]
resource_summary["block_encoding_toffoli"]
resource_summary["walk_toffoli"]
resource_summary["total_toffoli"]
resource_summary["logical_qubits"]
```

The class also exposes helper methods such as:

```python
solver.support_summary()
solver.lambda_eff()
solver.num_walk_steps()
solver.block_encoding_toffoli()
solver.walk_toffoli()
solver.toffoli_count()
solver.qubit_count()
```

These methods require resource-estimation parameters to be configured first, either through `estimate_resources(...)` or `configure_resource_estimates(...)`.

---

## Sparse support model

The current sparse-support model counts nonzero entries in the factorized SOS generator coefficients.

After the SDP is solved, the uppercase PSD matrices are factored as:

```text
S_G = Y_SF.T @ Y_SF
D   = Y_D.T  @ Y_D
Q   = Y_Q.T  @ Y_Q
```

The rows of the factor matrices define the lowercase generator coefficients:

```text
g, gbar, d, q
```

Then

```text
d_sparse = nnz(g) + nnz(gbar) + nnz(d) + nnz(q)
```

where

```text
nnz(x) = number of entries satisfying abs(x) > support_threshold
```

For spin-free runs, the code stores the spin-free `d` and `q` factors in the `d_up` and `q_up` fields, while `d_dn` and `q_dn` are empty placeholders.

For spin-explicit runs, the support is counted across both spin sectors:

```text
d_sparse = nnz(g) + nnz(gbar)
         + nnz(d_up) + nnz(d_dn)
         + nnz(q_up) + nnz(q_dn)
```

---

## Important conventions

### Spatial orbitals vs spin orbitals

`ResourceEstimates` expects `norb` to be the number of **spatial orbitals** when constructing the SDP problem.

```python
solver = ResourceEstimates(norb=nao, h1=h1_mo, h2=h2_mo)
```

The resource-estimation part expects `M_spin_orbitals` to be the number of **spin orbitals**.

```python
M_spin_orbitals = 2 * nao
```

If `M_spin_orbitals` is omitted from `configure_resource_estimates(...)` or `estimate_resources(...)`, the class defaults to `2 * norb`.

### Energy quantities

`E_SOS` is the SOS lower-bound energy returned by the SDP.

`E_0` is typically set to a reference upper-bound energy, such as the Hartree-Fock energy or another variational estimate.

The effective spectrum-amplified normalization is modeled as:

```text
lambda_eff = sqrt(E_gap * abs(2 Lambda - E_gap))
```

where the implementation defines

```text
E_gap = abs(abs(E_0) - E_SOS)
```

### Solver choices

Use `solver="SCS"` for a fully open-source run.

Use `solver="MOSEK"` if MOSEK is installed and licensed.

---

## Running the notebook

Start Jupyter:

```bash
jupyter notebook
```

Open:

```text
test_run_SOSSA.ipynb
```

Run the cells in order.

The notebook demonstrates:

1. Importing `ResourceEstimates` from `SOSSA_resource_estimates`.
2. Building molecular integrals with PySCF.
3. Solving the SOS SDP with `solver.solve()`.
4. Extracting `d_sparse` from `out["sos_support"]`.
5. Running resource estimates with `solver.estimate_resources(...)`.
6. Printing Toffoli and logical-qubit counts.

---

## Notes and caveats

This code is experimental.

The resource model is a research-level estimate, not a complete compiled circuit cost. The Toffoli and qubit estimates are intended for comparing scaling and relative costs under the assumptions encoded in the model.

The sparse support `d_sparse` is based on thresholded nonzero entries of the factorized SOS generator coefficients. Changing `support_threshold` may significantly change the resource estimates.

The SDP solve can be numerically sensitive. Check:

```python
out["diagnostics"]["status"]
out["diagnostics"]["residual_one_e_maxabs"]
out["diagnostics"]["g10_rows_dropped"]
```

before trusting resource estimates.

For larger systems, SDP memory and runtime can grow quickly. MOSEK is recommended when available.

---

## Citation

This repository implements prototype ideas inspired by spectrum-amplified sum-of-squares electronic-structure simulation and the spin-free SOS formalism described in:

```text
Fast quantum simulation of electronic structure by spectrum amplification
arXiv:2502.15882
```
