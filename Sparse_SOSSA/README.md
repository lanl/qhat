# SOSSA Resource Estimates

This repository contains prototype Python code for estimating quantum resources for spectrum-amplified sum-of-squares electronic-structure simulation, abbreviated here as **SOSSA**.

The workflow is:

1. Build molecular one- and two-electron integrals.
2. Solve a spin-free or spin-explicit sum-of-squares semidefinite program.
3. Factor the SDP Gram matrices into SOS generator coefficients.
4. Estimate a sparse support size `d_sparse` from the nonzero entries of the lowercase generator coefficients `g`, `gbar`, `d`, and `q`.
5. Use `d_sparse` in a sparse-qubitization-style resource model for Toffoli and logical-qubit estimates.

The implementation is intended for research prototyping and model comparison, not as a finalized fault-tolerant resource estimator.

---

## Repository contents

```text
.
├── one_norm_SOSSA_SDP.py
├── SOSSA_resource_count.py
├── test_run_SOSSA.ipynb
└── README.md
```

### `one_norm_SOSSA_SDP.py`

Contains the `ChemSOSSDP` class.

This class builds and solves the SOS semidefinite program based on the spin-free / spin-explicit algebra in Appendix G of the SOSSA paper.

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

### `SOSSA_resource_count.py`

Contains the `SOSSAResourceEstimates` class.

This class converts SOSSA quantities into rough Toffoli and logical-qubit estimates using a sparse-QROM/QROAM-inspired block-encoding model.

It supports two modes:

1. **Dense fallback mode**

   If `d_sparse=None`, the table size is estimated from a dense formula depending on the number of spin orbitals and the SOS rank.

2. **Generator-coefficient sparse mode**

   If `d_sparse` is provided, the class uses it directly as the sparse table size.

This second mode is the preferred mode for the current code.

### `test_run_SOSSA.ipynb`

Example notebook for building molecular integrals, running the SOS SDP, extracting `d_sparse`, and computing SOSSA resource estimates.

---

## Installation

Create a fresh Python environment:

```bash
conda create -n sossa python=3.10
conda activate sossa
```

Install core dependencies:

```bash
pip install numpy scipy cvxpy pyscf jupyter
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

from one_norm_SOSSA_SDP import ChemSOSSDP
from SOSSA_resource_count import SOSSAResourceEstimates
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

Solve the SOS SDP:

```python
solver = ChemSOSSDP(
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

Extract the sparse support:

```python
d_sparse = out["sos_support"]["d_sparse"]
```

Run the SOSSA resource estimator:

```python
epsilon = 1.0e-3
M_spin_orbitals = 2 * nao

estimates = SOSSAResourceEstimates(
    Lambda=out["Lambda"],
    E_SOS=out["ESOS"],
    E_0=E_hf,
    eps_qpe=epsilon,
    M=M_spin_orbitals,
    R=out["total_rank"],
    bits_precision=20,
    k_r=None,
    optimum=True,
    d_sparse=d_sparse,
)

print(estimates.summary())
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

The resource estimator summary contains fields like:

```python
estimates.summary()["lambda_eff"]
estimates.summary()["num_walk_steps"]
estimates.summary()["block_encoding_toffoli"]
estimates.summary()["walk_toffoli"]
estimates.summary()["total_toffoli"]
estimates.summary()["logical_qubits"]
```

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

`ChemSOSSDP` expects `norb` to be the number of **spatial orbitals**.

```python
solver = ChemSOSSDP(norb=nao, ...)
```

`SOSSAResourceEstimates` expects `M` to be the number of **spin orbitals**.

```python
M = 2 * nao
```

### Energy quantities

`E_SOS` is the SOS lower-bound energy returned by the SDP.

`E_0` is typically set to a reference upper-bound energy, such as the Hartree-Fock energy or another variational estimate.

The effective spectrum-amplified normalization is modeled as:

```text
lambda_eff = sqrt(|E_0 - E_SOS| * |2 Lambda - |E_0 - E_SOS||)
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

1. Building molecular integrals with PySCF.
2. Solving the SOS SDP.
3. Extracting `d_sparse`.
4. Running the resource estimator.
5. Printing Toffoli and logical-qubit counts.

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
