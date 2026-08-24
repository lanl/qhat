# Large-Matrix Low-Energy Diagonalization Workflow

This directory contains a workflow for:

1. generating QHAT/OpenFermion Hamiltonian pickle files from `.config` inputs,
2. computing a small number of the lowest-energy eigenvalues without constructing the full many-body Hamiltonian matrix,
3. running several diagonalizations in batch, and
4. collecting the results into a tab-separated table.

The central diagonalization script is `qhat_low_energies_fci.py`. It converts an OpenFermion-style `InteractionOperator` pickle back to spatial-orbital integrals and uses PySCF's direct FCI solver in a fixed `(N_alpha, N_beta)` electron sector.

> **Important method note:** large sparse Hermitian eigenproblems are commonly approached with the **Lanczos method**, and this README emphasizes that Krylov-subspace viewpoint. The supplied Python code itself currently uses PySCF's **Davidson** solver (`direct_spin1.FCI` with `davidson_only = True`), not a literal Lanczos recurrence. Davidson and Lanczos are related iterative subspace methods, but they are not the same algorithm. A Lanczos algorithm that requires a GPGPU is also in this directory but is of limited use as it requires more that 96GB of VRAM to get past 30 qubit diagonalizations, its name is qhat_pauli_matrix_free.py.

---

## 1. Why use an iterative method for large Hamiltonians?

A Hamiltonian acting on `n` qubits has a formal Hilbert-space dimension of

```text
2^n
```

so explicitly constructing and storing the entire matrix quickly becomes impractical.

For the molecular/Fermionic problem handled here, the calculation can instead be restricted to a fixed electron-number and spin sector. If there are `norb` spatial orbitals with `N_alpha` alpha electrons and `N_beta` beta electrons, the CI-vector dimension is

```text
C(norb, N_alpha) * C(norb, N_beta)
```

The supplied solver works directly with Hamiltonian-vector contractions in this sector. It does **not** construct the full `2^n x 2^n` Hamiltonian matrix.

This is exactly the kind of problem for which iterative eigensolvers such as Lanczos or Davidson are useful: usually only a few extremal eigenvalues, especially the ground state and a few low-lying excited states, are wanted.

---

## 2. Lanczos diagonalization: the key idea

For a large sparse Hermitian matrix `H`, the Lanczos method starts from a normalized vector `v1` and builds a Krylov subspace

```text
K_m(H, v1) = span{v1, H v1, H^2 v1, ..., H^(m-1) v1}.
```

Rather than diagonalizing `H` directly, Lanczos constructs a much smaller tridiagonal matrix representing `H` in this Krylov basis. Its eigenvalues are Ritz approximations to eigenvalues of the original problem. Extremal eigenvalues, including the lowest energies, often converge long before the full spectrum is known.

The ideal Lanczos recurrence is a short three-term relation,

```text
beta_(j+1) v_(j+1) = H v_j - alpha_j v_j - beta_j v_(j-1),
```

which means the main expensive operation is repeatedly applying `H` to a vector.

### Why that matters here

The essential benefit is that a matrix-free implementation only needs an operation equivalent to

```text
vector_out = H @ vector_in
```

without storing all elements of `H`.

That is also the central strategy of the supplied diagonalizer. PySCF contracts the Hamiltonian with CI vectors directly, avoiding construction of the full many-body sparse matrix.

### Davidson versus Lanczos in this repository

The supplied `qhat_low_energies_fci.py` uses the **Davidson method** through PySCF, not Lanczos. Davidson also builds a low-dimensional trial subspace and extracts Ritz eigenpairs, but it expands the subspace using preconditioned residual vectors rather than the strict Lanczos three-term recurrence.

For this workflow, the practical message is the same:

- do not form the enormous full Hamiltonian matrix,
- work in the physically relevant electron/spin sector,
- request only the lowest few roots,
- iteratively apply the Hamiltonian to trial vectors, and
- control convergence with the subspace size, tolerance, and iteration limit.

If a future implementation is changed to an explicit Lanczos solver, the surrounding config -> pickle -> batch -> table workflow can remain essentially the same.

---

## 3. Files used by the workflow

The supplied files are:

```text
qhat_low_energies_fci.py         # main low-energy Davidson/FCI diagonalizer
RunAllHamgens.cmd                # batch generation of Hamiltonian pickles
RunFCIAllDiagonalizationsK3.cmd  # batch of k=3 low-energy diagonalizations
proc3LowestEnergiesFromH.cmd     # builds a TSV table from diagonalization log files
proc3_lowest_energies.py         # extracts 3 energies from one output directory
requirements.txt                 # Python dependencies
```

The uploaded diagonalizer was named `qhat_low_energies_fci(1).py`, while `RunFCIAllDiagonalizationsK3.cmd` calls `qhat_low_energies_fci.py`. Rename the file if necessary:

```bash
mv 'qhat_low_energies_fci(1).py' qhat_low_energies_fci.py
```

### `.cmd` files are Bash scripts here

Despite their `.cmd` extension, the supplied batch files use Bash syntax such as `|&`, `time`, and `shopt`. Run them with Bash on Linux, a Linux cluster, WSL, or another Bash environment. Do not assume they are Windows `cmd.exe` batch files.

---

## 4. Python version and virtual environment

Use Python 3.10 or newer. The diagonalizer uses modern Python type syntax and should be run in a clean virtual environment.

### Create the environment

From the project directory:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Upgrade packaging tools:

```bash
python -m pip install --upgrade pip setuptools wheel
```

Install the dependencies:

```bash
python -m pip install -r requirements.txt
```

Verify the important packages:

```bash
python -c "import numpy, scipy, openfermion, pyscf; print('Python environment OK')"
```

When returning to the project later, reactivate the environment with:

```bash
source .venv/bin/activate
```

Leave it with:

```bash
deactivate
```

### Why OpenFermion is included

The main solver does not explicitly import OpenFermion for its normal solve path, but it unpickles an OpenFermion `InteractionOperator`-like object. Python usually needs the class's defining package available while unpickling. OpenFermion is also required by the optional `--validate-small` cross-check.

### Possible additional `hamgen.py` dependencies

`hamgen.py` itself was not included with the supplied files, so its complete dependency set cannot be verified here. If it imports additional project-specific or chemistry packages, add them to `requirements.txt` or install them into the same virtual environment.

---

## 5. Step 1: go from `.config` files to Hamiltonian `.pickle` files

`RunAllHamgens.cmd` runs commands of the form

```bash
time python hamgen.py ../../F-F_1.00_hgbs-5_as-010-010_JW.config \
    |& tee OutputF-F_1.00_hgbs-5_as-010-010_JW.txt
```

The supplied batch file contains six such cases.

### Expected setup

Run the batch file from a directory where:

- `hamgen.py` is available as `./hamgen.py`, and
- the config files are reachable using the relative paths written in `RunAllHamgens.cmd`, currently `../../<name>.config`.

Because the paths are relative to the **current working directory**, check your location before running:

```bash
pwd
ls hamgen.py
ls ../../F-F_1.00_hgbs-5_as-010-010_JW.config
```

Then run:

```bash
bash RunAllHamgens.cmd
```

Each `hamgen.py` run is timed and its terminal output is copied to an `Output*.txt` log by `tee`.

### Expected pickle names

The following diagonalization batch expects these pickle files:

```text
F-F_1.00_hgbs-5_as-010-010.pickle
F-F_1.00_hgbs-5_as-010-012.pickle
F-F_1.00_hgbs-5_as-012-012.pickle
F-F_1.00_hgbs-5_as-012-014.pickle
F-F_1.00_hgbs-5_as-014-014.pickle
F-F_1.00_hgbs-5_as-014-016.pickle
```

Confirm that Hamiltonian generation actually produced them:

```bash
ls -lh *.pickle
```

The attached files do not include `hamgen.py`, so this README cannot verify exactly how its config parser chooses the output filename or what every config field means. The config-to-pickle step above reflects the actual commands in `RunAllHamgens.cmd` and the pickle names expected by the diagonalization batch.

---

## 6. Step 2: inspect one Hamiltonian before solving

Before starting an expensive run, use the main script's `--inspect-only` option. For example:

```bash
python qhat_low_energies_fci.py \
    F-F_1.00_hgbs-5_as-012-012.pickle \
    -o test_inspect \
    --num-electrons 12 \
    --spin 0 \
    -k 3 \
    --threads 32 \
    --inspect-only
```

The script reports quantities including:

- number of spin orbitals/qubits,
- number of spatial orbitals,
- `N_alpha` and `N_beta`,
- CI-sector dimension,
- memory for one `float64` CI vector,
- requested number of roots,
- integral storage,
- PySCF memory limit, and
- thread count.

This is useful for detecting a calculation that is too large before entering the iterative solve.

---

## 7. Step 3: diagonalize one Hamiltonian

A representative command is:

```bash
python qhat_low_energies_fci.py \
    F-F_1.00_hgbs-5_as-012-012.pickle \
    -o FCI-F-F_1.00_hgbs-5_as-012-012-24q_lowest3 \
    --num-electrons 12 \
    --spin 0 \
    -k 3 \
    --threads 32 \
    --tol 1e-9
```

### Important options

`--num-electrons N`
: Total active electrons in the calculation.

`--spin S`
: `N_alpha - N_beta = 2*S_z`, **not** spin multiplicity. For an even-electron singlet sector, use `0`.

`-k`, `--num-roots K`
: Number of lowest-energy roots to request.

`--threads N`
: Sets the OpenMP/BLAS thread count used by the solver.

`--tol VALUE`
: Davidson convergence tolerance. The supplied batch uses `1e-9`.

`--max-cycle N`
: Maximum number of Davidson cycles. Default: `100`.

`--max-space N`
: Maximum Davidson trial-subspace size before restart. Default: `20`.

`--max-memory-mb N`
: Explicit PySCF memory limit. If omitted or zero, the script attempts to use about 80% of physical memory.

`--save-vectors`
: Also save each returned CI eigenvector. This can require substantial disk space.

`--inspect-only`
: Load and analyze the Hamiltonian and sector dimensions but do not solve.

`--validate-small`
: For small problems, compare the roots against an OpenFermion number-preserving sparse solve. The script restricts this check to at most 18 spin orbitals and currently supports `--spin 0` for the automatic comparison.

---

## 8. Output from one diagonalization

The output directory contains:

```text
eigenvalues.npy       # returned low-energy eigenvalues
lowest_energies.txt   # human-readable root index and energy table
diagonal_values.npy   # same returned low-energy values; NOT the full spectrum
metadata.json         # solver, sector, convergence, timing, and memory metadata
ci_vector_*.npy       # only when --save-vectors is requested
```

The script also writes progress and the final roots to standard output. The batch files pipe this output into `OutFCITrace-*.txt` logs.

A successful final section looks conceptually like:

```text
Lowest energies (Hartree):
     0  <energy>  converged
     1  <energy>  converged
     2  <energy>  converged
Solve time: <seconds> s
Output: <directory>
```

If one or more roots do not meet the convergence tolerance, the script warns and exits with status `2`. Typical responses are to request fewer roots, increase `--max-cycle`, or increase `--max-space` while watching memory use.

---

## 9. Step 4: run all diagonalizations

The supplied `RunFCIAllDiagonalizationsK3.cmd` runs six calculations with `k=3` and 32 threads.

Run it under Bash:

```bash
bash RunFCIAllDiagonalizationsK3.cmd
```

Each line follows the pattern

```bash
time python qhat_low_energies_fci.py INPUT.pickle \
    -o OUTPUT_DIRECTORY \
    --num-electrons N \
    --spin 0 \
    -k 3 \
    --threads 32 \
    --tol 1e-9 \
    |& tee OutFCITrace-CASE-k3.txt
```

### Check the electron counts carefully

The supplied batch uses:

```text
010-010 -> 10 electrons
010-012 -> 10 electrons
012-012 -> 12 electrons
012-014 -> 12 electrons
014-014 -> 14 electrons
014-016 -> 14 electrons
```

These values are part of the physical sector selection, so they should be verified against the intended active-space definition before a long production run.

### Note about output-directory labels

Several output-directory names in the supplied batch contain the repeated text

```text
...as-012-012-...q_lowest3
```

even when the input pickle is another active-space case. This does not change the input Hamiltonian passed to the solver, but it can make the resulting directory names misleading. Consider correcting those `-o` names before running a new production batch.

---

## 10. Step 5: turn completed runs into a table

There are two supplied `proc3*` approaches.

### Recommended: `proc3LowestEnergiesFromH.cmd`

This is the more complete table generator. It reads the `OutFCITrace-*.txt` files and prints these columns:

```text
Filename
Qubits
Orbitals
MaxMemoryMB
LowestEnergy0
LowestEnergy1
LowestEnergy2
SolveTime
```

After all six diagonalizations finish, run:

```bash
bash proc3LowestEnergiesFromH.cmd OutFCITrace-*.txt > lowest_energies_summary.tsv
```

**Do not quote the wildcard** in this particular invocation. Bash should expand `OutFCITrace-*.txt` into the individual filenames before passing them to the script.

Inspect the result:

```bash
cat lowest_energies_summary.tsv
```

If the `column` utility is installed, a convenient terminal view is:

```bash
column -t -s $'\t' lowest_energies_summary.tsv
```

This is the `proc3*` file to use when you want one combined table with run metadata and timings.

### Alternative: `proc3_lowest_energies.py`

This Python script reads one diagonalization output directory, loads its `eigenvalues.npy`, sorts the values, and prints one tab-separated row containing

```text
directory_name    lowest_0    lowest_1    lowest_2
```

For one directory:

```bash
python proc3_lowest_energies.py FCI-F-F_1.00_hgbs-5_as-012-012-24q_lowest3
```

To aggregate several output directories:

```bash
{
    printf 'Directory\tLowestEnergy0\tLowestEnergy1\tLowestEnergy2\n'
    shopt -s nullglob
    for d in FCI-*lowest3; do
        python proc3_lowest_energies.py "$d"
    done
} > lowest_energies_from_npy.tsv
```

This route is useful when the log files are missing but the output directories and `eigenvalues.npy` files are intact. It does **not** include qubit count, orbital count, memory, or solve time.

---

## 11. Suggested end-to-end workflow

Activate the environment:

```bash
source .venv/bin/activate
```

Generate all Hamiltonian pickles:

```bash
bash RunAllHamgens.cmd
```

Confirm the pickles exist:

```bash
ls -lh *.pickle
```

Optionally inspect the largest case first:

```bash
python qhat_low_energies_fci.py \
    F-F_1.00_hgbs-5_as-014-016.pickle \
    -o inspect_014_016 \
    --num-electrons 14 \
    --spin 0 \
    -k 3 \
    --threads 32 \
    --inspect-only
```

Run all diagonalizations:

```bash
bash RunFCIAllDiagonalizationsK3.cmd
```

Check for convergence warnings:

```bash
grep -H -E 'NOT CONVERGED|WARNING:|ERROR:' OutFCITrace-*.txt
```

A grep command that prints nothing is a good sign, but the final energies and metadata should still be checked.

Build the summary table:

```bash
bash proc3LowestEnergiesFromH.cmd OutFCITrace-*.txt > lowest_energies_summary.tsv
```

Display it:

```bash
column -t -s $'\t' lowest_energies_summary.tsv
```

---

## 12. Convergence and performance notes

### Request only the roots you need

The supplied production batch asks for only three roots:

```bash
-k 3
```

For very large CI spaces, asking for many roots increases both time and subspace/vector storage.

### Use `--inspect-only` before very large jobs

The script estimates the CI dimension and memory for a single vector. This is a useful first safety check because a matrix-free algorithm avoids the full Hamiltonian matrix but still needs CI vectors and Davidson work vectors.

### Threads

`--threads 32` sets common BLAS/OpenMP environment variables and asks PySCF to use that thread count. More threads do not always provide linear speedup; use a value appropriate to the node allocation.

### Davidson subspace

If convergence is difficult:

```bash
--max-cycle 200 --max-space 30
```

may help, but larger subspaces require more memory. Change one parameter at a time and keep the logs for comparison.

### Saved eigenvectors

Do not use `--save-vectors` unless the CI vectors are actually needed. A single vector has one `float64` value per determinant in the selected sector, so storage can become large.

---

## 13. Common problems

### `python: can't open file 'qhat_low_energies_fci.py'`

The uploaded file may still be named `qhat_low_energies_fci(1).py`. Rename it:

```bash
mv 'qhat_low_energies_fci(1).py' qhat_low_energies_fci.py
```

### Pickle fails to load because an OpenFermion module/class is missing

Make sure the virtual environment is active and OpenFermion is installed:

```bash
python -m pip install -r requirements.txt
```

### `hamgen.py` cannot find a config file

The config paths in `RunAllHamgens.cmd` are relative to the directory from which the script is launched. Check:

```bash
pwd
ls ../../*.config
```

and edit the batch file paths if the project was moved.

### Spin-block validation error

The diagonalizer checks whether the tensors look like the standard OpenFermion even/odd alpha/beta spin-orbital ordering. Do not bypass this validation with `--skip-spin-block-check` unless the tensor convention is independently known and verified.

### One or more roots are not converged

Look for `NOT CONVERGED` in the log. Consider increasing `--max-cycle`, increasing `--max-space` if memory allows, or reducing `-k`.

### Post-processing table is empty

For the Bash postprocessor, make sure the wildcard expands:

```bash
ls OutFCITrace-*.txt
bash proc3LowestEnergiesFromH.cmd OutFCITrace-*.txt
```

Do not use a quoted literal glob such as:

```bash
# Not recommended here:
bash proc3LowestEnergiesFromH.cmd 'OutFCITrace-*.txt'
```

because the script expects actual file arguments.

---

## 14. Reproducibility recommendations

For production calculations, keep together:

- the exact `.config` files,
- the generated `.pickle` files or a reproducible way to regenerate them,
- `hamgen.py`,
- `qhat_low_energies_fci.py`,
- the batch scripts,
- `requirements.txt`,
- all `Output*.txt` Hamiltonian-generation logs,
- all `OutFCITrace-*.txt` diagonalization logs,
- the FCI output directories and `metadata.json` files, and
- the final `.tsv` summary table.

It is also useful to save the environment snapshot:

```bash
python --version > python-version.txt
python -m pip freeze > requirements-lock.txt
```

`requirements.txt` remains the human-maintained dependency list, while `requirements-lock.txt` records the exact versions used for a particular run.
