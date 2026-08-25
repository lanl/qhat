# Trotter Error Baseline

## Purpose

Reference measurements of Trotterization error for the QHAT molecular Hamiltonian library.

Results cover all 9 molecules and 10 atoms in both basis sets (sto-6g and hgbs-5) and both fermion-to-qubit mappings (JW and BK).

For each Hamiltonian, **state error** is the distance between Trotterized and exact evolution under that Hamiltonian, starting from the Hartree-Fock (HF) initial state:

```
err = || psi_trotter(T) - psi_exact(T) ||_2
```

Both states start from the same HF determinant and evolve under the same generator, so no global phase enters and the number is pure Trotter error.

```
compute_trotter_state_error.py                             #compute the baseline
plot_state_errors.py                                       #order x delta_t grid; one figure per (molecule, T)
plot_state_errors_by_molecule.py                           #one panel per molecule, grouped by period
plot_state_errors_vs_T.py                                  #x = T; one figure per (molecule, active space)
state_errors/<basis>/<mapping>/state_errors_<species>.csv  #File structure
log_<molecule>.txt                                         #per-molecule run logs
```

## Config

**Systems.** Nine diatomics and ten single atoms:

```
molecules   H-H  He-H  He-He  Li-H  Li-Li  Be-H  Be-Be  B-H  B-B
atoms       H  He  Li  Be  B  C  N  O  F  Ne
```

**Bond lengths.** Ten points per diatomic, evenly spaced over `L in [0.6, 3.0] * L_eq`, with `L_eq` the sum of the Pyykko covalent radii. The range therefore runs from compressed, through equilibrium, out to dissociated. Atoms have no bond length and carry `L = nan`.

**Active spaces.** Grown from the smallest allowed up to 12 qubits.

**Generator.** `H_terms = pi * (H - shift*I) / normalization`, identity excluded, with `shift` the identity coefficient and `normalization = 2 * max(1, one_norm - |shift|)`.

**Initial state.** Lowest `n_occ` spin orbitals occupied, from the header key `number of active, occupied, single-occupancy orbitals`. JW stores that occupation directly; BK stores `(encoder_bk @ occ) mod 2`, as in `hamgen.get_initial_state`.

**Term ordering.** Ascending signed coefficient, ties broken by ascending lexicographic Pauli string, so the sequence is deterministic across runs and machines.

**Reference.** `scipy.sparse.linalg.expm_multiply` on the same generator.

**Grid.** `n_steps` is per U, where `U = exp(-i H_terms)` is one unit of evolution, so `delta_t = 1/n_steps` and `total_steps = T*n_steps`. `T` is swept and identical for every Hamiltonian, not a per-molecule lookup. Cells are total Trotter steps:

| T | Δt = 1 | Δt = 1/2 | Δt = 1/4 |
|---|---|---|---|
| 2⁶ = 64 | 64 | 128 | 256 |
| 2⁸ = 256 | 256 | 512 | 1024 |
| 2⁹ = 512 | 512 | 1024 | 2048 |

Orders 1, 2, 4 → **27 rows per (Hamiltonian, mapping)**. Both bases (hgbs-5, sto-6g) and both mappings (JW, BK) run by default.

Every CSV carries a `#` provenance header recording the conditions that produced its rows. Read the files with `pandas.read_csv(path, comment='#')`.

## How to run

### 1. Generate the Hamiltonians

Omitting the positional argument covers all 19 systems — the 9 diatomics plus atoms H through Ne:

```bash
python hamiltonian_generator/build_config_L_sweep.py \
    --L-steps 10 \
    --max-active 12 \
    --library library \
    --hamgen-dir hamiltonian_generator \
    --run
```

Equivalent explicit form, to restrict to a subset:

```bash
python hamiltonian_generator/build_config_L_sweep.py \
    H-H,He-H,He-He,Li-H,Li-Li,Be-H,Be-Be,B-H,B-B,H,He,Li,Be,B,C,N,O,F,Ne \
    --L-steps 10 --max-active 12 --library library \
    --hamgen-dir hamiltonian_generator --run
```

Diatomics land in `library/<mol>/<L>/<basis>/`, atoms in `library/atoms/<atom>/<basis>/`.

### 2. Compute the baseline
Run from the baseline directory.
```bash
# everything
python compute_trotter_state_error.py ../library/ --out-dir state_errors --n-L 10

# one species (atoms too: --molecule Ne)
python compute_trotter_state_error.py ../library/ --molecule Li-Li --out-dir state_errors

# the 9 diatomics in parallel; species write disjoint files, so one tree is safe
for m in H-H He-H He-He Li-H Li-Li Be-H Be-Be B-H B-B; do
    python compute_trotter_state_error.py ../library/ --molecule "$m" \
        --out-dir state_errors --n-L 10 > "log_$m.txt" 2>&1 &
done
wait

# the 10 atoms, same way. Point at library/atoms so the scan is quicker, and
# drop --n-L: atoms have no bond length, so there is nothing to subsample.
for a in H He Li Be B C N O F Ne; do
    python compute_trotter_state_error.py ../library/atoms --molecule "$a" \
        --out-dir state_errors > "log_$a.txt" 2>&1 &
done
wait
```

Resumable: existing rows are skipped, each row is flushed as computed, so an interrupted run loses nothing.

Groups are processed in sorted order — `hgbs-5/bk`, `hgbs-5/jw`, `sto-6g/bk`, `sto-6g/jw` — so missing directories mid-run are expected. He-He has no sto-6g Hamiltonians (no vacant orbitals in a minimal basis) and finishes with two CSVs.

### 3. Plot

Generates three different types of plots using the saved state error data. These are the commands that produced the committed `plots/` folder:

```bash
# active space on x, at fixed T; one figure per (molecule, T)
python plot_state_errors.py state_errors/hgbs-5/jw \
    --out-dir plots/active_space_vs_state_error_fixed_T_hgbs_jw
python plot_state_errors.py state_errors/hgbs-5/bk \
    --out-dir plots/active_space_vs_state_error_fixed_T_hgbs_bk

# one panel per molecule, grouped by period; colour = L / L_eq
python plot_state_errors_by_molecule.py state_errors/hgbs-5/jw \
    --order 1,2,4 --n-steps 1,2,4 \
    --out-dir plots/active_space_vs_state_error_all_molecule_hgbs_jw
python plot_state_errors_by_molecule.py state_errors/hgbs-5/bk \
    --order 1,2,4 --n-steps 1,2,4 \
    --out-dir plots/active_space_vs_state_error_all_molecule_hgbs_bk

# T on x, at fixed active space; one figure per (molecule, active space)
python plot_state_errors_vs_T.py state_errors/hgbs-5/jw \
    --out-dir plots/T_vs_state_error_fixed_as_hgbs_jw
python plot_state_errors_vs_T.py state_errors/hgbs-5/bk \
    --out-dir plots/T_vs_state_error_fixed_as_hgbs_bk
```


Narrower slices, and a single molecule:

```bash
python plot_state_errors.py state_errors/hgbs-5/jw/state_errors_Li-Li.csv --T 512
python plot_state_errors_vs_T.py state_errors/hgbs-5/jw --qubits 12 --out-dir figs
```

Common flags: `--T`, `--orders`, `--n-steps`, `--out-dir`, `--cmap`.

## Evaluating a new technique
 
To measure a new Trotter or time-evolution technique against this baseline, modify `compute_trotter_state_error.py` to use the new technique and write a CSV in the same format, then modify a plot script to show both series.
 
**Keep everything else fixed.** The comparison only means something if the generator, normalization, initial state, reference and grid are unchanged — a lower `err` under a different `delta_t` or a different reference is not an improvement, it is a different measurement. The provenance header at the top of each CSV records those conditions; the new run's header should match the baseline's except for the technique itself.
 
**Keep the CSV schema.** Same columns in the same order, so the two files can be concatenated and grouped. The plot scripts key on `molecule`, `L`, `basis`, `active`, `mapping`, `n_qubits`, `T`, `n_steps`, `order` and `err`.

 **Plot both.** Add a `source` column ("baseline" / "new") to each frame, concatenate them, and draw one line style per source in whichever plot script suits the comparison — `plot_state_errors.py` and `plot_state_errors_vs_T.py` for a single molecule across active spaces, `plot_state_errors_by_molecule.py` for coverage across the library.

 ## Future work
 
- Extend to larger active spaces, beyond the current 12-qubit cap.
- Add an energy error baseline: error in the eigenphase QPE recovers, not in the state vector.
- Add a Hamiltonian error baseline: operator distance between the Trotterized and exact unitaries, with commutator bounds.