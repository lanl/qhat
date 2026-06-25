# QHAT Standalone: Self-Contained Trotterization Tools

This directory contains a self-contained version of the quantum simulation code with no external dependencies (outside of Julia packages).

## Contents

- **`src/`** - Core simulation modules (copied from parent repo)
- **`Li-Li_jw/`** - Test Hamiltonians for Li-Li molecule (4-16 qubits, Jordan-Wigner encoding)
- **`Li-Li_bond/`** - Li-Li Hamiltonians across multiple bond lengths, bases, active spaces, and encodings
- **`trotter_expts.jl`** - Trotter error demonstration script
- **`trotter_ground_energy.jl`** - Compute exact metadata energy and Trotter-derived ground energies for one Hamiltonian
- **`bond_energy_curve.jl`** - Exact plus second-order Trotter energy-vs-bond-length curves
- **`first_order_bond_energy_curve.jl`** - Exact plus first-order Trotter energy-vs-bond-length curves
- **`Project.toml`** - Julia package dependencies

## Quick Start
Before the first run. Start julia REPL, hit `]` to get into package mode and run `activate .` and `instantiate` to have all the dependencies installed. Then outside the REP do the following: 

```bash
cd QHAT_standalone
julia --project=. trotter_expts.jl
```

This runs the demonstration on the smallest system (4 qubits). To use a different Hamiltonian:

```bash
julia --project=. trotter_expts.jl Li-Li_jw/Li-Li_2.90_hgbs-5_as-004-004_jw.dat  # 8 qubits
```

## Additional Scripts

### `trotter_ground_energy.jl`

Reads a single Hamiltonian, builds second-order Trotter product unitaries, and reports:
- metadata ground-state energy
- Trotter-derived ground-state energies for `nsteps = 1, 5, 10`

By default it uses a partial eigensolve on `U + U'`. You can switch to Arnoldi with the Hartree-Fock start vector:

```bash
julia --project=. trotter_ground_energy.jl
julia --project=. trotter_ground_energy.jl --arnoldi
julia --project=. trotter_ground_energy.jl --arnoldi Li-Li_jw/Li-Li_2.90_hgbs-5_as-004-004_jw.dat
```

### `bond_energy_curve.jl`

Builds a Li-Li potential energy curve across all bond lengths in `Li-Li_bond/` for a fixed:
- basis
- active occupied orbitals
- active vacant orbitals
- encoding

It prints rows as each bond length finishes:

```text
bond_length,exact_energy_hartree,trotter_1_step_hartree,trotter_5_step_hartree,exact_seconds,trotter_1_seconds,trotter_5_seconds,total_seconds,file
```

and saves a plot overlaying:
- exact ground-state energy
- second-order Trotter, 1 step
- second-order Trotter, 5 steps

Examples:

```bash
julia --project=. bond_energy_curve.jl
julia --project=. bond_energy_curve.jl --basis hgbs-5 --electrons 4 --vacant 8 --encoding jw
julia --project=. bond_energy_curve.jl --basis sto-6g --electrons 4 --vacant 8 --encoding jw
```

The default curve is `hgbs-5`, `as-004-004`, `jw`.

### `first_order_bond_energy_curve.jl`

Same interface as `bond_energy_curve.jl`, but the Trotter overlay uses first-order Trotter instead of second-order Trotter.

Examples:

```bash
julia --project=. first_order_bond_energy_curve.jl
julia --project=. first_order_bond_energy_curve.jl --basis hgbs-5 --electrons 4 --vacant 8 --encoding jw
julia --project=. first_order_bond_energy_curve.jl --basis sto-6g --electrons 4 --vacant 8 --encoding jw
```

## Bond-Length Data Layout

The bond scan data is organized as:

```text
Li-Li_bond/<bond_length>/<basis>/Li-Li_<bond>_<basis>_as-EEE-VVV_<encoding>.dat
```

where:
- `EEE` = active occupied single-occupancy orbitals
- `VVV` = active vacant single-occupancy orbitals
- `encoding` = `jw` or `bk`

## What the Demonstration Shows

The `trotter_expts.jl` script provides a comprehensive analysis of Trotterization errors:

### 1. System Information
- Number of qubits and Hamiltonian terms
- Ground state energy
- L1 norm and normalization factor

### 2. State Vector Errors
Compares Trotter evolution against Chebyshev reference:
- **1st order Trotter**: Basic splitting, O(dt) error per step
- **2nd order Trotter**: Strang splitting, O(dt²) error per step
- **4th order Trotter**: Recursive composition, O(dt⁴) error per step

### 3. Operator Errors
Compares full unitary matrices (spectral norm):
- `‖U_trotter - U_exact‖` for varying step counts
- Demonstrates proper O(dt³) scaling for 2nd order methods

### 4. Commutator Error Bounds
Theoretical upper bounds from Childs et al. (PRX 2021):
- Exact nested commutator computation
- Comparison of bounds vs actual errors
- Bound tightness analysis

### 5. Error Scaling Analysis
- Verifies theoretical scaling laws (e.g., 4× error reduction when doubling steps)
- Provides recommendations for achieving target accuracy

## Test Hamiltonians

Seven Li-Li molecule Hamiltonians with varying sizes generate by QHAT

The `Li-Li_bond/` directory contains bond-length scans for:
- bases `hgbs-5` and `sto-6g`
- active spaces `002-002`, `002-004`, `004-004`, `004-006`, `004-008`, `006-008`, `006-010`
- encodings `jw` and `bk`

## Implementation Details

- **Chebyshev order**: 10 (adjustable via `CHEBY_ORDER` constant)
- **Evolution time**: 1.0 (total time, split across Trotter steps)
- **Test step counts**: [10, 20, 50, 100]
- **Total runtime**: ~10-30 seconds for 4-qubit system

For the curve scripts:
- exact energies come from metadata if present, otherwise from sparse eigensolves of the Hamiltonian
- Trotter energies are extracted from the largest eigenvalues of `U + U'`
- the Trotter curve code uses matrix-free Arnoldi actions for `U + U'`, which is important for larger active spaces

