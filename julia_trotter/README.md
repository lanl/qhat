# QHAT Standalone: Self-Contained Trotterization Demonstration

This directory contains a self-contained version of the quantum simulation code with no external dependencies (outside of Julia packages).

## Contents

- **`src/`** - Core simulation modules (copied from parent repo)
- **`Li-Li_jw/`** - Test Hamiltonians for Li-Li molecule (4-16 qubits, Jordan-Wigner encoding)
- **`trotter_expts.jl`** - **Main demonstration script** (see below)
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

## Implementation Details

- **Chebyshev order**: 10 (adjustable via `CHEBY_ORDER` constant)
- **Evolution time**: 1.0 (total time, split across Trotter steps)
- **Test step counts**: [10, 20, 50, 100]
- **Total runtime**: ~10-30 seconds for 4-qubit system


