# QHAT Analysis

[[_TOC_]]

This directory contains the modules for performing analysis.  This is run
by
```
python driver.py
```
By default, this will load the setup from `config.py`, but this can be changed by passing the name
of a configuration file as a command line argument
```
python driver.py my_configuration_file.py
```

## Configuration Options

Configuration files are themselves Python scripts, allowing users to use control logic to build up
complex configuration files.

Configuration is broken down by several parts of the processing script.

### General

"General" configuration governs the behavior of the resource analysis script itself, including
logging and output file organization.

#### Output Directory

You can organize all output files in a dedicated directory by setting **`general.output_directory`**.
If set, all generated files (logfile, matrices, eigendecompositions, error analysis, numerical
simulation outputs) will be written to this directory. The directory is created automatically if
it doesn't exist.

**Examples**:
```python
# All outputs go to Be-H/ subdirectory
general.output_directory = "Be-H/"

# All outputs go to nested results/run1/ subdirectory
general.output_directory = "results/run1/"

# No output directory (default) - files written to current directory
general.output_directory = ""
```

**Path behavior**:
- Relative paths in filenames are joined with `output_directory`: 
  - `general.logfile = "logs/debug.log"` → `Be-H/logs/debug.log`
- Absolute paths override `output_directory`:
  - `general.logfile = "/tmp/debug.log"` → `/tmp/debug.log`
- Parent directories are created automatically

#### Logfile

You can configure the log file that the script will write to by setting **`general.logfile`** to
the name of the logfile you want to use.  The default is `analysis.log`. If `output_directory` is
set, the logfile will be written to that directory.

#### Log Level

The log level is set by calling one of the following functions. If you call multiple of these
functions, whichever is called last takes precedence.

- **`general.print_verbose()`** -- Calling this function increases the printouts to be more
  verbose, providing additional information as the resource analysis progresses.  This information
  may be of interest to users, depending on how much detail they want regarding the resource
  analysis.
- **`general.print_debug()`** -- Calling this function increases the printouts beyond even
  "verbose".  The additional information provided by "debug" printing is typically relevant to
  developers more than users.
- **`general.print_default()`** -- Calling this function resets the printout level back to the
  standard verbosity, removing "verbose" and "debug" printouts.

### Loading a Hamiltonian

This script can use Hamiltonians from a variety of sources.  Since each format is slightly
different, your own format may not currently be readable without help from the developers to add an
appropriate reader.  Trying to load more than one Hamiltonian will generate an error.  We currently
support the following inputs source:

- **Second-quantization Hamiltonians**: Second-quantized Hamiltonians can be loaded from HDF5 or
  NumPy files by calling **`hamiltonian.load_second_quantization(filename)`**. The file format is
  automatically detected based on the file extension:
  - `.h5` or `.hdf5`: HDF5 format (currently supports one-body term "1e" and two-body term "2e";
    constant terms and bosons not yet supported)
  - `.npy` or `.npz`: NumPy format with the following fields:
    - "constant" (optional): constant energy term
    - "one_body" (required): one-body tensor
    - "two_body" (required): two-body tensor
    - "bosonic\_scalar" (optional): scalar for bosonic terms
    - "fb\_interaction" (optional): fermion-boson interaction tensor

  If "bosonic\_scalar" and "fb\_interaction" are not present, this is a standard fermionic
  second-quantization Hamiltonian. The scripts in the `hamiltonian_generator` directory create
  Hamiltonians in this format, providing only the one-body and two-body tensors.

  This function accepts optional parameters for specifying fermion-to-qubit and boson-to-qubit
  transformations:
  - `fermion_to_qubit_transform`: Transformation method (default: "JW" for Jordan-Wigner)
  - `boson_to_qubit_transform`: Boson transformation method (default: "binary")
  - `max_bosons_per_state`: Maximum bosons per state (required for bosonic systems, default: None)

After loading a Hamiltonian, the user may have knowledge regarding the bounds on the eigenvalues of
the Hamiltonian, which can be used to optimize certain parts of the resource analysis.  In order to
specify these, use the functions **`hamiltonian.set_energy_lower_bound()`** and
**`hamiltonian.set_energy_upper_bound()`**.  These functions both take an argument specifying the
value of the energy (eigenvalue) bound.  There is an optional second argument specifying whether to
use this value exactly (`exact=True`) or only use it if the automatic bound estimate is a looser
constraint than the bound provided here (`exact=False`).  The last call to each function takes
precedence over prior calls to the same function.

Hamiltonians loaded in second-quantization format can use different fermion-to-qubit and/or
boson-to-qubit transformations. These transformations must be specified as optional parameters
when calling `load_second_quantization()`. The available parameters are:

- **`fermion_to_qubit_transform`**: Specify "JW" for Jordan-Wigner or "BK" for Bravyi-Kitaev.
  The default is "JW" (Jordan-Wigner).
- **`boson_to_qubit_transform`**: Currently the only available option is "binary" (the default),
  so this setting currently has no effect.
- **`max_bosons_per_state`**: Specifies the maximum number of bosons that can exist in a single
  bosonic state.  Formally an infinite number of bosons is permitted in each bosonic state, but
  this must be truncated to a finite value for computation.  For encodings such as "binary" that
  cannot provide an arbitrary upper limit on the number of bosons, this will be rounded up if
  necessary.  That is, the algorithm will be able to represent _at least_ this many bosons per
  bosonic state.  **This parameter is required for systems containing bosons** and has no default
  value -- the script will raise an error if bosons are present and this is not specified.

**Example usage**:
```python
# Load a NumPy file with custom transformation
hamiltonian.load_second_quantization("file.npz",
                                     fermion_to_qubit_transform="BK",
                                     boson_to_qubit_transform="binary",
                                     max_bosons_per_state=10)

# Load an HDF5 file with default Jordan-Wigner transformation
hamiltonian.load_second_quantization("data.h5")
```

### Encoding as a Unitary

Currently all of our applications involve encoding the Hamiltonian ($\hat{H}$) as a time-evolution
unitary
$$ \hat{U} = e^{i \hat{H} t / \hbar} $$
There are a variety of ways to encode a Hamiltonian into a unitary matrix.  Currently this script
supports

- Trotterization: The function **`unitary.encode_ramped_trotter()`** uses a Trotter formula to
  encode the Hamiltonian into a time-evolution unitary.  It takes the following arguments:
  - `timestep`: If not provided, QHAT will attempt to pick the timestep that provides the most
    efficient algorithm while still preventing aliasing of phases.  Providing a timestep will
    override this with a user-selected value.
  - `energy_error`: The maximum error allowed from the Trotterization process.  If not provided,
    the script will generate an error.
  - `error_scale`: This option is deprecated.
  - `phase_scale_factor`: (Optional, default: 1.01) Scales the energy range by this factor when
    computing the time evolution parameter. This ensures eigenvalue phases never hit exactly ±π,
    avoiding the aliasing ambiguity at the complex exponential's branch cut. Values slightly larger
    than 1.0 (e.g., 1.01) are recommended. The default maps phases to approximately [-3.11, 3.11]
    instead of [-π, π].
  - `trotter_implementation`: (Optional) Choose between two Trotter implementations:
    - `"flattened"` (default, recommended): Flattened QHAT implementation with flat expansion and
      optional term combining. Term combining reduces operation count (benefit varies by
      Hamiltonian).  Uses `CommutingPauliStringEvolution` internally, enabling future grouping of
      commuting terms.  This implementation has comprehensive test coverage (57 tests) and is
      recommended for all use cases.
    - `"original"`: Original QHAT implementation with nested bloq structure
      (RampedTrotterizedUnitary → RampedTrotterStep → TrotterRamp). Warning: This implementation
      has little to no test coverage and may not work correctly in all cases. Use only for exact
      reproduction of earlier results or debugging comparisons.
  - `trotter_combine_terms`: (Optional, only for `trotter_implementation="flattened"`)
    - `True` (default): Combine adjacent identical Pauli string terms for efficiency. The reduction
      in operation count varies depending on the Hamiltonian structure.
    - `False`: Keep all terms separate. Useful for comparing results with the original
      implementation. When disabled, produces the same gate counts as the original implementation.
  - `trotter_order`: (Optional) Explicitly specify which Trotter order to use:
    - `None` (default): Automatically select between first-order and second-order based on
      cost-effectiveness
    - `"first order"`: Force use of first-order Trotter
    - `"second order"`: Force use of second-order Trotter
    - `"fourth order"`: Force use of fourth-order Trotter using five-term Suzuki recursion (step
      count is estimates from second-order method, hence being excluded from automatic selection).
- Double-Factorization: The function **`unitary.encode_double_factorization()`** uses a
  double-factorized block-encoding of the Hamiltonian.  This model is preliminary: it has not been
  verified and is known to fail unexpectedly for some Hamiltonians.  It takes the following
  arguments:
  - `energy_error`: The maximum error allowed from the double-factorization process.  If not
    provided, the script will generate an error.

Calling more than one of these functions will generate an error.

### Generating an Algorithm

A variety of algorithms can be generated and analyzed based on the time-evolution unitary operator.
However, this section is still under development and currently has very limited options:

- Textbook Phase Estimation: Setting `algorithm.method` to "QPE: qualtran textbook" will embed the
  unitary encoding of the Hamiltonian into a phase estimation algorithm that uses the classic
  "textbook" method (see, for example, Nielson and Chuang's "Quantum Computation and Quantum
  Information").
- Qubitized Phase Estimation: Setting `algorithm.method` to "QPE: pyliqtr qubitized" will embed the
  unitary encoding of the Hamiltonian into a phase estimation algorithm that uses pyLIQTR's
  qubitized phase estimation.  This uses only a single ancilla qubit for the phase, with multiple
  measurements to extract the necessary number of bits of information.  This method only works with
  qubitized encodings such as double-factorization.  The integration of this method into the larger
  workflow has not yet been verified, so use this method with caution.
- Time Evolution: Setting `algorithm.method` to "time evolution" will return the time evolution
  unitary operator.  This is useful for analyzing the resource requirements of the unitary itself
  or for building custom algorithms.
- Controlled Time Evolution: Setting `algorithm.method` to "controlled time evolution" will return
  a singly-controlled version of the time evolution unitary.  This adds one control qubit to the
  algorithm and is useful for building custom quantum algorithms that require controlled time
  evolution operators (for example, when building phase estimation or iterative phase estimation
  algorithms manually).

When performing phase estimation, it is necessary to set the number of phase qubits (which, in the
qubitized method, translates to the number of measurements of the single phase qubit).  This can be
controlled directly by the user by setting `algorithm.num_phase_qubits`.  But this can also be
computed by the script by setting

- `algorithm.energy_error`: The maximum energy error permitted from phase estimation.
- `algorithm.probability_of_failure`: The maximum probability of measuring the wrong phase at the
  end of the algorithm.

### Analyzing an Algorithm

There are many details of the algorithm that may be worth analyzing. The available analyses are
discussed below.

#### Resource Estimation

- **pyLIQTR Resource Estimation**: Setting `analysis.resource_estimator` to "pyLIQTR" will use the
  resource estimation capability from pyLIQTR, which is in turn based on the resource estimation
  capability from Qualtran.
- **Cirq Resource Estimation**: Setting `analysis.resource_estimator` to "Cirq" is available but
  deprecated and may not work correctly.

#### Matrix and Eigendecomposition Output

**Full-Algorithm Unitary Matrix Output**: Setting `analysis.algorithm_matrix_output_file` to a
filename will compute and save the unitary matrix representation of the full algorithm. Supported
formats:
- `.npz`: NumPy compressed format
- `.h5` or `.hdf5`: HDF5 format with compression
- `.txt`, `.dat`, or `.csv`: Human-readable sparse text format

The matrix file includes metadata such as git hash, timestamp, unitarity error, and matrix norm.

**Example**:
```python
analysis.algorithm_matrix_output_file = "unitary_matrix.npz"
```

**Flexible Operator Output**: Use `analysis.save_operator_to_file()` to save exact or approximate
operators in either matrix or eigendecomposition form.

This is useful for:
- Validating approximate algorithms by comparing exact vs approximate eigenvalues
- Computing exact ground state energies for small systems
- Testing and debugging algorithm implementations
- Analyzing the full spectrum of eigenenergies

**Parameters**:
- `filename`: the name of the file to create, extension defines the file format
- `source`: `'exact'` (computed from the Hamiltonian) or `'approximate'` (computed from the
  approximate unitary encoding of the Hamiltonian)
- `operator_type`: `'hamiltonian'` for H or `'time_evolution'` for U = exp(-iHt/ℏ)
- `energy_shifted`: `False` for the physical energy scale or `True` to include the energy-shift
  applied to reduce the circuit depth
- `representation`: `'matrix'` or `'eigendecomposition'`

**Example (matrix output)**:
```python
analysis.save_operator_to_file(
    filename='H_approx.npz',
    source='approximate',
    operator_type='hamiltonian',
    energy_shifted=False,
    representation='matrix'
)
```

**Example (eigendecomposition output)**:
```python
analysis.save_operator_to_file(
    filename='H_exact_eig.npz',
    source='exact',
    operator_type='hamiltonian',
    energy_shifted=False,
    representation='eigendecomposition'
)
```

**Key features of eigendecomposition**:
- **Always computes full spectrum** (all eigenstates)
- **Only feasible for small systems**
- Output contains eigenvalues and eigenvectors

**Note**: For large systems, the exact matrix computation creates a matrix-free LinearOperator
that can be used with scipy sparse eigensolvers, but cannot be directly saved to a file. The
analysis will skip file output and note this in the results.

#### Error Analysis

- **Error Analysis**: Enables three types of error metrics comparing exact and approximate
  representations. Each error type is independently enabled by setting its corresponding
  configuration parameter.

  **Configuration parameters**:
  
  - **`enable_eigenvalue_errors`**: Enable eigenenergy comparison (default: False, disabled)
    - **Note**: Parameter name uses "eigenvalue" but specifically compares **eigenenergies**
      (unshifted H eigenvalues, not U eigenvalues, not energy-shifted)
    - Set to `True` to compute errors for ALL eigenstates (full spectrum)
    - Compares eigenenergies element-wise after both are sorted by energy
    - Ground state (exact) compared with ground state (approximate), etc.
    - Does not compare eigenvectors
    - **Best for**: Validating eigenenergy accuracy across the full spectrum

  - **`error_matrix_norms`**: Which matrix norms to compute (default: None, disabled)
    - Single string: `"frobenius"` or `"spectral"`
    - List: `["frobenius", "spectral"]` for both
    - **Frobenius norm**: Element-wise difference, fast to compute
      - ||H_exact - H_approx||_F = sqrt(sum of squares of all elements)
      - Good for quick comparisons
    - **Spectral norm**: Worst-case effect on any quantum state, physically meaningful
      - ||H_exact - H_approx||_2 = largest singular value
      - More expensive to compute, especially for large systems
    - **For large systems**: Uses matrix-free computation
      - Frobenius: Requires 2^Q matrix-vector products
      - Spectral: Uses power iteration (can take longer)
      - Progress warnings displayed during computation
    - **Best for**: Physical bounds on algorithm error

  - **`error_state_inputs`**: State files for state-dependent errors (default: None, disabled)
    - Single filename (string): `"ground_state.npy"`
    - Multiple filenames (list): `["state1.npy", "state2.npy"]`
    - Compares: ||H_exact|ψ⟩ - H_approx|ψ⟩||
    - **Best-scaling error metric** for large systems
      - Only requires O(2^Q) memory (state vectors), not O(2^(2Q)) (matrices)
      - Fast: just applies operators to states
    - **Best for**: Error on specific physically relevant states

  **Output file**: `error_analysis.npz` containing all computed error metrics

  **Examples**:
  ```python
  # Eigenvalue error comparison (compares all eigenvalues from eigendecomposition)
  analysis.enable_eigenvalue_errors = True

  # Matrix norm errors (physical bounds)
  analysis.error_matrix_norms = "frobenius"  # Fast
  # or
  analysis.error_matrix_norms = ["frobenius", "spectral"]  # Both

  # State-dependent errors (best scaling)
  analysis.error_state_inputs = "ground_state.npy"
  # or multiple states
  analysis.error_state_inputs = ["ground.npy", "excited.npy"]

  # All three error types together (comprehensive validation)
  analysis.enable_eigenvalue_errors = True
  analysis.error_matrix_norms = "frobenius"
  analysis.error_state_inputs = ["ground.npy", "excited.npy"]
  ```

  **Note**: Error analysis is currently only available for Trotter-based unitary encodings.
  If you request error analysis with other encoding methods (e.g., pauli-lcu, double-factorization),
  the analysis will be automatically disabled with a warning message.

  **When to use each error type**:
  - **Eigenvalue errors**: Use when you want to validate eigenvalues computed in the
    eigendecomposition
  - **Matrix norm errors**: Use for physical bounds, but expensive for large systems
  - **State errors**: Use for specific physically relevant states, scales best

#### Numerical Simulation (Approximate)

- **Numerical Simulation**: Setting `analysis.numerical_simulation_inputs` to one or more state
  vector files will apply the algorithm to the input state(s) via numerical simulation, producing
  output state(s).
  
  **Input format**: NumPy `.npy` format (compatible with `hamgen.py` output). Input can be:
  - Single filename (string): `"initial_state.npy"`
  - Multiple filenames (list): `["state1.npy", "state2.npy", "state3.npy"]`
  
  **Output naming**: Automatic suffix `_final` is added to input filename:
  - `initial_state.npy` → `initial_state_final.npy`
  
  **Example**:
  ```python
  # Single state simulation
  analysis.numerical_simulation_inputs = "initial_state.npy"
  
  # Multiple states
  analysis.numerical_simulation_inputs = [
      "ground_state.npy",
      "excited_state.npy",
      "superposition.npy"
  ]
  ```
  
  **Creating input states**: State vectors must be 1D complex NumPy arrays with dimension 2^Q
  (where Q is the number of qubits):
  ```python
  import numpy as np
  
  # Create 4-qubit state |0000⟩
  n_qubits = 4
  psi = np.zeros(2**n_qubits, dtype=complex)
  psi[0] = 1.0
  np.save("initial_state.npy", psi)
  ```

## Generated Files

The script will print a log both to the screen and to a logfile. It also generates output files
depending on the analyses requested. If `general.output_directory` is set, all output files are
written to that directory.

- **Log file**: Default `analysis.log`, configurable via `general.logfile`
- **TOML summary**: Hash-based filename (e.g., `12345678901234567890.toml`) containing inputs and
  results. Shows final results but not intermediate values.
- **Matrix file** (if `algorithm_matrix_output_file` specified): Unitary matrix in specified format
  (`.npz`, `.h5`, or `.txt`)
- **Matrix files** (if `save_operator_to_file` with `representation='matrix'` used): Exact and/or
  approximate Hamiltonian and/or time-evolution operator as a matrix
- **Eigendecomposition files** (if `save_operator_to_file` with `representation='eigendecomposition'`
  used): Eigendecomposition of corresponding matrices
- **Error analysis file** (if any error analysis enabled): `error_analysis.npz`
- **Final state files** (if `numerical_simulation_inputs` specified): Evolved quantum states with
  `_final` suffix (e.g., `initial_state.npy` → `initial_state_final.npy`)

**Note**: All file paths respect `general.output_directory` if set. This keeps outputs organized
when running multiple analyses.

The logfile is typically most useful for understanding the analysis process and intermediate
values.

## Examples

The `analysis/examples/` directory contains configuration files demonstrating various analysis
capabilities:

### Basic Example: `config.py`

The basic `config.py` file presents a simple configuration for generating resource estimates. It
loads data from the tensors file in the `examples` directory and demonstrates switching between
Trotterization-based and double-factorization-based analysis (showing that configuration files are
Python scripts, not just key-value lists).

### Complete Options Reference: `all_analyses.config`

**This is the definitive reference for ALL available options in the QHAT analysis tool.**

The `all_analyses.config` file demonstrates every available feature and configuration parameter:

1. **Resource Estimation**: Quantum gate counts, qubit requirements, circuit depth
2. **Flexible Operator Output**: Save any combination of operator forms (16 different combinations)
   - Matrix or eigendecomposition representation
   - Exact or approximate operators
   - Hamiltonian or time-evolution form
   - Unshifted (physical) or shifted (QPE) energy scales
3. **Full Algorithm Circuit**: Complete algorithm matrix (U_approx for time evolution, QPE circuit
   for QPE)
4. **Numerical Simulation**: Apply operators to quantum states
5. **Error Analysis**: Eigenvalue errors, matrix norm errors, state-dependent errors
6. **All Configuration Options**: Hamiltonian loading, encoding methods, algorithm selection

**NOTE**: This file produces ~18 output files and is computationally expensive. It is intended as a
comprehensive reference for developers and users learning the tool, not for production use.

For practical workflows, see `config.py` (basic) or `error_analysis_tutorial.config` (error
analysis tutorial).

### Error Analysis Tutorial: `error_analysis_tutorial.config`

**A pedagogical tutorial teaching error analysis from first principles.**

This tutorial configuration teaches you how to use QHAT's error analysis features to quantify
the accuracy of approximate quantum algorithms. It provides:

- **Clear learning objectives** and step-by-step progression
- **Three types of error analysis** with detailed explanations:
  - Eigenvalue errors (energy accuracy)
  - Matrix norm errors (global operator accuracy)
  - State-dependent errors (specific state accuracy)
- **Interpretation guidance** for understanding results
- **Experimentation suggestions** to learn the accuracy vs. cost trade-off

Recommended for: New users learning error analysis, anyone needing to validate algorithm accuracy
- Process multiple input states
- Structure configuration files with clear documentation
- Set up for future features (exact matrices, eigendecomposition, error analysis)

To run either example:
```bash
python3.11 -m qhat.analysis.driver config.py
# or
python3.11 -m qhat.analysis.driver examples/config_full_analysis.py
```

## L-sweep Trotter factorization benchmark

`benchmark_L_sweep_trotter.py` compares four schedules against the same
identity-free Jordan--Wigner Hamiltonian, initial state, evolution time,
formula order, step count, and error definitions:

- `jw_raw`: construct the full fermionic Hamiltonian, JW-map it, combine
  identical Pauli strings, remove identity, and exponentiate every remaining
  Pauli string in OpenFermion insertion order.
- `jw_coloring`: use the same final Pauli-string factors, ordered by greedy
  `largest_first` coloring of their noncommutation graph.
- `fermionic_coloring`: the backward-compatible, legacy
  **fermionic-induced JW ordering**. Complete Hermitian fermionic terms are
  colored first, but that order only induces a permutation of the final
  combined JW Pauli strings. Its factors are still individual Pauli strings.
- `fermionic_term_coloring`: color complete Hermitian fermionic terms, map
  each ordered complete term through JW, and retain its complete mapped dense
  Hermitian matrix as one factor. Its factor exponential is
  `expm(-1j * dt * JW(F_j))`, with the identity/global-phase component omitted
  consistently. The mapped Pauli summands inside a factor are never
  separately exponentiated.

The two pipelines are:

```text
Pauli level:
fermionic H -> JW -> combine Pauli keys -> remove identity
            -> [c_k P_k factors] -> raw or Pauli-graph order

Fermionic-term level:
fermionic H -> pair monomial + adjoint -> remove fermionic identity
            -> fermionic noncommutation graph -> color/order complete F_j
            -> JW(F_j) dense matrices -> remove global-phase components
            -> [complete mapped F_j factors]
```

Jordan--Wigner mapping alone does not determine factorization level. The
factor boundary does: `fermionic_term_coloring` establishes factors before
mapping and preserves each mapped sum as one exponential, whereas the other
methods exponentiate final Pauli strings.

The implementation validates monomial/adjoint pairing, Hermiticity of every
complete term and mapped matrix, reconstruction of the non-identity fermionic
Hamiltonian, and reconstruction of the exact identity-free JW matrix. It uses
the analytic sine/cosine exponential for Pauli factors and SciPy's dense
matrix exponential for complete fermionic factors. The same chronological
left-multiplication convention and the same second-order slice are used by
first-, second-, and five-stage fourth-order Suzuki formulas.

### CSV schema and resume behavior

New rows use schema version `2`. All historical fields remain in their
original order; the following semantic and comparison fields are appended:

```text
factorization_level
ordering_method
mapping
number_of_factors
factor_exponential_type
factor_reconstruction_error
graph_level
schema_version
operator_error_ratio_to_jw_coloring
state_infidelity_ratio_to_jw_coloring
operator_error_ratio_to_fermionic_coloring
state_infidelity_ratio_to_fermionic_coloring
```

`ordering` remains the canonical resume key component, so historical
`fermionic_coloring` rows cannot collide with `fermionic_term_coloring` rows.
On `--resume`, the benchmark reads both old and v2 rows. A CSV with the exact
historical header is migrated atomically to the appended v2 header; inferred
legacy semantics are recorded with `schema_version=1`, and old values are
preserved. A reordered or otherwise incompatible header is rejected instead
of being blindly appended. `jw_raw` is automatically added to an explicit
method subset because the legacy `*_ratio_to_jw_raw` columns are retained.

### Reproducible runs

A cheap run on the smallest discovered tensor case is:

```bash
python analysis/benchmark_L_sweep_trotter.py \
  --library hamiltonian_generator/library \
  --limit 1 \
  --orderings jw_raw jw_coloring fermionic_coloring fermionic_term_coloring \
  --formula-orders 1 2 \
  --steps 1 2 \
  --time 1 \
  --skip-operator-norm \
  --output /tmp/qhat_l_sweep_smoke.csv
```

Add `--formula-orders 4` for a fourth-order smoke run. Re-run the same command
with `--resume` to verify completion skipping. An explicit-method check is:

```bash
python analysis/benchmark_L_sweep_trotter.py \
  --library hamiltonian_generator/library --limit 1 \
  --orderings fermionic_term_coloring \
  --formula-orders 1 --steps 1 --time 1 \
  --output /tmp/qhat_l_sweep_partial.csv
```

The historical full matrix contains 153 cases, times 1/2/5, formula orders
1/2/4, and steps 10/100/1000. The resume-safe four-method rerun is:

```bash
for qhat_time in 1 2 5; do
  python analysis/benchmark_L_sweep_trotter.py \
    --library hamiltonian_generator/library \
    --orderings jw_raw jw_coloring fermionic_coloring fermionic_term_coloring \
    --formula-orders 1 2 4 \
    --steps 10 100 1000 \
    --time "$qhat_time" \
    --output "analysis/l_sweep_trotter_state_t${qhat_time}_v2.csv" \
    --resume
done
```

That matrix contains `153 * 3 * 3 * 3 * 4 = 16,524` successful rows when
every case completes.

Dense exact evolution stores `2^n x 2^n` matrices. General matrix
exponentials cost approximately `O(2^(3n))` time and `O(2^(2n))` memory per
factor, so the genuine fermionic-term method is substantially more expensive
than the analytic Pauli path. No case is silently omitted and no separate
qubit limit is imposed; use `--case-pattern` and `--limit` to stage large
runs. `--skip-operator-norm` avoids repeated spectral-norm calculations but
does not remove dense factor exponentials.

### Plotting old and new results

Both visualization commands infer missing factorization metadata in old
three-method CSV files and also recognize four-method v2 files:

```bash
python analysis/visualize_L_sweep_trotter.py \
  --csv /tmp/qhat_l_sweep_smoke.csv \
  --case-id H-H_0.70_sto-6g_as-001-001 \
  --output-dir /tmp/qhat_l_sweep_smoke_figures

python analysis/visualize_all_L_sweep_cases.py \
  --csv analysis/l_sweep_trotter_state_t1_v2.csv \
        analysis/l_sweep_trotter_state_t2_v2.csv \
        analysis/l_sweep_trotter_state_t5_v2.csv \
  --output-dir analysis/l_sweep_all_case_figures_v2
```

Curves and summaries tolerate method subsets. Legends state the factorization
level, and graph comparisons explicitly label fermionic-term graph vertices
versus Pauli-string graph vertices rather than treating them as the same
quantity.
