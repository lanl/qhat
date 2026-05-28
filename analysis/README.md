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

"General" configuration governs the behavior of the resource analysis script itself.  Currently the
only controls under this heading relate to logging progress of the script.

You can configure the log file that the script will write to by setting **`general.logfile`** to
the name of the logfile you want to use.  The default is `analysis.log`.

The other capability under general configuration is to set the log level, which is done by calling
the following functions.  If you call multiple of these functions, whichever is called last takes
precedence.

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
  - `trotter_implementation`: (Optional) Choose between two Trotter implementations:
    - `"flattened"` (default, recommended): Flattened QHAT implementation with flat expansion and
      optional term combining. Term combining reduces operation count (benefit varies by Hamiltonian).
      Uses `CommutingPauliStringEvolution` internally, enabling future grouping of commuting terms.
      This implementation has comprehensive test coverage (57 tests) and is recommended for all use
      cases.
    - `"original"`: Original QHAT implementation with nested bloq structure
      (RampedTrotterizedUnitary → RampedTrotterStep → TrotterRamp). Warning: This implementation
      has little to no test coverage and may not work correctly in all cases. Use only for exact
      reproduction of earlier results or debugging comparisons.
  - `trotter_combine_terms`: (Optional, only for `trotter_implementation="flattened"`)
    - `True` (default): Combine adjacent identical Pauli string terms for efficiency. The reduction
      in operation count varies depending on the Hamiltonian structure.
    - `False`: Keep all terms separate. Useful for comparing results with the original
      implementation. When disabled, produces the same gate counts as the original implementation.
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

There are many details of the algorithm that may be worth analyzing. The available analyses include:

#### Resource Estimation

- **pyLIQTR Resource Estimation**: Setting `analysis.resource_estimator` to "pyLIQTR" will use the
  resource estimation capability from pyLIQTR, which is in turn based on the resource estimation
  capability from Qualtran.
- **Cirq Resource Estimation**: Setting `analysis.resource_estimator` to "Cirq" is available but
  deprecated and may not work correctly.

#### Matrix Output

- **Unitary Matrix Output**: Setting `analysis.matrix_output_file` to a filename will compute and
  save the full unitary matrix representation of the algorithm. Supported formats:
  - `.npz`: NumPy compressed format
  - `.h5` or `.hdf5`: HDF5 format with compression
  - `.txt`, `.dat`, or `.csv`: Human-readable sparse text format
  
  The matrix file includes metadata such as git hash, timestamp, unitarity error, and matrix norm.
  
  **Example**:
  ```python
  analysis.matrix_output_file = "unitary_matrix.npz"
  ```

- **Exact Hamiltonian Matrix Output**: Setting `analysis.exact_matrix_output_file` to a filename
  will compute and save the exact matrix representation of the Hamiltonian **without any
  approximations** (no Trotter, no double-factorization). Supported formats are the same as for
  unitary matrix output.
  
  This is useful for:
  - Validating approximate algorithms by comparing exact vs approximate eigenvalues
  - Computing exact ground state energies for small systems
  - Testing and debugging algorithm implementations
  
  **System size considerations**:
  - **Small systems (≤15 qubits)**: Full dense matrix is computed and saved to file
  - **Large systems (>15 qubits)**: Matrix-free operator is created but not saved (too large)
  
  **Example**:
  ```python
  analysis.exact_matrix_output_file = "exact_hamiltonian.npz"
  ```
  
  **Note**: For large systems, the exact matrix computation creates a matrix-free LinearOperator
  that can be used with scipy sparse eigensolvers, but cannot be directly saved to a file. The
  analysis will skip file output and note this in the results.

#### Numerical Simulation

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
  
  **Creating input states**: State vectors must be 1D complex NumPy arrays with dimension 2^n
  (where n is the number of qubits):
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
depending on the analyses requested:

- **Log file**: Default `analysis.log`, configurable via `general.logfile`
- **TOML summary**: Hash-based filename (e.g., `12345678901234567890.toml`) containing inputs and
  results. Shows final results but not intermediate values.
- **Matrix file** (if `matrix_output_file` specified): Unitary matrix in specified format (`.npz`,
  `.h5`, or `.txt`)
- **Final state files** (if `numerical_simulation_inputs` specified): Evolved quantum states with
  `_final` suffix (e.g., `initial_state.npy` → `initial_state_final.npy`)

The logfile is typically most useful for understanding the analysis process and intermediate values.

## Examples

The `analysis/examples/` directory contains configuration files demonstrating various analysis capabilities:

### Basic Example: `config.py`

The basic `config.py` file presents a simple configuration for generating resource estimates. It loads data from the tensors file in the `examples` directory and demonstrates switching between Trotterization-based and double-factorization-based analysis (showing that configuration files are Python scripts, not just key-value lists).

### Comprehensive Example: `config_full_analysis.py`

The `config_full_analysis.py` file demonstrates **ALL currently available analysis capabilities**:

1. **Resource Estimation**: Quantum gate counts, qubit requirements, circuit depth
2. **Matrix Output**: Save the unitary matrix representation to various formats (.npz, .h5, .txt)
3. **Numerical Simulation**: Apply the unitary to one or more input quantum states

This comprehensive example serves as a complete reference showing how to:
- Configure all available analyses in one file
- Use different output formats
- Process multiple input states
- Structure configuration files with clear documentation
- Set up for future features (exact matrices, eigendecomposition, error analysis)

To run either example:
```bash
python3.11 -m qhat.analysis.driver examples/config.py
# or
python3.11 -m qhat.analysis.driver examples/config_full_analysis.py
```
