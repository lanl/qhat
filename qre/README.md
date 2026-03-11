# QHAT Quantum Resource Estimation

[[_TOC_]]

This directory contains the modules for performing quantum resource estimation (QRE).  This is run
by
```
python qre_driver.py
```
By default, this will load the setup from `config.py`, but this can be changed by passing the name
of a configuration file as a command line argument
```
python qre_driver.py my_configuration_file.py
```

## Configuration Options

Configuration files are themselves Python scripts, allowing users to use control logic to build up
complex configuration files.

Configuration is broken down by several parts of the processing script.

### General

"General" configuration governs the behavior of the resource analysis script itself.  Currently the
only controls under this heading relate to logging progress of the script.

You can configure the log file that the script will write to by setting **`general.logfile`** to
the name of the logfile you want to use.  The default is `qre.log`.

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

- HDF5: HDF5 files can be loaded by calling the function **`hamiltonian.load_hdf5()`**, which takes
  as argument the name of the HDF5 file.  Currently the format assumes that there is no constant
  term, the one-body term is called "1e", and the two-body term is called "2e".
- NumPy: NumPy-generated files (.npy or .npz) can be loaded by calling the function
  **`hamiltonian.load_numpy()`**, which takes as argument the name of the NPY or NPZ file.
  Currently the format is assumed to include:
  - a constant called "constant" (optional)
  - a one-body tensor called "one-body" (required)
  - a two-body tensor called "two-body" (required)
  - a scalar for bosons called "bosonic\_scalar" (optional)
  - a tensor for fermion-boson interactions called "fb\_interaction" (optional)

  If the "bosonic\_scalar" and "fb\_interaction" options are not present, this is the standard
  fermionic second-quantization format.  The scripts in the `hamiltonian_generator` directory
  create Hamiltonians in this format, but provide only the one-body and two-body tensors.

After loading a Hamiltonian, the user may have knowledge regarding the bounds on the eigenvalues of
the Hamiltonian, which can be used to optimize certain parts of the resource analysis.  In order to
specify these, use the functions **`hamiltonian.set_energy_lower_bound()`** and
**`hamiltonian.set_energy_upper_bound()`**.  These functions both take an argument specifying the
value of the energy (eigenvalue) bound.  There is an optional second argument specifying whether to
use this value exactly (`exact=True`) or only use it if the automatic bound estimate is a looser
constraint than the bound provided here (`exact=False`).  The last call to each function takes
precedence over prior calls to the same function.

Certain formats load Hamiltonians in second-quantization format and can use different
fermion-to-qubit and/or boson-to-qubit transformations.  These can be set by the following options:

- Set **`hamiltonian.fermion_to_qubit_transform`** to "JW" for Jordan-Wigner or "BK" for
  Bravyi-Kitaev.  The default is Jordan-Wigner.
- The **`hamiltonian.boson_to_qubit_transform`** option also exists, but currently the only option
  available (which is the default) is "binary", so this option currently has no effect.
- The **`hamiltonian.max_bosons_per_state`** option specifies the maximum number of bosons that can
  exist in a single bosonic state.  Formally an infinite number of bosons is permitted in each
  bosonic state, but this must be truncated to a finite value for computation.  For encodings such
  as "binary" that cannot provide an arbitrary upper limit on the number of bosons, this will be
  rounded up if necessary.  That is, the circuit will be able to represent _at least_ this many
  bosons per bosonic state.  There is no default value, so if your system contains bosons this
  option is required.

### Encoding as a Unitary

Currently all of our applications involve encoding the Hamiltonian ($\hat{H}$) as a time-evolution
unitary
$$ \hat{U} = e^{i \hat{H} t / \hbar} $$
There are a variety of ways to encode a Hamiltonian into a unitary matrix.  Currently this script
supports

- Trotterization: The function **`unitary.encode_ramped_trotter()`** uses a Trotter formula to
  encode the Hamiltonian into a time-evolution unitary.  It takes the following arguments:
  - `timestep`: If not provided, QHAT will attempt to pick the timestep that provides the most
    efficient circuit while still preventing aliasing of phases.  Providing a timestep will
    override this with a user-selected value.
  - `energy_error`: The maximum error allowed from the Trotterization process.  If not provided,
    the script will generate an error.
  - `error_scale`: This option is deprecated.
- Double-Factorization: The function **`unitary.encode_double_factorization()`** uses a
  double-factorized block-encoding of the Hamiltonian.  This model is preliminary: it has not been
  verified and is known to fail unexpectedly for some Hamiltonians.  It takes the following
  arguments:
  - `energy_error`: The maximum error allowed from the double-factorization process.  If not
    provided, the script will generate an error.

Calling more than one of these functions will generate an error.

### Generating a Circuit

A variety of circuits can be generated and analyzed based on the time-evolution unitary operator.
However, this section is still under development and currently has very limited options:

- Textbook Phase Estimation: Setting `circuit.method` to "qualtran textbook" will embed the unitary
  encoding of the Hamiltonian into a phase estimation circuit that uses the classic
  "textbook" method (see, for example, Nielson and Chuang's "Quantum Computation and Quantum
  Information").
- Qubitized Phase Estimation: Setting `circuit.method` to "pyliqtr qubitized" will embed the
  unitary encoding of the Hamiltonian into a phase estimation circuit that uses pyLIQTR's qubitized
  phase estimation.  This uses only a single ancilla qubit for the phase, with multiple
  measurements to extract the necessary number of bits of information.  This method only works with
  qubitized encodings such as double-factorization.  The integration of this method into the larger
  workflow has not yet been verified, so use this method with caution.

When performing phase estimation, it is necessary to set the number of phase qubits (which, in the
qubitized method, translates to the number of measurements of the single phase qubit).  This can be
controlled directly by the user by setting `circuit.num_phase_qubits`.  But this can also be
computed by the script by setting

- `circuit.energy_error`: The maximum energy error permitted from phase estimation.
- `circuit.probability_of_failure`: The maximum probability of measuring the wrong phase at the end
  of the circuit.

### Analyzing a Circuit

There are many details of the circuit that may be worth analyzing.  At this point, this script is
focused on resource estimation.

- pyLIQTR Resource Estimation: Setting `analysis.resource_estimator` to "pyLIQTR" will use the
  resource estimation capability from pyLIQTR, which is in turn based on the resource estimation
  capability from Qualtran.
- Cirq Resource Estimation: Setting `analysis.resource_estimator` to "Cirq"is available but
  deprecated and may not work correctly.

## Generated Files

The script will print a log both to the screen and to a logfile.  It also generates a TOML file
that summarizes the inputs and results.  The TOML file is based on a hash, so the filename will
likely be a long string of numbers with the `.toml` extension.  The TOML file only shows the final
results and does not report the intermediate values, so the logfile will typically be more useful.

## Example

The provided `config.py` file presents an example configuration file that can be used to generate
resource estimates.  It loads data from the tensors file in the `examples` directory, and uses the
options specified in the configuration file to generate resource estimates.  The first two lines of
the provided `config.py` file allow the user to easily switch between a Trotterization-based
analysis or a double-factorization-based analysis (additionally demonstrating that configuration
files are themselves Python scripts rather than simple key-value lists).
