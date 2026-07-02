# QHAT Hamiltonian Generation

[[_TOC_]]

The Hamiltonian generation tool takes a description of a simple molecule and generates a
Hamiltonian that can be used with other QHAT analysis tools.  This is run by
```
python hamgen.py
```
By default, this will load the setup from `config.py`, but this can be changed by passing the name
of a configuration file as a command line argument
```
python hamgen.py my_configuration_file.py
```
This script is a wrapper around the Python package `pySCF`.

## Configuration Options

Configuration files are themselves Python script, allowing users to use control logic to build up
complex configuration files.

### General

"General" configuration governs the behavior of the Hamiltonian generation script itself.

This script generates a number of files, and all of them are named based on a common "file stub".
This is set by **`general.file_stub`**, which has no default and must be set.

The final result of the Hamiltonian generation is to express the Hamiltonian in terms of Pauli
strings.  There are a variety of formats for writing Pauli strings, and this script supports two
formats:
- Set **`general.file_format`** to "hamlib" to get a format that matches what is used by HamLib.
- Set **`general.file_format`** to "default" to get the default Pauli string format.

This script write progress notes to a log file.  The name can be customized by setting
**`general.logfile`**.  By default, the log file is called "hamgen.log".

Users can set the log level, which is done by calling the following functions.  If you call
multiple of these functions, whichever is called last takes precedence.

- **`general.print_verbose()`** -- Calling this function increases the printouts to be more
  verbose, providing additional information as the Hamiltonian generation progresses.  This
  information may be of interest to users, depending on how much detail they want.
- **`general.print_debug()`** -- Calling this function increases the printouts beyond even
  "verbose".  The additional information provided by "debug" printing is typically relevant to
  developers more than users.
- **`general.print_default()`** -- Calling this function resets the printout level back to the
  standard verbosity, removing "verbose" and "debug" printouts.

### Describing the Hamiltonian

Users need to describe the system for which a Hamiltonian will be generated, as well as provide
some details about how to represent the Hamiltonian.

The geometry of the molecule is managed by the **`hamiltonian.add_atom()`** function, which takes
as arguments the atom to be placed followed by the x-, y-, and z-coordinates of the atom.  This
script uses the Python package `mendeleev` to interpret the identity of a given element, so formats
supported by `mendeleev` are also supported by this script.

The atomic orbital basis functions used to construct molecular orbitals are specified by setting
**`hamiltonian.basis`**.  The default is "sto-3g".  The pySCF package provides some basis sets
automatically.  If your Python environment includes the `basis-set-exchange` package, then you have
access to any basis set from the online Basis Set Exchange tool.

In order to specify your active space, set **`hamiltonian.num_active_occupied`** and
**`hamiltonian.num_active_vacant`**.  Because this script leverages pySCF for computing
Hamiltonians, there are restrictions on the sizes of active space.  This script will attempt to
provide helpful error messages, but see the pySCF documentation for more detail.

The fermion-to-qubit transform can be set through the **`hamiltonian.f2q_mapping`** setting.
Currently allowed values are
- for Jordan-Wigner: "Jordan-Wigner", "Jordan Wigner", "JW", "jordan-wigner", "jordan\_wigner",
  "jordan wigner", "jw"
- for Bravyi-Kitaev: "Bravyi-Kitaev", "Bravyi Kitaev", "BK", "bravyi-kitaev", "bravyi\_kitaev",
  "bravyi kitaev", or "bk"

### SymUCCSD: Symmetry Reduction (Optional)

QHAT supports symmetry-reduced Hamiltonian generation using molecular point-group symmetry, based on the SymUCCSD method from Cao et al. (Physical Review A 105, 062452, 2022). This can significantly reduce the number of operators in the Hamiltonian for molecules with symmetry, enabling analysis of larger systems.

**To enable symmetry reduction:**

1. Enable PySCF symmetry detection:
   ```python
   hamiltonian.enable_symmetry = True
   ```

2. Enable symmetry-based filtering:
   ```python
   hamiltonian.apply_symmetry_reduction = True
   ```

**Example:**
```python
# BeH2 with D2h symmetry
general.file_stub = "BeH2_sym"

hamiltonian.add_atom("Be", 0.0, 0.0, 0.0)
hamiltonian.add_atom("H", 0.0, 0.0, -1.3264)
hamiltonian.add_atom("H", 0.0, 0.0, 1.3264)
hamiltonian.basis = "sto-3g"
hamiltonian.num_active_occupied = 6
hamiltonian.num_active_vacant = 8

# Enable symmetry reduction
hamiltonian.enable_symmetry = True
hamiltonian.apply_symmetry_reduction = True

hamiltonian.f2q_mapping = "Jordan-Wigner"
```

**Supported Point Groups:**
- **Abelian groups** (fully supported): C1, Cs, Ci, C2, C2v, C2h, D2, D2h
- **Non-Abelian groups** (mapped to Abelian subgroups): D∞h/Dooh → D2h, C∞v/Coov → C2v, and others

**Benefits:**
- Reduces Hamiltonian tensor elements by filtering non-symmetry-preserving terms
- Enables resource estimation for larger molecules
- Files are tagged with point group (e.g., `BeH2_as-006-008_symD2h_jw.dat`)

**Limitations:**
- Works best with high-symmetry molecules
- Currently supports only closed-shell (spin-restricted) systems
- Reduction amount depends on molecular symmetry and active space

## Generated Files

This script will generate various files for intermediate and final stages of the process.  These
files are useful when generating a number of related Hamiltonians, as the script will attempt to
load intermediate results and re-use them rather than performing the calculations from the
beginning.

The first stage of the Hamiltonian generation process is to leverage pySCF to perform a
Hartree-Fock calculation of the molecule.  The results of this step will be saved as
"[filestub].pickle".

The second stage is to apply the active space, freezing some electrons in low-lying orbital and
locking some high-energy orbitals as being vacant.  The results of this step will be saved as
"[filestub]\_[astag].pickle", where the "[astag]" is a shorthand notation indicating the number of
active occupied orbitals and the number of active vacant orbitals.  If symmetry reduction is enabled,
the filename will include a symmetry tag: "[filestub]\_[astag]\_sym[PG].pickle" (e.g., "BeH2_as-006-008_symD2h.pickle").
The one-body and two-body tensors from this stage of the calculation will also be saved in a file called
"[filestub]\_[astag].tensors.npz", which can be loaded into the resource estimation software using
the `load_numpy` function.

The third stage is to transform the fermionic creation and annihilation operators to qubit
operators, yielding some metadata and a set of Pauli strings.  This information will be written to
a plaintext file called "[filestub]\_[astag]\_[f2q].dat", where "[f2q]" will be "jw" or "bk".  
With symmetry reduction, the filename will be "[filestub]\_[astag]\_sym[PG]\_[f2q].dat".

A logfile will also be written to, recording the progress of the calculation.  The name for this is
set by the `general.logfile` option.

## Examples

The `config.py` file included provides a simple example that will generate the Hamiltonian for a
single molecule.  Running the `hamgen.py` script with no arguments will load this file and generate
the described Hamiltonian.

The `build_config.py` script is an example of how to build a large suite of configuration files to
build Hamiltonians for many related molecules.  It does not run the Hamiltonian generation script,
but only builds the configuration files that can then be fed into the Hamiltonian generate script.
