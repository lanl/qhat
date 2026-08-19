"""
QHAT Analysis Configuration: basic example

This configuration demonstrates the use of command-line parameters to select between different
quantum simulation methods without editing the file. This is useful for CI pipelines and automated
testing.

Usage:
  python driver.py                                  # Uses Trotter (default)
  python driver.py -p method=pauli-lcu              # Uses Pauli-LCU
  python driver.py -p method=double-factorization   # Uses double-factorization

Each method automatically outputs to a separate directory (Be-H-trotter/, Be-H-pauli-lcu/, etc.)
to avoid overwriting results. You can override this with -p output_directory=custom_dir
"""

# Method selection: override via command-line with -p method=<name>
my_method = params.get("method", "Trotter")

# this configuration file assumes equipartition of energy between Trotterization error and phase
# estimation error (see usage of energy_error below)
energy_error = meV_to_Hartree(1e4) # 0.01 keV

# _________________________________________________________________________________________________
# general

general.print_verbose()

# Set output directory based on method to keep results separated
# Override with: -p output_directory=custom_dir
default_output_dir = {
    "Trotter": "Be-H-trotter",
    "pauli-lcu": "Be-H-pauli-lcu",
    "double-factorization": "Be-H-double-factorization"
}.get(my_method, "Be-H")

general.output_directory = params.get("output_directory", default_output_dir)

general.logfile = "Be-H.log"

# _________________________________________________________________________________________________
# hamiltonian

hamiltonian.load_second_quantization("examples/Be-H_1.30_sto-6g_as-003-003.tensors.npz")

# _________________________________________________________________________________________________
# unitary encoding

if my_method == "Trotter":
    unitary.encode_ramped_trotter(
            energy_error = 0.5 * energy_error,
            error_scale = 1.0,
            trotter_implementation = "flattened",
            trotter_combine_terms = True,
            ordering_method = "lexicographical"
            )
elif my_method == "pauli-lcu":
    unitary.use_library = "pyLIQTR"
#    unitary.use_library = "qualtran"
    unitary.encode_pauli_lcu(energy_error=1.0e-4)
elif my_method == "double-factorization":
    unitary.encode_double_factorization(energy_error=1.0e-4)
else:
    raise ValueError(f"Invalid value for `my_method`: \"{my_method}\"")

# _________________________________________________________________________________________________
# algorithm

if my_method == "Trotter":
    algorithm.method = "QPE: qualtran textbook"
    algorithm.energy_error = 0.5 * energy_error
    algorithm.probability_of_failure = 0.01
elif my_method == "pauli-lcu":
    algorithm.method = "QPE: pyLIQTR qubitized"
#    algorithm.method = "QPE: qualtran qubitization"
    algorithm.num_phase_qubits = 12
elif my_method == "double-factorization":
    algorithm.method = "QPE: pyLIQTR qubitized"
    algorithm.num_phase_qubits = 12
else:
    raise ValueError(f"Invalid value for `my_method`: \"{my_method}\"")

# _________________________________________________________________________________________________
# analysis

analysis.resource_estimator = "pyLIQTR"
#analysis.resource_estimator = "qualtran"

