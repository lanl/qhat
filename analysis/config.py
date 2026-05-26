my_method = "Trotter"
#my_method = "double-factorization"

# this configuration file assumes equipartition of energy between Trotterization error and phase
# estimation error (see usage of energy_error below)
energy_error = meV_to_Hartree(1e4) # 0.01 keV

# _________________________________________________________________________________________________
# general

general.print_verbose()
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
elif my_method == "double-factorization":
    unitary.encode_double_factorization(energy_error=1.0e-4)
else:
    raise ValueError(f"Invalid value for `my_method`: \"{my_method}\"")

# _________________________________________________________________________________________________
# phase estimation algorithm

if my_method == "Trotter":
    algorithm.method = "QPE: qualtran textbook"
    algorithm.energy_error = 0.5 * energy_error
    algorithm.probability_of_failure = 0.01
elif my_method == "double-factorization":
    algorithm.method = "QPE: pyLIQTR qubitized"
    algorithm.num_phase_qubits = 12
else:
    raise ValueError(f"Invalid value for `my_method`: \"{my_method}\"")

# _________________________________________________________________________________________________
# analysis

analysis.resource_estimator = "pyLIQTR"

