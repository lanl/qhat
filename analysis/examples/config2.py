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

unitary.encode_ramped_trotter(
        energy_error = 0.5 * energy_error,
        trotter_implementation = "flattened",
        trotter_combine_terms = True,
        ordering_method = "lexicographical"
        )

# _________________________________________________________________________________________________
# algorithm

algorithm.method = "time evolution"

# _________________________________________________________________________________________________
# analysis

analysis.resource_estimator = "pyLIQTR"
analysis.matrix_output_file = "my_algorithm_matrix.npz"

