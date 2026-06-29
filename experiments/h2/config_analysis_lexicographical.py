energy_error = meV_to_Hartree(1e4)

general.print_verbose()
general.logfile = "experiments/h2/analysis_lexicographical.log"

hamiltonian.load_second_quantization(
    "experiments/h2/h2_0.7414_sto3g_as-002-002.tensors.npz",
    fermion_to_qubit_transform="JW",
)

unitary.encode_ramped_trotter(
    energy_error=0.5 * energy_error,
    trotter_implementation="flattened",
    trotter_combine_terms=False,
    ordering_method="lexicographical",
)

unitary.error_coeff_mode = "exact"

algorithm.method = "time evolution"

analysis.resource_estimator = "pyLIQTR"
