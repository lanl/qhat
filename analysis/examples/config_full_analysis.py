"""
Comprehensive Configuration Example for QHAT Analysis

This configuration demonstrates ALL currently available analysis capabilities:
1. Resource estimation (quantum gate counts, qubit requirements)
2. Matrix output (save unitary matrix representation)
3. Numerical simulation (apply unitary to input state vectors)
4. Exact matrix computation (compute exact Hamiltonian matrix)
5. Eigendecomposition analysis (compute eigenvalues and eigenvectors)
6. Error analysis (eigenvalue errors, matrix norm errors, state-dependent errors)
7. Exact numerical simulation (apply exact Hamiltonian to states)

Usage:
    python3.11 -m qhat.analysis.driver examples/config_full_analysis.py
"""

# Energy error budget (equipartition between Trotterization and phase estimation)
energy_error = meV_to_Hartree(1e4)  # 0.01 keV = 10 eV

# =================================================================================================
# GENERAL CONFIGURATION
# =================================================================================================

general.print_verbose()  # Options: print_default(), print_verbose(), print_debug()
general.logfile = "Be-H_full_analysis.log"

# =================================================================================================
# HAMILTONIAN SPECIFICATION
# =================================================================================================

hamiltonian.load_second_quantization("analysis/examples/Be-H_1.30_sto-6g_as-003-003.tensors.npz")

# Optional: Set energy bounds if known (for validation purposes)
# hamiltonian.set_energy_lower_bound(-5.0, exact=False)
# hamiltonian.set_energy_upper_bound(0.0, exact=False)

# =================================================================================================
# UNITARY ENCODING
# =================================================================================================

unitary.encode_ramped_trotter(
    energy_error=0.5 * energy_error,  # Split error budget
    trotter_implementation="flattened",  # Options: "flattened", "recursive"
    trotter_combine_terms=True,
    ordering_method="lexicographical"  # Options: "lexicographical", "random", None
)

# Alternative encoding methods:
# unitary.encode_double_factorization(energy_error=0.5 * energy_error)

# =================================================================================================
# ALGORITHM SELECTION
# =================================================================================================

algorithm.method = "time evolution"  # Options: "time evolution", "QPE: qualtran textbook"

# For QPE-based algorithms:
# algorithm.num_phase_qubits = 10
# algorithm.probability_of_failure = 0.01

# =================================================================================================
# ANALYSIS CONFIGURATION
# =================================================================================================

# -------------------------------------------------------------------------------------------------
# 1. RESOURCE ESTIMATION
# -------------------------------------------------------------------------------------------------

analysis.resource_estimator = "pyLIQTR"  # Options: "pyLIQTR", "cirq", "qualtran"

# -------------------------------------------------------------------------------------------------
# 2. MATRIX OUTPUT
# -------------------------------------------------------------------------------------------------

analysis.matrix_output_file = "Be-H_unitary_matrix.npz"

# Supported formats: .npz (recommended), .h5/.hdf5, .txt (small matrices only)

# -------------------------------------------------------------------------------------------------
# 3. NUMERICAL SIMULATION
# -------------------------------------------------------------------------------------------------

analysis.numerical_simulation_inputs = "analysis/examples/initial_state.npy"

# For multiple states:
# analysis.numerical_simulation_inputs = [
#     "analysis/examples/initial_state.npy",
#     "analysis/examples/excited_state.npy"
# ]

# -------------------------------------------------------------------------------------------------
# 4. EXACT MATRIX COMPUTATION
# -------------------------------------------------------------------------------------------------

analysis.exact_matrix_output_file = "Be-H_exact_hamiltonian.npz"

# Optional: Configure memory threshold for dense/sparse selection
# analysis.matrix_memory_threshold_gb = 16.0  # Default: 16 GB

# -------------------------------------------------------------------------------------------------
# 5. EIGENDECOMPOSITION ANALYSIS
# -------------------------------------------------------------------------------------------------

analysis.num_eigenvalues = 5  # Integer or "all" for full eigendecomposition
analysis.eigendecomposition_matrices = "both"  # Options: "exact", "approximate", "both"
analysis.which_eigenvalues = "smallest"  # Options: "smallest", "largest", "both"

# Configuration options:
#   num_eigenvalues:
#     - 0 (default): Disabled
#     - Positive int: Compute that many eigenvalues
#     - "all": Full eigendecomposition (only for small systems ≤10 qubits)
#
#   eigendecomposition_matrices:
#     - "exact": Only exact Hamiltonian
#     - "approximate": Only unitary matrix
#     - "both": Both (for comparison)
#
#   which_eigenvalues:
#     - "smallest": Ground state + low-lying excited states
#     - "largest": High-energy states
#     - "both": Both ends of spectrum (returns 2k eigenvalues)

# -------------------------------------------------------------------------------------------------
# 6. ERROR ANALYSIS
# -------------------------------------------------------------------------------------------------

analysis.enable_eigenvalue_errors = True  # Compare eigenvalues from eigendecomposition
analysis.error_matrix_norms = "frobenius"  # Options: "frobenius", "spectral", ["frobenius", "spectral"]
analysis.error_state_inputs = "analysis/examples/initial_state.npy"  # String or list of .npy files

# Error types computed:
#   1. Eigenvalue errors: Compare all eigenvalues (requires eigendecomposition_matrices="both")
#   2. Matrix norm errors: Frobenius and/or spectral norms
#   3. State-dependent errors: ||H_exact|ψ⟩ - H_approx|ψ⟩|| for specific states

# -------------------------------------------------------------------------------------------------
# 7. EXACT NUMERICAL SIMULATION
# -------------------------------------------------------------------------------------------------

analysis.exact_simulation_inputs = "analysis/examples/initial_state.npy"

# This applies the exact Hamiltonian to input state(s)
# Output: input.npy → input_exact_final.npy

# Compare with approximate simulation:
#   approximate: input_final.npy
#   exact: input_exact_final.npy

# =================================================================================================
