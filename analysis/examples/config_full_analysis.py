"""
Comprehensive Configuration Example for QHAT Analysis

This configuration demonstrates ALL currently available analysis capabilities:
1. Resource estimation (quantum gate counts, qubit requirements)
2. Flexible operator output (NEW: save any combination of operator forms)
3. Numerical simulation (apply unitary to input state vectors)
4. Error analysis (eigenvalue errors, matrix norms, state-dependent errors)

Usage:
    python3.11 -m qhat.analysis.driver examples/config_full_analysis.py
"""

# =================================================================================================
# GENERAL CONFIGURATION
# =================================================================================================

general.print_verbose()  # Options: print_default(), print_verbose(), print_debug()
general.logfile = "Be-H_full_analysis.log"

# =================================================================================================
# HAMILTONIAN SPECIFICATION
# =================================================================================================
# Load a second-quantized Hamiltonian from file
# Supports: .npz (NumPy), .h5/.hdf5 (HDF5)

hamiltonian.load_second_quantization("examples/Be-H_1.30_sto-6g_as-003-003.tensors.npz")

# Optional: Set energy bounds if known (improves phase qubit estimation)
# hamiltonian.set_energy_lower_bound(-5.0, exact=False)
# hamiltonian.set_energy_upper_bound(0.0, exact=False)

# =================================================================================================
# UNITARY ENCODING
# =================================================================================================
# How to encode the Hamiltonian as a quantum circuit
# Currently supports: ramped_trotter, pauli_lcu, double_factorization

# Define energy error budget (used for both Trotter and phase estimation)
energy_error = 0.001  # 1 meV in natural units

unitary.encode_ramped_trotter(
    energy_error=0.5 * energy_error,  # Split error budget between Trotter and phase estimation
    trotter_implementation="flattened",  # Options: "flattened", "recursive"
    trotter_combine_terms=True,  # Combine commuting terms for efficiency
    ordering_method="lexicographical"  # Options: "lexicographical", "random", None
)

# Alternative encoding methods:
# unitary.encode_pauli_lcu(energy_error=0.5 * energy_error)
# unitary.encode_double_factorization(energy_error=0.5 * energy_error)

# =================================================================================================
# ALGORITHM SELECTION
# =================================================================================================
# The quantum algorithm to use
# Options: "time evolution", "qualtran textbook"

algorithm.method = "time evolution"

# For QPE-based algorithms, you can also specify:
# algorithm.num_phase_qubits = 10
# algorithm.probability_of_failure = 0.01

# =================================================================================================
# ANALYSIS CONFIGURATION
# =================================================================================================
# Request various analyses of the quantum algorithm

# -------------------------------------------------------------------------------------------------
# 1. RESOURCE ESTIMATION
# -------------------------------------------------------------------------------------------------
# Estimate quantum computing resources required (gate counts, qubit counts, circuit depth)

analysis.resource_estimator = "pyLIQTR"

# Expected output in results:
#   - Clifford_count: Number of Clifford gates
#   - T_count: Number of T gates (dominant cost for fault-tolerant quantum computing)
#   - qubit_count: Number of logical qubits required
#   - circuit_depth: Depth of the quantum circuit

# -------------------------------------------------------------------------------------------------
# 2. FLEXIBLE OPERATOR OUTPUT (NEW API)
# -------------------------------------------------------------------------------------------------
# Save any combination of operator forms using the new flexible API
#
# Each operator has 4 independent choices:
#   operator: 'exact' or 'approximate'
#   form: 'hamiltonian' or 'time_evolution'
#   shift: 'unshifted' or 'shifted'
#   representation: matrix or eigendecomposition
#
# This gives 2 × 2 × 2 × 2 = 16 possible output files per operator type (matrix vs eigendecomp)

# ALL PARAMETERS ARE REQUIRED - be explicit about what you want!

# =================================================================================================
# 2A. MATRIX OUTPUTS
# =================================================================================================
# Call save_matrix_to_file() for each matrix you want to save

# --- EXACT OPERATOR: All 4 matrix forms ---

# 1. Exact Hamiltonian, unshifted (physical energy scale)
analysis.save_matrix_to_file(
    filename='H_exact_physical.npz',
    operator='exact',
    form='hamiltonian',
    shift='unshifted'
)

# 2. Exact Hamiltonian, shifted (QPE energy scale)
analysis.save_matrix_to_file(
    filename='H_exact_QPE.npz',
    operator='exact',
    form='hamiltonian',
    shift='shifted'
)

# 3. Exact time-evolution, unshifted (from physical H)
analysis.save_matrix_to_file(
    filename='U_exact_physical.npz',
    operator='exact',
    form='time_evolution',
    shift='unshifted'
)

# 4. Exact time-evolution, shifted (from QPE H)
analysis.save_matrix_to_file(
    filename='U_exact_QPE.npz',
    operator='exact',
    form='time_evolution',
    shift='shifted'
)

# --- APPROXIMATE OPERATOR: All 4 matrix forms ---

# 5. Approximate Hamiltonian, unshifted (physical scale, extracted from Trotter/LCU)
analysis.save_matrix_to_file(
    filename='H_approx_physical.npz',
    operator='approximate',
    form='hamiltonian',
    shift='unshifted'
)

# 6. Approximate Hamiltonian, shifted (QPE scale, extracted from Trotter/LCU)
analysis.save_matrix_to_file(
    filename='H_approx_QPE.npz',
    operator='approximate',
    form='hamiltonian',
    shift='shifted'
)

# 7. Approximate time-evolution, unshifted (Trotter/LCU on physical scale)
analysis.save_matrix_to_file(
    filename='U_approx_physical.npz',
    operator='approximate',
    form='time_evolution',
    shift='unshifted'
)

# 8. Approximate time-evolution, shifted (Trotter/LCU as computed)
analysis.save_matrix_to_file(
    filename='U_approx_QPE.npz',
    operator='approximate',
    form='time_evolution',
    shift='shifted'
)

# =================================================================================================
# 2B. EIGENDECOMPOSITION OUTPUTS
# =================================================================================================
# Call save_eigendecomposition_to_file() for each eigendecomposition you want to save

# --- EXACT OPERATOR: All 4 eigendecomposition forms ---

# 1. Exact Hamiltonian eigendecomp, unshifted (physical eigenenergies)
analysis.save_eigendecomposition_to_file(
    filename='H_exact_physical_eig.npz',
    operator='exact',
    form='hamiltonian',
    shift='unshifted'
)

# 2. Exact Hamiltonian eigendecomp, shifted (QPE eigenenergies)
analysis.save_eigendecomposition_to_file(
    filename='H_exact_QPE_eig.npz',
    operator='exact',
    form='hamiltonian',
    shift='shifted'
)

# 3. Exact time-evolution eigendecomp, unshifted (phases from physical H)
analysis.save_eigendecomposition_to_file(
    filename='U_exact_physical_eig.npz',
    operator='exact',
    form='time_evolution',
    shift='unshifted'
)

# 4. Exact time-evolution eigendecomp, shifted (phases from QPE H)
analysis.save_eigendecomposition_to_file(
    filename='U_exact_QPE_eig.npz',
    operator='exact',
    form='time_evolution',
    shift='shifted'
)

# --- APPROXIMATE OPERATOR: All 4 eigendecomposition forms ---

# 5. Approximate Hamiltonian eigendecomp, unshifted (physical scale)
analysis.save_eigendecomposition_to_file(
    filename='H_approx_physical_eig.npz',
    operator='approximate',
    form='hamiltonian',
    shift='unshifted'
)

# 6. Approximate Hamiltonian eigendecomp, shifted (QPE scale)
analysis.save_eigendecomposition_to_file(
    filename='H_approx_QPE_eig.npz',
    operator='approximate',
    form='hamiltonian',
    shift='shifted'
)

# 7. Approximate time-evolution eigendecomp, unshifted (phases, physical scale)
analysis.save_eigendecomposition_to_file(
    filename='U_approx_physical_eig.npz',
    operator='approximate',
    form='time_evolution',
    shift='unshifted'
)

# 8. Approximate time-evolution eigendecomp, shifted (phases, QPE scale)
analysis.save_eigendecomposition_to_file(
    filename='U_approx_QPE_eig.npz',
    operator='approximate',
    form='time_evolution',
    shift='shifted'
)

# =================================================================================================
# FLEXIBLE OUTPUT API: Understanding the Parameters
# =================================================================================================
#
# operator: Which quantum operator?
#   'exact' = true Hamiltonian from Pauli strings (no approximations)
#   'approximate' = algorithm output (Trotter, LCU, double-factorization, etc.)
#
# form: How to represent the operator?
#   'hamiltonian' = Hamiltonian matrix H (eigenvectors are energy eigenstates)
#   'time_evolution' = Time-evolution operator U = exp(-i*H*t) (eigenvectors same, eigenvalues on unit circle)
#
# shift: Which energy scale?
#   'unshifted' = physical energy scale (eigenvalues can be negative)
#   'shifted' = QPE energy scale (H' = H + E*I, all eigenvalues positive)
#     - The energy shift E = |min(eigenvalues)| is applied automatically
#     - QPE requires positive eigenvalues to extract phases correctly
#     - For error analysis, usually compare on unshifted (physical) scale
#
# representation: (determined by method called)
#   save_matrix_to_file() = dense matrix (2^n × 2^n array)
#   save_eigendecomposition_to_file() = eigenvalues + eigenvectors
#
# Key insights:
#   - All 16 forms are mathematically equivalent (just different views)
#   - OperatorRepresentation handles all conversions automatically
#   - Most common for error analysis: 'unshifted' (physical scale) comparisons
#   - QPE algorithm works on 'shifted' scale but we typically report physical energies

# =================================================================================================
# COMMON USE CASES (comment out the above and use these instead)
# =================================================================================================

# --- Use Case 1: Basic error analysis (most common) ---
# Compare exact vs approximate on physical scale
#
# analysis.save_matrix_to_file(
#     filename='U_exact.npz',
#     operator='exact',
#     form='time_evolution',
#     shift='unshifted'
# )
# analysis.save_matrix_to_file(
#     filename='U_approx.npz',
#     operator='approximate',
#     form='time_evolution',
#     shift='unshifted'
# )
# analysis.save_eigendecomposition_to_file(
#     filename='H_exact_eig.npz',
#     operator='exact',
#     form='hamiltonian',
#     shift='unshifted'
# )
# analysis.save_eigendecomposition_to_file(
#     filename='H_approx_eig.npz',
#     operator='approximate',
#     form='hamiltonian',
#     shift='unshifted'
# )

# --- Use Case 2: Debug energy shift corrections ---
# Verify shift/unshift operations are working correctly
#
# analysis.save_eigendecomposition_to_file(
#     filename='H_exact_physical_eig.npz',
#     operator='exact',
#     form='hamiltonian',
#     shift='unshifted'
# )
# analysis.save_eigendecomposition_to_file(
#     filename='H_exact_QPE_eig.npz',
#     operator='exact',
#     form='hamiltonian',
#     shift='shifted'
# )
# # Eigenvalues should differ by exactly E (the energy shift)

# --- Use Case 3: Verify H ↔ U conversion ---
# Same operator in both representations should have matching eigenvalues
#
# analysis.save_eigendecomposition_to_file(
#     filename='exact_as_H_eig.npz',
#     operator='exact',
#     form='hamiltonian',
#     shift='unshifted'
# )
# analysis.save_eigendecomposition_to_file(
#     filename='exact_as_U_eig.npz',
#     operator='exact',
#     form='time_evolution',
#     shift='unshifted'
# )
# # After converting U eigenvalues (phases) to energies, should match H eigenvalues

# -------------------------------------------------------------------------------------------------
# 3. NUMERICAL SIMULATION (APPROXIMATE UNITARY)
# -------------------------------------------------------------------------------------------------
# Apply the approximate unitary matrix to one or more input quantum states

# Single input state:
analysis.numerical_simulation_inputs = "examples/Be-H_1.30_sto-6g_as-003-003_jw.npy"

# Multiple input states (processes each independently):
# analysis.numerical_simulation_inputs = [
#     "examples/initial_state.npy",
#     "examples/ground_state.npy",
#     "examples/excited_state.npy"
# ]

# Input states must be:
#   - 1D complex arrays
#   - Dimension 2^N (matching system size)
#   - Saved as NumPy .npy files
#   - Normalized (||ψ|| = 1) for physical quantum states
#
# Output:
#   - For each input "state.npy", creates "state_final.npy"
#   - Reports final state norm (should be 1.0 for unitary evolution)

# Optional: Apply exact Hamiltonian evolution to states (for comparison)
# analysis.exact_simulation_inputs = "examples/Be-H_1.30_sto-6g_as-003-003_jw.npy"
# Output: For each input "state.npy", creates "state_exact_final.npy"

# -------------------------------------------------------------------------------------------------
# 4. ERROR ANALYSIS (OLD API - still supported)
# -------------------------------------------------------------------------------------------------
# Compare exact and approximate representations using multiple error metrics

# Note: With the new flexible output API, you might not need the old error analysis API
# since you can directly save and compare any operator forms you want. However, the
# error analysis API still provides convenient aggregate error metrics.

analysis.enable_eigenvalue_errors = True  # Compare all eigenvalues
analysis.error_matrix_norms = ["frobenius", "spectral"]  # Norm-based errors
analysis.error_state_inputs = "examples/Be-H_1.30_sto-6g_as-003-003_jw.npy"  # State-dependent errors

# This computes three independent error types:
#
# 1. EIGENVALUE ERRORS:
#    - Compares eigenenergies: λ(H_exact) vs λ(H_approx)
#    - Reports absolute errors, relative errors, max error, RMS error
#
# 2. MATRIX NORM ERRORS:
#    - Compares time-evolution operators: ||U_exact - U_approx||
#    - Frobenius norm (fast) and/or spectral norm (physically meaningful)
#
# 3. STATE-DEPENDENT ERRORS:
#    - Compares time-evolved states: ||U_exact|ψ⟩ - U_approx|ψ⟩||
#    - Most relevant for practical applications

# Output: error_analysis.npz containing all computed errors with metadata

# =================================================================================================
# SUMMARY OF ALL OPTIONS IN THIS CONFIG
# =================================================================================================
# This configuration demonstrates EVERY available analysis feature:
#
# ✓ Resource estimation (pyLIQTR)
#
# ✓ Flexible operator outputs (NEW API):
#     8 matrix files covering all exact operator forms
#     8 matrix files covering all approximate operator forms
#     8 eigendecomposition files covering all exact operator forms
#     8 eigendecomposition files covering all approximate operator forms
#     = 32 output files total (all 16 possible combinations × 2 representations)
#
# ✓ Numerical simulation with approximate unitary
#
# ✓ Error analysis (eigenvalue, matrix norm, state-dependent)
#
# For production runs, you would typically use only a small subset:
#   - Basic validation: 2-4 files (U_exact, U_approx, maybe eigendecompositions)
#   - Debugging: 4-8 files (compare shifted vs unshifted, H vs U forms)
#   - Comprehensive: This full config (all 32 files for complete analysis)
#
# The flexible API gives you complete control over exactly what gets saved,
# ensuring you get the precise operator forms needed for your analysis workflow.
# =================================================================================================
