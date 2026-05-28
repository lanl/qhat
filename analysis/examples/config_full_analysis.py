"""
Comprehensive Configuration Example for QHAT Analysis

This configuration demonstrates ALL currently available analysis capabilities:
1. Resource estimation (quantum gate counts, qubit requirements)
2. Matrix output (save unitary matrix representation)
3. Numerical simulation (apply unitary to input state vectors)

This serves as a baseline example showing the full power of the analysis framework
before adding exact matrix computation, eigendecomposition, and error analysis features.

Usage:
    python3.11 -m qhat.analysis.driver config_full_analysis.py
"""

# This configuration file assumes equipartition of energy between Trotterization error and phase
# estimation error (see usage of energy_error below)
energy_error = meV_to_Hartree(1e4)  # 0.01 keV

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

# Optional: Set energy bounds if known
# hamiltonian.set_energy_lower_bound(-5.0, exact=False)
# hamiltonian.set_energy_upper_bound(0.0, exact=False)

# =================================================================================================
# UNITARY ENCODING
# =================================================================================================
# How to encode the Hamiltonian as a quantum circuit
# Currently supports: ramped_trotter, double_factorization

unitary.encode_ramped_trotter(
    energy_error=0.5 * energy_error,  # Split error budget between Trotter and phase estimation
    trotter_implementation="flattened",  # Options: "flattened", "recursive"
    trotter_combine_terms=True,  # Combine commuting terms for efficiency
    ordering_method="lexicographical"  # Options: "lexicographical", "random", None
)

# Alternative: Double factorization encoding (for molecular Hamiltonians)
# unitary.encode_double_factorization(
#     energy_error=0.5 * energy_error
# )

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
# Options: "pyLIQTR", "cirq"

analysis.resource_estimator = "pyLIQTR"

# Expected output:
#   - Clifford_count: Number of Clifford gates
#   - T_count: Number of T gates (dominant cost for fault-tolerant quantum computing)
#   - qubit_count: Number of logical qubits required

# -------------------------------------------------------------------------------------------------
# 2. MATRIX OUTPUT
# -------------------------------------------------------------------------------------------------
# Generate and save the unitary matrix representation of the algorithm
# This materializes the full 2^N × 2^N unitary matrix
# Practical for small systems (≤15 qubits, dimension ≤ 32768)

analysis.matrix_output_file = "Be-H_unitary_matrix.npz"

# Supported formats (auto-detected from file extension):
#   - .npz: NumPy compressed format (recommended, includes metadata)
#   - .h5 / .hdf5: HDF5 format (for compatibility with other tools)
#   - .txt: Human-readable text (only for very small matrices)

# The saved file includes:
#   - Unitary matrix U
#   - Metadata: git hash, timestamp, unitarity check ||U†U - I||_F, matrix norm

# Alternative output files for different formats:
# analysis.matrix_output_file = "Be-H_unitary_matrix.h5"    # HDF5 format
# analysis.matrix_output_file = "Be-H_unitary_matrix.txt"   # Text format (small matrices only)

# -------------------------------------------------------------------------------------------------
# 3. NUMERICAL SIMULATION
# -------------------------------------------------------------------------------------------------
# Apply the unitary matrix to one or more input quantum states
# Useful for testing algorithm behavior on specific initial states

# Single input state:
analysis.numerical_simulation_inputs = "examples/initial_state.npy"

# Multiple input states (processes each independently):
# analysis.numerical_simulation_inputs = [
#     "examples/ground_state.npy",
#     "examples/excited_state_1.npy",
#     "examples/excited_state_2.npy",
#     "examples/superposition_state.npy"
# ]

# Input states must be:
#   - 1D complex arrays
#   - Dimension 2^N (power of 2)
#   - Saved as NumPy .npy files
#   - Normalized (||ψ|| = 1) for physical quantum states

# Output:
#   - For each input "state.npy", creates "state_final.npy"
#   - Reports final state norm (should be 1.0 for unitary evolution)
#   - Can be used to verify algorithm correctness or study state evolution

# To create input state files:
#   import numpy as np
#   # Example: 2-qubit ground state |00⟩
#   state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
#   np.save("ground_state.npy", state)

# =================================================================================================
# 4. EXACT MATRIX COMPUTATION (Branch 1: Now Available!)
# =================================================================================================
# Compute the exact matrix representation of the Hamiltonian without any approximations
# (no Trotter, no double-factorization)

# analysis.exact_matrix_output_file = "Be-H_exact_hamiltonian.npz"

# This computes the true Hamiltonian matrix H directly from the Pauli string representation.
# Useful for:
#   - Validating approximate algorithms (compare exact vs approximate eigenvalues)
#   - Small-scale testing and verification
#   - Computing exact ground state energies

# Implementation notes:
#   - For small systems (≤15 qubits): Materializes full dense matrix
#   - For large systems (>15 qubits): Creates matrix-free operator (LinearOperator)
#   - Matrix-free operators work with scipy eigensolvers but cannot be directly saved

# Supported formats (auto-detected from extension):
#   - .npz: NumPy compressed format (recommended, includes metadata)
#   - .h5 / .hdf5: HDF5 format
#   - .txt: Human-readable text (only for very small matrices)

# Note: For systems >15 qubits, the exact matrix is too large to save. The analysis
#       will create a matrix-free operator but skip the file output.

# =================================================================================================
# FUTURE ANALYSES (Coming in subsequent branches)
# =================================================================================================

# Branch 2: Eigendecomposition
# analysis.num_eigenvalues = 5  # Compute 5 lowest eigenvalues
# analysis.eigendecomposition_matrices = "both"  # Analyze both exact and approximate
# analysis.which_eigenvalues = "smallest"  # Ground state and low-lying excited states

# Branch 3: Error Metrics
# analysis.error_num_eigenvalues = 1  # Compare ground state energy
# analysis.error_matrix_norms = ["frobenius", "spectral"]  # Matrix difference norms
# analysis.error_state_inputs = ["examples/ground_state.npy"]  # State-dependent errors

# Branch 4: Exact Numerical Simulation
# analysis.exact_simulation_inputs = "examples/initial_state.npy"
# Simulate using exact Hamiltonian evolution (for comparison with approximate)

# =================================================================================================
