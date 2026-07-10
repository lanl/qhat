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

# -------------------------------------------------------------------------------------------------
# 4. EXACT MATRIX COMPUTATION
# -------------------------------------------------------------------------------------------------
# Compute the exact matrix representation of the Hamiltonian without any approximations

analysis.exact_matrix_output_file = "Be-H_exact_hamiltonian.npz"

# This computes the true Hamiltonian matrix H directly from the Pauli string representation.
# Useful for:
#   - Validating approximate algorithms (compare exact vs approximate eigenvalues)
#   - Small-scale testing and verification
#   - Computing exact ground state energies

# Implementation notes:
#   - Uses memory-based threshold to choose between dense and sparse representations
#   - Default threshold: 16 GB (allows dense matrices up to ~15 qubits)
#   - Dense matrices: Full materialization, can be saved to file
#   - Sparse matrices (above threshold): Matrix-free operator, works with scipy eigensolvers
#   - Matrix-free operators cannot be directly saved to file

# Optional: Configure memory threshold for dense/sparse selection
# analysis.matrix_memory_threshold_gb = 16.0  # Default: 16 GB
# Examples:
#   - 1.0 GB threshold allows dense up to ~13 qubits
#   - 16.0 GB threshold allows dense up to ~15 qubits
#   - 64.0 GB threshold allows dense up to ~16 qubits

# Supported formats (auto-detected from extension):
#   - .npz: NumPy compressed format (recommended, includes metadata)
#   - .h5 / .hdf5: HDF5 format
#   - .txt: Human-readable text (only for very small matrices)

# Note: For systems exceeding the memory threshold, the exact matrix computation will use
#       a matrix-free operator (enabling other analyses) but skip the file output.

# -------------------------------------------------------------------------------------------------
# 5. EIGENDECOMPOSITION ANALYSIS
# -------------------------------------------------------------------------------------------------
# Compute eigenvalues and eigenvectors of exact and/or approximate matrices.  May not be practical
# or even possible for large systems.

# Eigendecomposition: always computes full spectrum, sorted by eigenenergy
analysis.eigendecomposition_matrices = "both"  # Options: None (disabled), "exact", "approximate", "both"

# This analysis computes full-spectrum eigendecompositions (all eigenstates, sorted by energy):
#   - Always computes ALL eigenvalues and eigenvectors
#   - Results are automatically sorted by eigenenergy (ground state = index 0)
#   - Only feasible for small systems

# Configuration option:
#   eigendecomposition_matrices:
#     - None (default): Eigendecomposition disabled
#     - "exact": Only eigendecompose exact Hamiltonian matrix
#     - "approximate": Only eigendecompose approximate unitary matrix
#     - "both": Eigendecompose both (required for eigenvalue error analysis)

# Output files:
#   - exact_eigendecomposition.npz (if "exact" or "both")
#   - approximate_eigendecomposition.npz (if "approximate" or "both")

# Examples:

# Full spectrum for exact Hamiltonian only
# analysis.eigendecomposition_matrices = "exact"

# Full spectrum for approximate unitary only
# analysis.eigendecomposition_matrices = "approximate"

# Both (required for comparing exact vs approximate eigenenergies)
# analysis.eigendecomposition_matrices = "both"

# =================================================================================================
# 5. EIGENDECOMPOSITION ANALYSIS
# =================================================================================================
# Compute eigenvalues and eigenvectors of exact and/or approximate matrices

# Eigendecomposition: always computes full spectrum, sorted by eigenenergy
analysis.eigendecomposition_matrices = "both"  # Options: None (disabled), "exact", "approximate", "both"

# This analysis computes full-spectrum eigendecompositions (all eigenstates, sorted by energy):
#   - Always computes ALL eigenvalues and eigenvectors
#   - Results are automatically sorted by eigenenergy (ground state = index 0)
#   - Only feasible for small systems

# Configuration options:
#   num_eigenvalues:
#     - 0 (default): Eigendecomposition disabled
#     - Positive int (e.g., 5): Compute that many eigenvalues (recommended for large systems)
#     - "all": Full eigendecomposition (only for systems ≤10 qubits)
#
#   eigendecomposition_matrices:
#     - "approximate" (default): Only eigendecompose unitary matrix
#     - "exact": Only eigendecompose exact Hamiltonian matrix
#     - "both": Eigendecompose both (useful for comparison)
#
#   which_eigenvalues (ignored for "all"):
#     - "smallest" (default): Algebraically smallest (most negative, ground state)
#     - "largest": Algebraically largest (most positive)
#     - "both": Compute smallest AND largest (returns 2k eigenvalues)

# Output files:
#   - exact_eigendecomposition.npz (if "exact" or "both")
#   - approximate_eigendecomposition.npz (if "approximate" or "both")

# Examples:

# Ground state energy comparison (most common use case for large systems)
# analysis.num_eigenvalues = 1
# analysis.eigendecomposition_matrices = "both"
# analysis.which_eigenvalues = "smallest"

# Low-lying spectrum (5 lowest-energy states)
# analysis.num_eigenvalues = 5
# analysis.eigendecomposition_matrices = "exact"
# analysis.which_eigenvalues = "smallest"

# High and low energy states (3 smallest + 3 largest = 6 total)
# analysis.num_eigenvalues = 3
# analysis.eigendecomposition_matrices = "both"
# analysis.which_eigenvalues = "both"

# Full spectrum for small system
# analysis.num_eigenvalues = "all"
# analysis.eigendecomposition_matrices = "exact"

# =================================================================================================
# 6. ERROR ANALYSIS
# =================================================================================================
# Compare exact and approximate representations using three types of error metrics

analysis.enable_eigenvalue_errors = True  # Compare all eigenvalues from eigendecomposition
analysis.error_matrix_norms = "frobenius"  # Options: "frobenius", "spectral", or ["frobenius", "spectral"]
analysis.error_state_inputs = "examples/initial_state.npy"  # Can be string or list

# This analysis computes three independent error types:
#   1. Eigenvalue errors: Compare all eigenvalues computed in the eigendecomposition
#   2. Matrix norm errors: Frobenius and/or spectral norms of difference matrix
#   3. State-dependent errors: ||H_exact|ψ⟩ - H_approx|ψ⟩|| for specific states

# Configuration options:
#   enable_eigenvalue_errors:
#     - False (default): Eigenvalue errors disabled
#     - True: Compute errors for ALL eigenvalues from eigendecomposition
#     - Requires eigendecompositions (set eigendecomposition_matrices = "both")
#     - Compares ALL eigenvalues (full spectrum) element-wise after sorting
#
#   error_matrix_norms:
#     - None (default): Matrix norm errors disabled
#     - "frobenius": Fast element-wise difference norm
#     - "spectral": Physically meaningful worst-case norm (slower)
#     - ["frobenius", "spectral"]: Both norms
#
#   error_state_inputs:
#     - None (default): State-dependent errors disabled
#     - String or list of .npy files containing quantum states
#     - Fast: just applies operators to states

# Output: error_analysis.npz containing all computed errors

# Examples:

# Eigenvalue error comparison (compare all eigenvalues from eigendecomposition)
# analysis.enable_eigenvalue_errors = True

# Matrix norms (physical bounds, but expensive for large systems)
# analysis.error_matrix_norms = "frobenius"  # Fast
# analysis.error_matrix_norms = ["frobenius", "spectral"]  # Both (slower)

# State-dependent errors (best scaling, physically relevant)
# analysis.error_state_inputs = "examples/initial_state.npy"
# analysis.error_state_inputs = ["examples/ground.npy", "examples/excited.npy"]

# All three error types (comprehensive validation for small-medium systems)
# analysis.enable_eigenvalue_errors = True
# analysis.error_matrix_norms = "frobenius"
# analysis.error_state_inputs = ["examples/ground.npy", "examples/excited.npy"]

# Minimal setup for large systems
# analysis.enable_eigenvalue_errors = True  # Compare eigenvalues
# analysis.error_state_inputs = "examples/ground.npy"  # Scales well
# # Skip error_matrix_norms for large systems (too slow)

# =================================================================================================
# 5. EIGENDECOMPOSITION ANALYSIS (Branch 2: Now Available!)
# =================================================================================================
# Compute eigenvalues and eigenvectors of exact and/or approximate matrices

# Eigendecomposition: always computes full spectrum, sorted by eigenenergy
analysis.eigendecomposition_matrices = "both"  # Options: None (disabled), "exact", "approximate", "both"

# This analysis computes full-spectrum eigendecompositions (all eigenstates, sorted by energy):
#   - Always computes ALL eigenvalues and eigenvectors
#   - Results are automatically sorted by eigenenergy (ground state = index 0)
#   - Only feasible for small systems (~15-18 qubits, dimension ≤256k)

# Configuration options:
#   num_eigenvalues:
#     - 0 (default): Eigendecomposition disabled
#     - Positive int (e.g., 5): Compute that many eigenvalues (recommended for large systems)
#     - "all": Full eigendecomposition (only for systems ≤10 qubits)
#
#   eigendecomposition_matrices:
#     - "approximate" (default): Only eigendecompose unitary matrix
#     - "exact": Only eigendecompose exact Hamiltonian matrix
#     - "both": Eigendecompose both (useful for comparison)
#
#   which_eigenvalues (ignored for "all"):
#     - "smallest" (default): Algebraically smallest (most negative, ground state)
#     - "largest": Algebraically largest (most positive)
#     - "both": Compute smallest AND largest (returns 2k eigenvalues)
#
# System size guidance:
#   - ≤10 qubits: "all" feasible
#   - 10-12 qubits: "all" possible but expensive
#   - ≥12 qubits: Use partial (specify k=5-10)
#   - 20+ qubits: MUST use partial
#   - 25-30 qubits: Partial still scales well

# Output files:
#   - exact_eigendecomposition.npz (if "exact" or "both")
#   - approximate_eigendecomposition.npz (if "approximate" or "both")

# Examples:

# Ground state energy comparison (most common use case for large systems)
# analysis.num_eigenvalues = 1
# analysis.eigendecomposition_matrices = "both"
# analysis.which_eigenvalues = "smallest"

# Low-lying spectrum (5 lowest-energy states)
# analysis.num_eigenvalues = 5
# analysis.eigendecomposition_matrices = "exact"
# analysis.which_eigenvalues = "smallest"

# High and low energy states (3 smallest + 3 largest = 6 total)
# analysis.num_eigenvalues = 3
# analysis.eigendecomposition_matrices = "both"
# analysis.which_eigenvalues = "both"

# Full spectrum for small system
# analysis.num_eigenvalues = "all"
# analysis.eigendecomposition_matrices = "exact"

# =================================================================================================
# FUTURE ANALYSES (Coming in subsequent branches)
# =================================================================================================
# Apply the exact Hamiltonian matrix (no approximations) to input state(s)

# Branch 3: Error Metrics
# Eigenvalue errors: compares all eigenstates (full spectrum, sorted by energy)
analysis.error_matrix_norms = ["frobenius", "spectral"]  # Matrix difference norms
analysis.error_state_inputs = ["examples/ground_state.npy"]  # State-dependent errors

# Branch 4: Exact Numerical Simulation
# analysis.exact_simulation_inputs = "examples/initial_state.npy"

# This applies the exact Hamiltonian H to the input state(s), producing exact time evolution.
# Useful for:
#   - Comparing exact evolution with approximate algorithm results
#   - Validating algorithm correctness on small-scale problems
#   - Understanding exact quantum dynamics

# Configuration:
#   exact_simulation_inputs:
#     - None (default): Exact simulation disabled
#     - String: Single input state file
#     - List: Multiple input state files
#
#   Output files:
#     - input.npy → input_exact_final.npy (note: _exact_final suffix)
#     - Distinguishes from approximate simulation (_final suffix)

# System size guidance:
#   - ≤15 qubits: Fast (uses dense exact matrix)
#   - 16-30 qubits: Uses matrix-free operator (scales well)
#   - Limited by state vector size (16 GB at 30 qubits)

# Examples:

# Single state exact simulation
# analysis.exact_simulation_inputs = "examples/initial_state.npy"

# Multiple states
# analysis.exact_simulation_inputs = [
#     "examples/ground_state.npy",
#     "examples/excited_state.npy",
#     "examples/superposition_state.npy"
# ]

# Compare exact vs approximate evolution (run both!)
# analysis.numerical_simulation_inputs = "examples/initial_state.npy"  # Approximate
# analysis.exact_simulation_inputs = "examples/initial_state.npy"      # Exact
# # Produces both:
# #   initial_state_final.npy        (approximate result)
# #   initial_state_exact_final.npy  (exact result)
# # Then you can manually compare the output states

# Full validation workflow (all analyses together)
# This demonstrates using all available QHAT analysis capabilities:
# analysis.resource_estimator = "pyLIQTR"                  # Gate counts
# analysis.matrix_output_file = "unitary.npz"              # Approximate matrix
# analysis.exact_matrix_output_file = "exact.npz"          # Exact matrix
# analysis.num_eigenvalues = 5                              # Eigendecomposition
# analysis.eigendecomposition_matrices = "both"             # Both exact and approximate
# analysis.which_eigenvalues = "smallest"                   # Ground + excited states
# analysis.error_num_eigenvalues = 1                        # Ground state error
# analysis.error_state_inputs = "examples/ground.npy"       # State-dependent error
# analysis.numerical_simulation_inputs = "examples/initial_state.npy"   # Approximate evolution
# analysis.exact_simulation_inputs = "examples/initial_state.npy"       # Exact evolution

# =================================================================================================
