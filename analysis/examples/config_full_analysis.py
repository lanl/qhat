"""
Comprehensive Configuration Example for QHAT Analysis

This configuration demonstrates ALL currently available analysis capabilities:
1. Resource estimation (quantum gate counts, qubit requirements)
2. Matrix output (save unitary matrix representation)
3. Numerical simulation (apply unitary to input state vectors)
4. Exact Hamiltonian matrix computation
5. Eigendecomposition analysis (exact and/or approximate)
6. Error analysis (eigenvalue errors, matrix norms, state-dependent errors)

Usage:
    python3.11 -m qhat.analysis.driver config_full_analysis.py
"""

# Energy error budget: split between Trotterization and phase estimation
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

# Optional: Set energy bounds if known (improves phase qubit estimation)
# hamiltonian.set_energy_lower_bound(-5.0, exact=False)
# hamiltonian.set_energy_upper_bound(0.0, exact=False)

# =================================================================================================
# UNITARY ENCODING
# =================================================================================================
# How to encode the Hamiltonian as a quantum circuit
# Currently supports: ramped_trotter, pauli_lcu, double_factorization

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
# Options: "pyLIQTR", "cirq"

analysis.resource_estimator = "pyLIQTR"

# Expected output in results:
#   - Clifford_count: Number of Clifford gates
#   - T_count: Number of T gates (dominant cost for fault-tolerant quantum computing)
#   - qubit_count: Number of logical qubits required
#   - circuit_depth: Depth of the quantum circuit

# -------------------------------------------------------------------------------------------------
# 2. UNITARY MATRIX OUTPUT
# -------------------------------------------------------------------------------------------------
# Generate and save the unitary matrix representation of the approximate algorithm
# This materializes the full 2^N × 2^N unitary matrix
# Practical for small systems (≤15 qubits, dimension ≤ 32768)

analysis.matrix_output_file = "Be-H_unitary_matrix.npz"

# Supported formats (auto-detected from file extension):
#   - .npz: NumPy compressed format (recommended, includes metadata)
#   - .h5 / .hdf5: HDF5 format (for compatibility with other tools)
#   - .txt: Human-readable text (only for very small matrices)

# The saved file includes:
#   - Unitary matrix U
#   - Metadata: timestamp, matrix shape, matrix norm
#   - Unitarity error: ||U†U - I||_F (should be ~1e-15 for correct implementations)

# -------------------------------------------------------------------------------------------------
# 3. EXACT HAMILTONIAN MATRIX OUTPUT
# -------------------------------------------------------------------------------------------------
# Compute and save the exact matrix representation of the Hamiltonian (no approximations)

analysis.exact_matrix_output_file = "Be-H_exact_hamiltonian.npz"

# This computes the true Hamiltonian matrix H directly from the Pauli string representation.
# Useful for:
#   - Validating approximate algorithms (compare exact vs approximate eigenvalues)
#   - Small-scale testing and verification
#   - Computing exact ground state energies

# Implementation notes:
#   - Uses memory threshold to choose between dense and sparse representations
#   - Default threshold: 16 GB (allows dense matrices up to ~15 qubits)
#   - Dense matrices: Full materialization, can be saved to file
#   - Sparse matrices (above threshold): Matrix-free operator, enables analysis but no file output

# Optional: Configure memory threshold for dense/sparse selection
# analysis.matrix_memory_threshold_gb = 16.0  # Default: 16 GB
# Examples:
#   - 1.0 GB allows dense up to ~13 qubits
#   - 16.0 GB allows dense up to ~15 qubits
#   - 64.0 GB allows dense up to ~16 qubits

# Supported formats (same as unitary matrix output):
#   - .npz: NumPy compressed format (recommended)
#   - .h5 / .hdf5: HDF5 format
#   - .txt: Human-readable text (only for very small matrices)

# The saved file includes:
#   - Hamiltonian matrix H
#   - Metadata: timestamp, matrix shape, matrix norm
#   - Hermiticity error: ||H - H†||_F (should be ~1e-15 for Hermitian matrices)

# -------------------------------------------------------------------------------------------------
# 4. NUMERICAL SIMULATION (APPROXIMATE)
# -------------------------------------------------------------------------------------------------
# Apply the approximate unitary matrix to one or more input quantum states
# Useful for testing algorithm behavior on specific initial states

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

# Output:
#   - For each input "state.npy", creates "state_final.npy"
#   - Reports final state norm (should be 1.0 for unitary evolution)
#   - Can be used to verify algorithm correctness or study state evolution

# -------------------------------------------------------------------------------------------------
# 5. EIGENDECOMPOSITION ANALYSIS
# -------------------------------------------------------------------------------------------------
# Compute eigenvalues and eigenvectors of exact Hamiltonian and/or approximate unitary

analysis.eigendecomposition_matrices = "both"  # Options: None, "exact", "approximate", "both"

# This analysis performs full eigendecomposition (all eigenvalues and eigenvectors):
#   - Exact matrix: Diagonalizes Hamiltonian H to find eigenenergies and eigenstates
#   - Approximate matrix: Diagonalizes unitary U and converts phases to energies
#   - Results are sorted by eigenenergy (ground state = index 0)
#   - Only feasible for small-medium systems (~15-18 qubits, dimension ≤ 256k)

# Configuration options:
#   - None (default): Eigendecomposition disabled
#   - "exact": Only eigendecompose exact Hamiltonian matrix
#   - "approximate": Only eigendecompose approximate unitary matrix
#   - "both": Eigendecompose both (required for eigenvalue error analysis)

# Output files:
#   - exact_eigendecomposition.npz (if "exact" or "both")
#       Contains: eigenenergies, eigenvectors, energy_shift applied
#   - approximate_eigendecomposition.npz (if "approximate" or "both")
#       Contains: eigenenergies (corrected for energy shift), eigenvectors, timestep used

# Important notes:
#   - Energy shift correction: The approximate time evolution uses a shifted Hamiltonian
#     H̃ = H - E_min·I to reduce phase qubit requirements. The eigendecomposition analysis
#     automatically corrects approximate eigenenergies to the same scale as exact eigenenergies
#     for proper comparison.
#   - Requires timestep for approximate case (automatically extracted from unitary config)

# -------------------------------------------------------------------------------------------------
# 6. ERROR ANALYSIS
# -------------------------------------------------------------------------------------------------
# Compare exact and approximate representations using multiple error metrics

analysis.enable_eigenvalue_errors = True  # Compare all eigenvalues from eigendecomposition
analysis.error_matrix_norms = "frobenius"  # Options: "frobenius", "spectral", ["frobenius", "spectral"]
analysis.error_state_inputs = "examples/Be-H_1.30_sto-6g_as-003-003_jw.npy"  # String or list of state files

# This analysis computes three independent error types:
#
# 1. EIGENVALUE ERRORS:
#    - Compares all eigenvalues computed in the eigendecomposition
#    - Requires eigendecomposition_matrices = "both"
#    - Reports absolute errors, relative errors, max error, RMS error
#    - Accounts for energy shift correction automatically
#
# 2. MATRIX NORM ERRORS:
#    - Computes ||H_exact - H_approx||_F (Frobenius norm)
#    - Computes ||H_exact - H_approx||_2 (spectral norm) if requested
#    - Frobenius: Fast, measures total element-wise difference
#    - Spectral: Slower, measures worst-case error (physically meaningful)
#
# 3. STATE-DEPENDENT ERRORS:
#    - Computes ||H_exact|ψ⟩ - H_approx|ψ⟩|| for specified states
#    - Fast: Just applies operators to states
#    - Most relevant for practical applications (how accurate for states you care about?)

# Configuration options:
#   enable_eigenvalue_errors:
#     - False (default): Eigenvalue errors disabled
#     - True: Compute errors for ALL eigenvalues from eigendecomposition
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
#     - Same format as numerical_simulation_inputs

# Output: error_analysis.npz containing all computed errors with metadata

# System size guidance:
#   - ≤15 qubits: All error types feasible
#   - 16-20 qubits: eigenvalue_errors + error_state_inputs recommended
#   - 20+ qubits: error_state_inputs only (best scaling)

# Examples for different use cases:
#
# Comprehensive validation (small systems):
# analysis.enable_eigenvalue_errors = True
# analysis.error_matrix_norms = ["frobenius", "spectral"]
# analysis.error_state_inputs = ["examples/ground.npy", "examples/excited.npy"]
#
# Standard validation (medium systems):
# analysis.enable_eigenvalue_errors = True
# analysis.error_matrix_norms = "frobenius"
# analysis.error_state_inputs = "examples/ground.npy"
#
# Minimal validation (large systems):
# analysis.enable_eigenvalue_errors = False  # Too expensive
# analysis.error_matrix_norms = None  # Too expensive
# analysis.error_state_inputs = "examples/ground.npy"  # Scales well

# =================================================================================================
# SUMMARY OF ALL ANALYSES ENABLED IN THIS CONFIG
# =================================================================================================
# This configuration enables every available analysis feature:
#   ✓ Resource estimation (pyLIQTR)
#   ✓ Unitary matrix output (Be-H_unitary_matrix.npz)
#   ✓ Exact Hamiltonian matrix output (Be-H_exact_hamiltonian.npz)
#   ✓ Numerical simulation with approximate unitary
#   ✓ Eigendecomposition of both exact and approximate matrices
#   ✓ Error analysis: eigenvalue errors, Frobenius norm, state-dependent errors
#
# This serves as a complete example of QHAT's analysis capabilities.
# For production runs, you may want to disable some analyses to improve performance.
# =================================================================================================
