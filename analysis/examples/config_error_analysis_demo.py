"""
QHAT Error Analysis Tutorial

This configuration demonstrates how to use error analysis to quantify the
accuracy of approximate quantum algorithms.

Error analysis answers the question: "How good is my approximation?"
It provides three independent measures:
1. Eigenvalue errors - How accurate are the computed energies?
2. Matrix norm errors - How close is the approximate unitary to the exact one?
3. State-dependent errors - How accurately is my specific input state evolved?

Usage:
    python3.11 -m qhat.analysis.driver examples/config_error_analysis_demo.py
"""

energy_error = meV_to_Hartree(1e4)  # 0.01 keV

# =================================================================================================
# GENERAL CONFIGURATION
# =================================================================================================

general.print_verbose()  # Print detailed progress information
general.logfile = "error_analysis_demo.log"

# =================================================================================================
# HAMILTONIAN SPECIFICATION
# =================================================================================================

hamiltonian.load_second_quantization("examples/Be-H_1.30_sto-6g_as-003-003.tensors.npz")

# =================================================================================================
# UNITARY ENCODING - Trotter Decomposition
# =================================================================================================
# We use Trotter decomposition which introduces approximation errors.
# Error analysis will quantify how good this approximation is.

unitary.encode_ramped_trotter(
    energy_error=energy_error,
    trotter_implementation="flattened",
    trotter_combine_terms=True,
    ordering_method="lexicographical"
)

# =================================================================================================
# ALGORITHM - Time Evolution Only
# =================================================================================================
# Use time evolution (not QPE) so we can directly compare U_exact vs U_approx
# This is the cleanest setup for demonstrating error analysis

algorithm.method = "time evolution"

# =================================================================================================
# ANALYSIS CONFIGURATION
# =================================================================================================

# The "width" of an algorithm (number of qubits) and "length" of an algorithm (a proxy for the
# runtime).
analysis.resource_estimator = "pyLIQTR"

# -------------------------------------------------------------------------------------------------
# Matrix and Eigendecomposition Output - Required for Error Analysis
# -------------------------------------------------------------------------------------------------
# Error analysis needs matrices and eigendecompositions from both exact and approximate operators

# EXACT OPERATOR: Use flexible API to save exact Hamiltonian
analysis.save_matrix_to_file(
    filename='error_demo_exact_hamiltonian.npz',
    operator='exact',
    form='hamiltonian',
    shift='unshifted'
)
# The exact Hamiltonian (no approximations)
# Will be automatically converted to U_exact = exp(-i*H*t) for error comparison

# EXACT EIGENDECOMPOSITION: Save exact eigenvalues for eigenvalue error analysis
analysis.save_eigendecomposition_to_file(
    filename='error_demo_exact_eigendecomp.npz',
    operator='exact',
    form='hamiltonian',
    shift='unshifted'
)
# Eigenenergies and eigenvectors of H_exact (ground state + excited states)

# APPROXIMATE ALGORITHM: Save the full algorithm circuit matrix
analysis.algorithm_matrix_output_file = "error_demo_algorithm_output.npz"
# The unitary matrix produced by the Trotter circuit
# For time evolution: this is U_approx
# For QPE: this would be the full QPE circuit

# APPROXIMATE EIGENDECOMPOSITION: Save approximate eigenvalues for eigenvalue error analysis
analysis.save_eigendecomposition_to_file(
    filename='error_demo_approx_eigendecomp.npz',
    operator='approximate',
    form='time_evolution',
    shift='shifted'
)
# Eigenenergies extracted from U_approx eigenphases
# Energy shifts are automatically handled for fair comparison with exact energies

# -------------------------------------------------------------------------------------------------
# Error Analysis - Three Independent Error Types
# -------------------------------------------------------------------------------------------------

# 1. EIGENVALUE ERRORS: Compare λ(H_exact) vs λ(H_approx)
#    - Reports errors for ALL eigenstates (ground + excited)
#    - Absolute errors: |λ_exact - λ_approx|
#    - Relative errors: |λ_exact - λ_approx| / |λ_exact|
#    - Typical values: 0.01-0.1% for good approximations

analysis.enable_eigenvalue_errors = True

# 2. MATRIX NORM ERRORS: ||U_exact - U_approx||
#    - Frobenius norm: √(Σᵢⱼ |Uᵢⱼ_exact - Uᵢⱼ_approx|²)
#    - Spectral norm: max singular value of (U_exact - U_approx)
#    - Measures overall approximation quality across all matrix elements
#    - Typical Frobenius values: 0.01-0.1 for good approximations

analysis.error_matrix_norms = ["frobenius", "spectral"]  # Compute both norms

# 3. STATE-DEPENDENT ERRORS: ||U_exact|ψ⟩ - U_approx|ψ⟩||
#    - Most practically relevant: how accurate for the states you care about?
#    - Absolute error: ||U_exact|ψ⟩ - U_approx|ψ⟩||
#    - Relative error: ||U_exact|ψ⟩ - U_approx|ψ⟩|| / ||U_exact|ψ⟩||
#    - Tells you if your specific input state is evolved accurately
#    - Typical relative values: 0.1-1% for good approximations

analysis.error_state_inputs = "examples/Be-H_1.30_sto-6g_as-003-003_jw.npy"

# =================================================================================================
# HOW ERROR ANALYSIS WORKS
# =================================================================================================
#
# Error analysis compares the exact and approximate operators to quantify accuracy.
# All necessary conversions happen automatically behind the scenes.
#
# INPUT:
#   - H_exact: The exact Hamiltonian (from Pauli string representation)
#   - U_approx: The approximate time-evolution operator (from Trotter/LCU/DF algorithm)
#
# AUTOMATIC CONVERSIONS:
#   1. H_exact → U_exact: Matrix exponential computes U_exact = exp(-i * H * t)
#   2. Energy shift handling: Algorithms may use shifted energies internally;
#      this is automatically corrected for fair comparison
#   3. U_approx → H_approx: Eigenvalues are extracted and converted to energies
#
# OUTPUT:
#   - Eigenvalue errors: Compare all eigenenergies (ground + excited states)
#   - Matrix norms: Overall approximation quality (Frobenius and spectral norms)
#   - State errors: Accuracy for your specific input quantum state
#
# EFFICIENCY:
#   - Converted operators are cached and reused across different error types
#   - No redundant computations (e.g., matrix exponentials only computed once)
#
# You just configure what analyses you want - the conversions are automatic!
# =================================================================================================

# =================================================================================================
# EXPECTED OUTPUT
# =================================================================================================
# This configuration will produce:
#
# 1. Log file (error_analysis_demo.log) with:
#    - Operator conversion progress
#    - Unitarity checks (should be ~1e-15 for both exact and approximate)
#    - Error analysis results
#
# 2. Output files:
#    - error_demo_exact_hamiltonian.npz: The exact Hamiltonian matrix
#    - error_demo_exact_eigendecomp.npz: Eigenenergies and eigenvectors of H_exact
#    - error_demo_algorithm_output.npz: The full algorithm circuit (U_approx for time evolution)
#    - error_demo_approx_eigendecomp.npz: Eigenenergies extracted from U_approx
#    - error_analysis.npz: All error analysis results in one file
#
# 3. Error analysis results:
#    - Eigenvalue errors for all eigenstates (ground state usually most accurate)
#    - Matrix Frobenius norm: ~0.01-0.1 for good approximations
#    - Matrix spectral norm: worst-case error measure
#    - State-dependent errors: ~0.1-1% relative error for good approximations
#
# INTERPRETATION:
#    - Eigenvalue errors tell you if computed energies are accurate
#    - Matrix norms tell you overall approximation quality
#    - State errors tell you if YOUR specific quantum state is well-approximated
#    - All three should be small and consistent for a good approximation
# =================================================================================================
