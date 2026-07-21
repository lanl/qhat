"""
QHAT Error Analysis Demo - Phase 2 Operator Framework

This configuration demonstrates the error analysis capabilities with focus on
the Phase 2 operator framework improvements (2026-07).

Key improvements in Phase 2:
1. Correct comparison: U_exact vs U_approx (not H_exact vs U_approx)
2. Automatic operator conversion via OperatorRepresentation class
3. Energy shift handling (approximate algorithms use H_shifted for phase reduction)
4. Lazy evaluation and caching for efficiency

Expected error values (for good Trotter approximations):
- Eigenvalue errors: ~0.01-0.1% relative error
- Frobenius norm: ~0.01-0.1 (down from ~25 in Phase 0 bug)
- State relative error: ~0.1-1% (down from ~147% in Phase 0 bug)

Usage:
    python3.11 -m qhat.analysis.driver examples/config_error_analysis_demo.py
"""

# Energy error budget: split between Trotterization and phase estimation
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
    energy_error=0.5 * energy_error,  # Split error budget
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

# Resource estimation is fast and always useful
analysis.resource_estimator = "pyLIQTR"

# -------------------------------------------------------------------------------------------------
# Matrix Output - Required for Error Analysis
# -------------------------------------------------------------------------------------------------
# Error analysis needs both exact Hamiltonian and approximate unitary matrices

analysis.exact_matrix_output_file = "error_demo_exact_hamiltonian.npz"
# This is the exact H_exact (no approximations)
# Phase 2: OperatorRepresentation will convert this to U_exact via exp(-i*H*t)

analysis.matrix_output_file = "error_demo_approx_unitary.npz"
# This is the approximate U_approx from Trotter decomposition
# Note: This is actually U_s,approx (energy-shifted) internally
# Phase 2: OperatorRepresentation will convert this to U_approx (unshifted) for comparison

# -------------------------------------------------------------------------------------------------
# Eigendecomposition - Required for Eigenvalue Errors
# -------------------------------------------------------------------------------------------------
# Compute eigenvalues of both exact and approximate operators

analysis.eigendecomposition_matrices = "both"
# This produces:
#   - Eigenenergies of H_exact (ground state, excited states)
#   - Eigenenergies of H_approx (extracted from U_approx eigenphases)
# Phase 2: Automatic energy shift correction ensures proper comparison

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
#    - Phase 2 fix: Now correctly compares U vs U (not H vs U!)
#    - Typical Frobenius values: 0.01-0.1 (was ~25 in Phase 0 bug)

analysis.error_matrix_norms = ["frobenius", "spectral"]  # Compute both norms

# 3. STATE-DEPENDENT ERRORS: ||U_exact|ψ⟩ - U_approx|ψ⟩||
#    - Most practically relevant: how accurate for states you care about?
#    - Absolute error: ||U_exact|ψ⟩ - U_approx|ψ⟩||
#    - Relative error: ||U_exact|ψ⟩ - U_approx|ψ⟩|| / ||U_exact|ψ⟩||
#    - Phase 2 fix: Now correctly evolves with U (not H!)
#    - Typical relative values: 0.1-1% (was ~147% in Phase 0 bug)

analysis.error_state_inputs = "examples/Be-H_1.30_sto-6g_as-003-003_jw.npy"

# =================================================================================================
# UNDERSTANDING THE OPERATOR CONVERSIONS (Phase 2 Internals)
# =================================================================================================
#
# When error analysis runs, it performs these conversions internally:
#
# 1. INPUT OPERATORS:
#    - H_exact: Exact Hamiltonian (unshifted), from exact_matrix_output
#    - U_s,approx: Approximate time-evolution (energy-shifted), from algorithm
#
# 2. WHAT NEEDS TO BE COMPARED:
#    - For eigenvalue errors: λ(H_exact) vs λ(H_approx)
#    - For matrix norm errors: U_exact vs U_approx (both unshifted)
#    - For state errors: U_exact|ψ⟩ vs U_approx|ψ⟩ (both unshifted)
#
# 3. CONVERSIONS PERFORMED AUTOMATICALLY:
#    a. H_exact → U_exact:
#       - Compute U_exact = exp(-i * H_exact * t) via matrix exponential
#       - Verify unitarity: ||U†U - I||_F ≈ 1e-15
#
#    b. U_s,approx → U_approx:
#       - Remove energy shift: U_approx = exp(-i*E*t) * U_s,approx
#       - Verify unitarity: ||U†U - I||_F ≈ 1e-15
#
#    c. U_approx → H_approx (for eigenvalue errors):
#       - Extract eigenvalues of U_approx: exp(-i*λ_H*t)
#       - Convert to energies: λ_H = i*log(λ_U)/t
#
# 4. CACHING FOR EFFICIENCY:
#    - Converted operators are cached
#    - Matrix norm errors and state errors share the same U_exact and U_approx
#    - No redundant matrix exponential computations
#
# 5. WHY ENERGY SHIFTS?
#    - Approximate algorithms use H_shifted = H - E_min·I
#    - This keeps all eigenphases in [0, 2π), reducing phase qubit requirements
#    - Error analysis automatically handles this, comparing unshifted operators
#
# All of this happens automatically - you just configure what analyses you want!
# =================================================================================================

# =================================================================================================
# EXPECTED OUTPUT
# =================================================================================================
# This configuration will produce:
#
# 1. Log file with detailed progress and operator conversion info:
#    - "Converting H_exact → U_exact (unshifted)"
#    - "Converting U_s,approx → U_approx (unshifted)"
#    - "U_exact unitarity check: ||U†U - I||_F = ~1e-15"
#    - "U_approx unitarity check: ||U†U - I||_F = ~1e-15"
#    - "Frobenius norm ||U_exact - U_approx||_F: ~0.01-0.1"
#    - "State error ||U_exact|ψ⟩ - U_approx|ψ⟩||: relative ~0.1-1%"
#
# 2. Matrix files:
#    - error_demo_exact_hamiltonian.npz (H_exact)
#    - error_demo_approx_unitary.npz (U_s,approx, energy-shifted)
#    - exact_eigendecomposition.npz (eigenenergies and eigenvectors of H_exact)
#    - approximate_eigendecomposition.npz (eigenenergies of H_approx, corrected)
#
# 3. Error analysis results (in results dict and log):
#    - Eigenvalue errors for all states (ground state usually most accurate)
#    - Matrix Frobenius norm: measures total error across all matrix elements
#    - Matrix spectral norm: measures worst-case error (max singular value)
#    - State-dependent errors: how accurate for this specific quantum state
#
# The error values will be physically meaningful and consistent across all three types!
# =================================================================================================
