# =============================================================================
# trotter_expts.jl
#
# Comprehensive Trotterization Error Analysis Demonstration
#
# Demonstrates:
#   1. State vector errors (1st, 2nd, 4th order Trotter vs Chebyshev)
#   2. Operator errors (unitary matrix comparison)
#   3. Commutator bounds (theoretical error predictions)
#
# Usage: julia --project=. trotter_expts.jl [hamiltonian_file.dat]
#        If no file specified, uses default Li-Li 4-qubit system
# =============================================================================

# Includes
using LinearAlgebra, SparseArrays, Printf, Arpack
include("src/parser.jl")
include("src/quantum_utils.jl")
include("src/hamiltonian_utils.jl")
include("src/statevector_simulators.jl")
include("src/densesimulators.jl")
include("src/tensorsimulators.jl")
include("src/error_bounds.jl")


# Constants
const T_TOTAL = 1.0           # Total evolution time
const NSTEPS_LIST = [10, 20, 50, 100]  # Trotter step counts to test
const CHEBY_ORDER = 10        # Chebyshev expansion order

function main(filepath::String)
    println("\n" * "="^70)
    println("  TROTTERIZATION ERROR ANALYSIS DEMONSTRATION")
    println("="^70)

    # ========== Part 1: Parse Hamiltonian ==========
    println("\n[1] PARSING HAMILTONIAN")
    println("-"^70)

    meta, ham = parse_hamiltonian_file(filepath)

    nqubits = parse(Int, meta["number of qubits"])
    nterms = parse(Int, meta["number of terms in sum of Pauli strings"])
    nelectrons = parse(Int, meta["number of active, occupied, single-occupancy orbitals"])

    println("  File: ", basename(filepath))
    println("  Qubits: ", nqubits)
    println("  Hamiltonian terms: ", nterms)
    println("  Active electrons: ", nelectrons)
    println("  Ground state energy: ", meta["smallest eigenvalue"], " Ha")

    # Normalization
    norm_result = normalize_hamiltonian(meta, ham)
    println("  L1 norm (excl. identity): ", @sprintf("%.4f", norm_result.norm_bound))
    println("  Normalization factor: ", @sprintf("%.4f", norm_result.normalization))

    # ========== Part 2: State Vector Errors ==========
    println("\n[2] STATE VECTOR ERRORS")
    println("-"^70)
    println("Computing Trotter evolution vs Chebyshev reference...")

    # Build terms
    ham_result = build_hamiltonian_terms(ham, meta; normalize=true, scale_by_pi=true)
    H_terms = ham_result.H_terms

    # Initial state (Hartree-Fock)
    stateHF = construct_hf_state(nqubits, nelectrons)

    # Reference evolution (Chebyshev)
    println("  Computing reference (Chebyshev order=$(CHEBY_ORDER))...")
    psi_exact, _, _ = chebyshev_simulation(
        stateHF, T_TOTAL,
        x -> hamiltonian_matvec(x, H_terms),
        π, [];
        order=CHEBY_ORDER,
        returnstate=true
    )

    # Compute errors for each Trotter order
    println("\n  State errors: ‖ψ_trotter - ψ_exact‖")
    println("  " * "-"^66)
    println(@sprintf("  %-8s  %-10s  %-12s  %-12s  %-12s",
                     "nsteps", "dt", "1st order", "2nd order", "4th order"))
    println("  " * "-"^66)

    for nsteps in NSTEPS_LIST
        dt = T_TOTAL / nsteps

        psi_1 = first_order_trotter_statevec(H_terms, stateHF, dt, nsteps)
        err_1 = norm(psi_1 - psi_exact)

        psi_2 = second_order_trotter_statevec(H_terms, stateHF, dt, nsteps)
        err_2 = norm(psi_2 - psi_exact)

        psi_4 = fourth_order_trotter_statevec(H_terms, stateHF, dt, nsteps)
        err_4 = norm(psi_4 - psi_exact)

        println(@sprintf("  %-8d  %-10.4f  %-12.2e  %-12.2e  %-12.2e",
                         nsteps, dt, err_1, err_2, err_4))
    end

    # ========== Part 3: Operator Errors ==========
    println("\n[3] OPERATOR ERRORS")
    println("-"^70)
    println("Computing unitary matrix errors ‖U_trotter - U_exact‖...")

    # Build full H matrix
    H_full = build_sparse_hamiltonian(ham, nqubits)
    H_normalized = (Hermitian(H_full) - norm_result.shift * I) / norm_result.normalization
    H_qpe = (H_normalized + 0.5*I) * π

    # Exact unitary via Chebyshev
    println("  Computing U_exact (Chebyshev)...")
    U_exact = chebyshev_exponentiation(H_qpe, CHEBY_ORDER, 1.0, π)

    println("\n  Operator errors (spectral norm):")
    println("  " * "-"^40)
    println(@sprintf("  %-8s  %-12s", "nsteps", "2nd order"))
    println("  " * "-"^40)

    op_errors = Float64[]
    for nsteps in NSTEPS_LIST
        # Build Trotter unitary
        U_trotter = reference_trotter_unitary(
            ham, norm_result.normalization, nqubits;
            numsteps=nsteps, time=π
        )

        # Compute operator norm of difference
        E_op = U_trotter - U_exact
        if nqubits <= 10
            # Small system: use dense norm
            op_err = opnorm(Matrix(E_op))
        else
            # Large system: use Arpack
            vals, _ = Arpack.eigs(E_op; nev=1, which=:LM)
            op_err = abs(vals[1])
        end

        push!(op_errors, op_err)
        println(@sprintf("  %-8d  %-12.2e", nsteps, op_err))
    end

    # ========== Part 4: Commutator Bounds ==========
    println("\n[4] COMMUTATOR ERROR BOUNDS")
    println("-"^70)
    println("Computing theoretical bounds (Childs et al. PRX 2021)...")

    # Compute bounds
    bounds, method_used = commutator_error_bounds(
        ham, norm_result.normalization, nqubits;
        nsteps_list=NSTEPS_LIST,
        time=π
    )

    println("  Method used: ", method_used)
    println("\n  Bounds vs Actual Operator Errors:")
    println("  " * "-"^55)
    println(@sprintf("  %-8s  %-15s  %-15s  %-10s",
                     "nsteps", "Actual", "Bound", "Ratio"))
    println("  " * "-"^55)

    for (i, nsteps) in enumerate(NSTEPS_LIST)
        ratio = bounds[i] / op_errors[i]
        println(@sprintf("  %-8d  %-15.2e  %-15.2e  %-10.2f",
                         nsteps, op_errors[i], bounds[i], ratio))
    end

    # ========== Part 5: Summary ==========
    println("\n[5] SUMMARY & ANALYSIS")
    println("-"^70)

    # Verify error scaling
    println("\nError Scaling Analysis (2nd order Trotter):")
    if length(NSTEPS_LIST) >= 2
        ratio_10_20 = op_errors[1] / op_errors[2]
        expected_ratio = (NSTEPS_LIST[2] / NSTEPS_LIST[1])^2  # dt² scaling
        println("  Observed ratio (10/20 steps): ", @sprintf("%.2f", ratio_10_20))
        println("  Expected (dt² scaling):       ", @sprintf("%.2f", expected_ratio))

    end

    # Bound quality
    println("\nCommutator Bound Quality:")
    avg_ratio = sum(bounds ./ op_errors) / length(bounds)
    println("  Average bound/actual ratio: ", @sprintf("%.2f", avg_ratio))


end

# CLI: Default to smallest Hamiltonian if no args
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) >= 1
        filepath = ARGS[1]
    else
        # Default to smallest system
        filepath = "Li-Li_jw/Li-Li_2.90_hgbs-5_as-002-002_jw.dat"
        println("\nNo file specified. Using default: ", filepath)
        println("Usage: julia --project=. trotter_expts.jl [hamiltonian_file.dat]\n")
    end

    main(filepath)
end
