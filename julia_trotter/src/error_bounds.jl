# error_bounds.jl
#
# Commutator-based error bounds for Trotter-Suzuki formulas.
# Based on Childs et al., "Theory of Trotter Error with Commutator Scaling",
# PRX Quantum 2, 010323 (2021).

using ArnoldiMethod
using LinearAlgebra, SparseArrays

include("quantum_utils.jl")  # For OP_from_string

# ----------------------------------------------------------------------
# Strang (2nd-order) commutator-bound helper (fast upper bound)
# ----------------------------------------------------------------------

"""
    commutator_bound_prefactor_fast(abs_coeffs; nqubits, norm_mode=:l1)

Compute prefactor C for Strang splitting error bound using inequality-based upper bounds.

Given abs_coeffs[j] = ‖H_j‖ for a decomposition H = Σ_j H_j, returns a
prefactor `C` such that the *single-step* Strang error obeys:

    ‖S2(dt) - exp(-i dt H)‖ ≤ (dt^3) * C

This is obtained from Prop. 16 (Eq. 152) in Childs et al. (PRX 2021)
by upper-bounding the nested-commutator norms using:

    ‖[A,[A,B]]‖ ≤ 4 ‖A‖^2 ‖B‖

# Arguments
- `abs_coeffs`: Vector of term norms ‖H_j‖
- `nqubits`: Number of qubits
- `norm_mode`: `:l1` (sum of norms) or `:fro` (Frobenius norm scaling)

# Returns
Prefactor C for single-step error bound
"""
function commutator_bound_prefactor_fast(
    abs_coeffs::Vector{Float64};
    nqubits::Int,
    norm_mode::Symbol = :l1,
)
    error("There is a better bound - use exact method for accurate estimates")
    m = length(abs_coeffs)
    m == 0 && return 0.0

    # Suffix bounds for ‖R_j‖ = ‖Σ_{k=j+1}^m H_k‖
    if norm_mode == :l1
        suffix = zeros(Float64, m + 1)
        for j in m:-1:1
            suffix[j] = suffix[j + 1] + abs_coeffs[j]
        end
        Rnorm = j -> suffix[j + 1]
    elseif norm_mode == :fro
        # ‖Σ c_k P_k‖ ≤ ‖·‖_F = √(2^n) √(Σ |c_k|^2)
        suffix2 = zeros(Float64, m + 1)
        for j in m:-1:1
            suffix2[j] = suffix2[j + 1] + abs_coeffs[j]^2
        end
        fro_factor = sqrt(2.0^nqubits)
        Rnorm = j -> fro_factor * sqrt(suffix2[j + 1])
    else
        error("Unknown norm_mode=$norm_mode (use :l1 or :fro)")
    end

    # Prop. 16 structure:
    #   (dt^3/12) Σ_j ‖[R_j,[R_j,H_j]]‖ + (dt^3/24) Σ_j ‖[H_j,[H_j,R_j]]‖
    # Upper bound:
    #   ‖[R,[R,Hj]]‖ ≤ 4‖R‖^2‖Hj‖
    #   ‖[Hj,[Hj,R]]‖ ≤ 4‖Hj‖^2‖R‖
    sum1 = 0.0
    sum2 = 0.0
    for j in 1:m
        hj = abs_coeffs[j]
        rj = Rnorm(j)
        sum1 += 4.0 * (rj^2) * hj
        sum2 += 4.0 * (hj^2) * rj
    end

    return (sum1 / 12.0) + (sum2 / 24.0)
end


# ----------------------------------------------------------------------
# Exact Prop. 16 prefactor via explicit sparse nested-commutators
# (Expensive but tight bounds, limited to ~12 qubits)
# ----------------------------------------------------------------------

@inline comm(A, B) = A * B - B * A

"""
    hermitian_opnorm(A; scale=1E2)

Compute operator norm of Hermitian matrix via largest eigenvalue.

# Arguments
- `A`: Sparse matrix (will be symmetrized)
- `scale`: Numerical scaling factor for eigensolve stability

# Returns
Operator norm ‖A‖
"""
function hermitian_opnorm(A::SparseMatrixCSC{ComplexF64, Int}; scale::Float64 = 1E2)::Float64
    # Numerically symmetrize before wrapping
    H = scale * Hermitian((A + A') / 2)
    n = size(A, 1)
    maxdim = min(40, n - 1)  # Ensure maxdim < matrix dimension
    vals, _ = partialschur(H; nev=1, which=:LM, maxdim=maxdim)
    return abs(real(vals.R[1])) / scale
end

"""
    commutator_bound_prefactor_exact(term_mats)

Compute exact Prop. 16 prefactor via explicit nested commutators.

Returns prefactor C such that:
    ‖S2(dt) - exp(-i dt Σ H_j)‖ ≤ dt^3 * C

Explicitly constructs:
    [R_j,[R_j,H_j]]  and  [H_j,[H_j,R_j]]
and computes their operator norms via eigensolves.

# Arguments
- `term_mats`: Vector of Hermitian sparse matrices H_j

# Returns
Prefactor C for single-step error bound

# Warning
Expensive for large systems (scales as number of terms × matrix operations).
Practical limit ~12 qubits.
"""
function commutator_bound_prefactor_exact(
    term_mats::Vector{SparseMatrixCSC{ComplexF64, Int}}
)::Float64
    m = length(term_mats)
    m == 0 && return 0.0
    N = size(term_mats[1], 1)
    R = spzeros(ComplexF64, N, N)
    sum1 = 0.0
    sum2 = 0.0

    for j in m:-1:1
        Hj = term_mats[j]
        # R is Σ_{k=j+1}^m H_k
        inner1 = comm(R, Hj)          # anti-Hermitian
        C1     = comm(R, inner1)      # Hermitian
        sum1  += hermitian_opnorm(C1)

        inner2 = comm(Hj, R)          # anti-Hermitian
        C2     = comm(Hj, inner2)     # Hermitian
        sum2  += hermitian_opnorm(C2)

        R += Hj
    end
    return (sum1 / 12.0) + (sum2 / 24.0)
end


"""
    commutator_error_bounds(ham, normalization, nqubits; nsteps_list, time=π)

Compute commutator-based upper bounds for full Strang evolution error.

For each nsteps in nsteps_list, computes:
    ‖(S2(dt))^nsteps - exp(-i time/normalization * H0)‖ ≤ nsteps * ‖S2(dt) - exp(-i dt H0)‖

using telescoping inequality for unitaries.

# Arguments
- `ham`: Hamiltonian dictionary (Pauli string => coefficient)
- `normalization`: Scaling factor (dt = time / (nsteps * normalization))
- `nqubits`: Number of qubits
- `nsteps_list`: List of Trotter step counts
- `time`: Total evolution time (default π)

# Environment Variables
- `COMM_METHOD`: "exact" (default) or "fast"
  - exact: Compute nested commutators explicitly (≤12 qubits)
  - fast: Use inequality-based upper bound (looser but cheap)
- `COMM_NORM_MODE`: "l1" or "fro" (only for fast method)

# Returns
- `bounds`: Vector of error bounds (one per nsteps)
- `method`: String indicating which method was used
"""
function commutator_error_bounds(
    ham::Dict,
    normalization::Real,
    nqubits::Int;
    nsteps_list::Vector{Int} = [1],
    time = π,
)
    id_key = "I"^nqubits
    H_terms = [(k, v) for (k, v) in ham if k != id_key]

    # Method selection
    comm_method = lowercase(get(ENV, "COMM_METHOD", "exact"))

    pref = 0.0
    comm_method_used = comm_method

    if comm_method == "exact"
        # Guard rails: exact method limited to small systems
        if nqubits > 12
            @warn "nqubits=$nqubits > 12, switching to fast method"
            comm_method_used = "fast"
        else
            term_mats = SparseMatrixCSC{ComplexF64, Int}[]
            for (k, v) in H_terms
                @assert imag(v) ≈ 0.0 "Coefficient not real for $k"
                push!(term_mats, real(v) * sparse(OP_from_string(k)))
            end
            pref = commutator_bound_prefactor_exact(term_mats)
        end
    end

    if comm_method_used == "fast"
        coeffs = Float64[]
        for (k, v) in H_terms
            @assert imag(v) ≈ 0.0 "Coefficient not real for $k"
            push!(coeffs, abs(real(v)))
        end
        norm_mode = Symbol(get(ENV, "COMM_NORM_MODE", "l1"))
        pref = commutator_bound_prefactor_fast(coeffs; nqubits=nqubits, norm_mode=norm_mode)
    end

    bounds = Float64[]
    for nsteps in nsteps_list
        dt = time / (nsteps * normalization)
        local_bound = (dt^3) * pref
        total_bound = nsteps * local_bound
        push!(bounds, total_bound)
    end

    return bounds, comm_method_used
end


"""
    reference_trotter_unitary(ham, normalization, nqubits; numsteps=10, time=π)

Construct reference second-order Trotter unitary matrix for benchmarking.

Computes: exp(-i time/2) * [S2(dt)]^numsteps
where dt = time / (numsteps * normalization) and S2 is Strang splitting.

# Arguments
- `ham`: Hamiltonian dictionary (Pauli string => coefficient)
- `normalization`: Scaling factor
- `nqubits`: Number of qubits
- `numsteps`: Number of Trotter steps
- `time`: Total evolution time

# Returns
Unitary matrix (sparse)

# Note
This constructs the full unitary matrix. For large systems, use
statevector methods instead.
"""
function reference_trotter_unitary(
    ham::Dict,
    normalization::Real,
    nqubits::Int;
    numsteps::Int=10,
    time=π
)
    U_s      = sparse(I, 2^nqubits, 2^nqubits) .+ 0.0im
    Identity = sparse(I, 2^nqubits, 2^nqubits) .+ 0.0im
    id_key   = "I"^nqubits

    # Ignore identity term in Trotter product
    H_terms = [(k, v) for (k, v) in ham if k != id_key]

    dt = time / (numsteps * normalization)

    for _step in 1:numsteps
        # Forward
        for (k, v) in H_terms
            θ = real(v) * dt / 2
            U_s = (cos(θ) * Identity - im * sin(θ) * sparse(OP_from_string(k))) * U_s
        end
        # Backward
        for (k, v) in reverse(H_terms)
            θ = real(v) * dt / 2
            U_s = (cos(θ) * Identity - im * sin(θ) * sparse(OP_from_string(k))) * U_s
        end
    end

    # Global phase from +0.5I shift (common QPE convention)
    return exp(-im * time / 2) * U_s
end
