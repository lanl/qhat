# error_bounds.jl
#
# Commutator-based error bounds for Trotter-Suzuki formulas.
# Based on Childs et al., "Theory of Trotter Error with Commutator Scaling",
# Phys. Rev. X 11, 011020 (2021).

using ArnoldiMethod
using LinearAlgebra, SparseArrays

if !isdefined(@__MODULE__, :ordered_hamiltonian_terms)
    include("hamiltonian_utils.jl")
end

function pauli_anticommuting_suffix_weights(
    pauli_strings::Vector{String},
    abs_coeffs::Vector{Float64},
)
    length(pauli_strings) == length(abs_coeffs) || throw(DimensionMismatch(
        "Pauli strings and coefficients must have the same length",
    ))
    isempty(pauli_strings) && return Float64[]

    nqubits = length(first(pauli_strings))
    nchunks = cld(nqubits, 64)
    x_masks = zeros(UInt64, nchunks, length(pauli_strings))
    z_masks = zeros(UInt64, nchunks, length(pauli_strings))

    for (term, pstr) in enumerate(pauli_strings)
        length(pstr) == nqubits || throw(DimensionMismatch(
            "All Pauli strings must have the same length",
        ))
        for (site, pauli) in enumerate(pstr)
            chunk = (site - 1) ÷ 64 + 1
            bit = UInt64(1) << ((site - 1) % 64)
            if pauli == 'X'
                x_masks[chunk, term] |= bit
            elseif pauli == 'Y'
                x_masks[chunk, term] |= bit
                z_masks[chunk, term] |= bit
            elseif pauli == 'Z'
                z_masks[chunk, term] |= bit
            elseif pauli != 'I'
                throw(ArgumentError("Invalid Pauli '$pauli' in '$pstr'"))
            end
        end
    end

    anticommuting_weights = zeros(Float64, length(pauli_strings))
    for j in 1:length(pauli_strings)-1
        for k in j+1:length(pauli_strings)
            parity = false
            @inbounds for chunk in 1:nchunks
                crossings = (x_masks[chunk, j] & z_masks[chunk, k]) ⊻
                            (z_masks[chunk, j] & x_masks[chunk, k])
                parity = xor(parity, isodd(count_ones(crossings)))
            end
            parity && (anticommuting_weights[j] += abs_coeffs[k])
        end
    end
    return anticommuting_weights
end

# ----------------------------------------------------------------------
# Strang (2nd-order) commutator-bound helper (fast upper bound)
# ----------------------------------------------------------------------

"""
    commutator_bound_prefactor_fast(pauli_strings, abs_coeffs;
                                    nqubits, norm_mode=:l1)

Compute prefactor C for Strang splitting error bound using inequality-based upper bounds.

Given `abs_coeffs[j] = ‖H_j‖` and the corresponding Pauli strings for a
decomposition `H = Σ_j H_j`, returns a prefactor `C` such that the
*single-step* Strang error obeys:

    ‖S2(dt) - exp(-i dt H)‖ ≤ (dt^3) * C

This is obtained from Proposition 10, Eq. (121), in Childs et al. (PRX 2021)
by upper-bounding the nested-commutator norms using:

    ‖[A,[A,B]]‖ ≤ 4 ‖A‖^2 ‖B‖

# Arguments
- `abs_coeffs`: Vector of term norms ‖H_j‖
- `pauli_strings`: Pauli string for each coefficient
- `nqubits`: Number of qubits
- `norm_mode`: `:l1` (sum of norms) or `:fro` (Frobenius norm scaling)

# Returns
Prefactor C for single-step error bound
"""
function commutator_bound_prefactor_fast(
    pauli_strings::Vector{String},
    abs_coeffs::Vector{Float64};
    nqubits::Int,
    norm_mode::Symbol = :l1,
)
    m = length(abs_coeffs)
    m == 0 && return 0.0
    anticommuting_weights = pauli_anticommuting_suffix_weights(
        pauli_strings, abs_coeffs
    )

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

    # Proposition 10 structure:
    #   (dt^3/12) Σ_j ‖[R_j,[R_j,H_j]]‖ + (dt^3/24) Σ_j ‖[H_j,[H_j,R_j]]‖
    # If Aj is the coefficient weight of suffix terms that anticommute with
    # Hj, then ‖[R,Hj]‖ ≤ 2‖Hj‖Aj. Consequently,
    #   ‖[R,[R,Hj]]‖ ≤ 4‖R‖‖Hj‖Aj
    #   ‖[Hj,[Hj,R]]‖ ≤ 4‖Hj‖²Aj.
    sum1 = 0.0
    sum2 = 0.0
    for j in 1:m
        hj = abs_coeffs[j]
        rj = Rnorm(j)
        aj = anticommuting_weights[j]
        sum1 += 4.0 * rj * hj * aj
        sum2 += 4.0 * (hj^2) * aj
    end

    return (sum1 / 12.0) + (sum2 / 24.0)
end


# ----------------------------------------------------------------------
# Exact first- and second-order prefactors via sparse commutators
# (Expensive but tight bounds, limited to ~12 qubits)
# ----------------------------------------------------------------------

@inline comm(A, B) = A * B - B * A

function first_order_commutator_bound_prefactor_fast(
    pauli_strings::Vector{String},
    abs_coeffs::Vector{Float64},
)
    anticommuting_weights = pauli_anticommuting_suffix_weights(
        pauli_strings, abs_coeffs
    )
    return sum(abs_coeffs .* anticommuting_weights)
end

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
    # Small commutator matrices often have large degeneracies. A dense solve
    # is deterministic here and avoids an iterative solve missing an extremal
    # invariant subspace, which would underestimate a purported bound.
    if n <= 64
        return maximum(abs, eigvals(Hermitian(Matrix(H)))) / scale
    end

    # ArnoldiMethod requires its default `mindim` to be no larger than
    # `maxdim`.  Letting the Krylov space span a small matrix keeps that
    # invariant and is also exact for the small systems used in validation.
    maxdim = min(40, n)
    vals, _ = partialschur(H; nev=1, which=:LM, maxdim=maxdim)
    return abs(real(vals.R[1])) / scale
end

"""
    first_order_commutator_bound_prefactor_exact(term_mats)

Compute the first-order Lie-Trotter prefactor

    C₁ = (1/2) ∑ⱼ ‖[∑ₖ₌ⱼ₊₁ Hₖ, Hⱼ]‖,

so that one step satisfies `‖S₁(dt) - exp(-im*dt*H)‖ ≤ dt²*C₁`.
This is Proposition 9, Eq. (120), of Childs et al., Phys. Rev. X 11,
011020 (2021), DOI: 10.1103/PhysRevX.11.011020.
"""
function first_order_commutator_bound_prefactor_exact(
    term_mats::Vector{SparseMatrixCSC{ComplexF64, Int}}
)::Float64
    isempty(term_mats) && return 0.0
    N = size(term_mats[1], 1)
    R = spzeros(ComplexF64, N, N)
    commutator_sum = 0.0

    for j in length(term_mats):-1:1
        Hj = term_mats[j]
        # [R,Hj] is anti-Hermitian, so im*[R,Hj] is Hermitian with the
        # same operator norm.
        commutator_sum += hermitian_opnorm(im * comm(R, Hj))
        R += Hj
    end

    return commutator_sum / 2
end

"""
    first_order_commutator_error_bounds(ham, normalization, nqubits;
                                        nsteps_list=[1], time=π,
                                        term_ordering=:magnitude)

Compute rigorous additive-error bounds for first-order Lie-Trotter evolution.
For `nsteps` over total time `time`, the single-step bound is applied with
`dt = time/(nsteps*normalization)` and accumulated using the unitary
telescoping inequality.
`term_ordering` must match the ordering used by the product formula.
"""
function first_order_commutator_error_bounds(
    ham::Dict,
    normalization::Real,
    nqubits::Int;
    nsteps_list::Vector{Int}=[1],
    time::Real=π,
    term_ordering=:magnitude,
)
    all(>(0), nsteps_list) || throw(ArgumentError("nsteps values must be positive"))
    normalization > 0 || throw(ArgumentError("normalization must be positive"))

    term_mats = SparseMatrixCSC{ComplexF64,Int}[]
    pauli_strings = String[]
    abs_coeffs = Float64[]
    for (pauli, coefficient) in ordered_hamiltonian_terms(
        ham, nqubits; term_ordering=term_ordering
    )
        @assert isapprox(imag(coefficient), 0.0) "Coefficient not real for $pauli"
        push!(pauli_strings, pauli)
        push!(abs_coeffs, abs(real(coefficient)))
        if nqubits <= 12
            push!(term_mats, real(coefficient) * sparse(OP_from_string(pauli)))
        end
    end

    prefactor = nqubits <= 12 ?
                first_order_commutator_bound_prefactor_exact(term_mats) :
                first_order_commutator_bound_prefactor_fast(pauli_strings, abs_coeffs)
    return [
        nsteps * (time / (nsteps * normalization))^2 * prefactor
        for nsteps in nsteps_list
    ]
end

"""
    commutator_bound_prefactor_exact(term_mats)

Compute the exact Proposition 10 prefactor via explicit nested commutators.

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
    commutator_error_bounds(ham, normalization, nqubits; nsteps_list, time=π,
                            term_ordering=:magnitude)

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
- `term_ordering`: Ordering accepted by `ordered_hamiltonian_terms`; it must
  match the product formula

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
    term_ordering=:magnitude,
)
    H_terms = ordered_hamiltonian_terms(
        ham, nqubits; term_ordering=term_ordering
    )

    # Method selection
    comm_method = lowercase(get(ENV, "COMM_METHOD", "exact"))
    @assert comm_method in ("fast", "exact") "COMM_METHOD must be 'fast' or 'exact', got '$comm_method'"

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
        pauli_strings = String[]
        coeffs = Float64[]
        for (k, v) in H_terms
            @assert imag(v) ≈ 0.0 "Coefficient not real for $k"
            push!(pauli_strings, k)
            push!(coeffs, abs(real(v)))
        end
        norm_mode = Symbol(get(ENV, "COMM_NORM_MODE", "l1"))
        pref = commutator_bound_prefactor_fast(
            pauli_strings, coeffs; nqubits=nqubits, norm_mode=norm_mode
        )
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
    reference_trotter_unitary(ham, normalization, nqubits; numsteps=10, time=π,
                              term_ordering=:magnitude)

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
    time=π,
    term_ordering=:magnitude,
)
    U_s      = sparse(I, 2^nqubits, 2^nqubits) .+ 0.0im
    Identity = sparse(I, 2^nqubits, 2^nqubits) .+ 0.0im
    H_terms = ordered_hamiltonian_terms(
        ham, nqubits; term_ordering=term_ordering
    )

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
