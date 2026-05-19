# hamiltonian_utils.jl
#
# Utilities for processing and preparing Hamiltonians for simulation.
# Handles normalization, term construction, and state preparation.

using SparseArrays, LinearAlgebra

include("quantum_utils.jl")  # For OP_from_string

"""
    normalize_hamiltonian(meta, ham)

Extract energy shift and compute normalization factor from metadata.

The normalization convention is:
    normalization = 2 * max(1.0, norm_bound)
where norm_bound is the L1 norm minus the identity coefficient.

# Arguments
- `meta`: Metadata dictionary from parsed Hamiltonian file
- `ham`: Hamiltonian dictionary (Pauli string => coefficient)

# Returns
Named tuple with fields:
- `shift`: Energy offset (coefficient of identity term)
- `norm_bound`: L1 norm excluding identity
- `normalization`: Scaling factor (2 * max(1, norm_bound))
"""
function normalize_hamiltonian(meta::Dict{String,String}, ham::Dict{String,ComplexF64})
    nqubits = parse(Int, meta["number of qubits"])
    id_key = "I"^nqubits

    shift = haskey(ham, id_key) ? real(ham[id_key]) : 0.0
    norm_bound = parse(Float64, meta["one-norm of sum of Pauli strings"]) - abs(shift)
    normalization = 2 * (norm_bound < 1.0 ? 1.0 : norm_bound)

    return (shift=shift, norm_bound=norm_bound, normalization=normalization)
end

"""
    build_hamiltonian_terms(ham, meta; normalize=true, scale_by_pi=true)

Build list of (coefficient, Pauli_matrix) tuples for state-vector simulation.

This function:
1. Extracts all non-identity terms from the Hamiltonian
2. Optionally normalizes coefficients
3. Optionally scales by π (common for QPE convention)
4. Converts Pauli strings to sparse matrices
5. Sorts terms by coefficient for determinism

# Arguments
- `ham`: Hamiltonian dictionary (Pauli string => coefficient)
- `meta`: Metadata dictionary (needed for normalization)
- `normalize`: If true, divide coefficients by normalization factor
- `scale_by_pi`: If true, multiply coefficients by π

# Returns
Named tuple with fields:
- `H_terms`: Vector of `(coeff::Float64, P::SparseMatrixCSC)` tuples
- `normalization`: The normalization factor used
- `shift`: The energy shift (identity coefficient)
"""
function build_hamiltonian_terms(
    ham::Dict{String,ComplexF64},
    meta::Dict{String,String};
    normalize::Bool=true,
    scale_by_pi::Bool=true
)
    nqubits = parse(Int, meta["number of qubits"])
    id_key = "I"^nqubits

    # Get normalization if needed
    if normalize
        norm_info = normalize_hamiltonian(meta, ham)
        shift = norm_info.shift
        normalization = norm_info.normalization
    else
        shift = haskey(ham, id_key) ? real(ham[id_key]) : 0.0
        normalization = 1.0
    end

    # Build terms
    H_terms = Tuple{Float64, SparseMatrixCSC{ComplexF64, Int}}[]
    for (k, v) in ham
        if k != id_key
            coeff = real(v)
            if normalize
                coeff /= normalization
            end
            if scale_by_pi
                coeff *= π
            end
            push!(H_terms, (coeff, sparse(OP_from_string(k))))
        end
    end

    # Sort by coefficient for determinism
    sort!(H_terms, by=x -> x[1])

    return (H_terms=H_terms, normalization=normalization, shift=shift)
end

"""
    construct_hf_state(nqubits, nelectrons)

Construct Hartree-Fock state in Jordan-Wigner encoding.

The HF state has the first `nelectrons` orbitals occupied:
    |11...100...0⟩  (nelectrons ones, then zeros)

In the computational basis, this corresponds to:
    kint = Σᵢ 2^(nqubits - nelectrons + i - 1) + 1   for i ∈ [1, nelectrons]

# Arguments
- `nqubits`: Number of qubits
- `nelectrons`: Number of occupied orbitals

# Returns
Normalized state vector with one nonzero entry at HF configuration
"""
function construct_hf_state(nqubits::Int, nelectrons::Int)
    kint = sum(2^ii for ii in nqubits-nelectrons:nqubits-1) + 1
    state = zeros(ComplexF64, 2^nqubits)
    state[kint] = 1.0
    return state
end

"""
    resolve_hamiltonian_path(filepath)

Resolve Hamiltonian file path handling cluster/local environment differences.

Tries multiple strategies:
1. Direct path (if file exists)
2. Relative to repo root (one level up from src/)
3. Extract path after /model_lib/ marker and search from repo root

# Arguments
- `filepath`: Path to Hamiltonian file (possibly incomplete)

# Returns
Resolved absolute path

# Throws
Error if file cannot be found by any strategy
"""
function resolve_hamiltonian_path(filepath::AbstractString)
    if isfile(filepath)
        return filepath
    end

    # Try relative to repo root (assuming we're in src/ or expts/)
    repo_root_candidate = normpath(joinpath(@__DIR__, "..", filepath))
    if isfile(repo_root_candidate)
        return repo_root_candidate
    end

    # Try extracting path after /model_lib/ marker
    marker = "/model_lib/"
    idx = findfirst(marker, filepath)
    if idx !== nothing
        relpath = filepath[last(idx)+1:end]
        candidate = normpath(joinpath(@__DIR__, "..", "model_lib", relpath))
        if isfile(candidate)
            return candidate
        end
    end

    error("Hamiltonian file not found: $filepath")
end

"""
    build_sparse_hamiltonian(ham, nqubits)

Build full sparse Hamiltonian matrix from Pauli string representation.

# Arguments
- `ham`: Hamiltonian dictionary (Pauli string => coefficient)
- `nqubits`: Number of qubits

# Returns
Sparse matrix representing full Hamiltonian
"""
function build_sparse_hamiltonian(ham::Dict{String,ComplexF64}, nqubits::Int)
    H = spzeros(ComplexF64, 2^nqubits, 2^nqubits)
    id_key = "I"^nqubits

    for (k, v) in ham
        @assert imag(v) ≈ 0.0 "Coefficient not real for $k"
        if k == id_key
            H += real(v) * sparse(I, 2^nqubits, 2^nqubits)
        else
            H += real(v) * sparse(OP_from_string(k))
        end
    end

    return H
end
