# trotter_spectral_energy.jl
#
# Shared helpers for extracting Trotter-derived ground energies from the
# largest eigenvalues of U + U'.

using ArnoldiMethod, Arpack, LinearAlgebra, LinearMaps, SparseArrays

include("hamiltonian_utils.jl")
include("error_bounds.jl")
include("statevector_simulators.jl")

const TROTTER_SPECTRAL_DEFAULT_NEV = 4

function cosine_value_to_energy(
    cosine_eigenvalue::Real,
    shift::Real,
    normalization::Real,
    time::Real=pi
)
    theta = acos(clamp(cosine_eigenvalue / 2, -1.0, 1.0))
    return shift + normalization * (theta / time - 0.5)
end

function trotter_candidate_count(n::Int; nev::Int=TROTTER_SPECTRAL_DEFAULT_NEV)
    return min(nev, n - 1)
end

function trotter_krylov_dimension(n::Int, nev::Int)
    return min(40, max(nev + 2, min(n - 1, 2 * nev + 1)))
end

function build_trotter_terms(ham::Dict{String,ComplexF64}, nqubits::Int)
    id_key = "I"^nqubits
    terms = Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}[]

    for (k, v) in ham
        k == id_key && continue
        @assert isapprox(imag(v), 0.0) "Coefficient not real for $k"
        push!(terms, (real(v), sparse(OP_from_string(k))))
    end

    return terms
end

function apply_first_order_trotter_step!(
    psi::Vector{ComplexF64},
    terms::Vector{Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}},
    dt::Real,
    normalization::Real
)
    normalized_dt = dt / normalization

    for (coeff, P) in terms
        apply_pauli_rotation!(psi, P, coeff, normalized_dt)
    end

    return psi
end

function apply_first_order_trotter_adjoint_step!(
    psi::Vector{ComplexF64},
    terms::Vector{Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}},
    dt::Real,
    normalization::Real
)
    normalized_dt = -dt / normalization

    for (coeff, P) in reverse(terms)
        apply_pauli_rotation!(psi, P, coeff, normalized_dt)
    end

    return psi
end

function apply_second_order_trotter_step!(
    psi::Vector{ComplexF64},
    terms::Vector{Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}},
    dt::Real,
    normalization::Real
)
    half_dt = dt / (2 * normalization)

    for (coeff, P) in terms
        apply_pauli_rotation!(psi, P, coeff, half_dt)
    end

    for (coeff, P) in reverse(terms)
        apply_pauli_rotation!(psi, P, coeff, half_dt)
    end

    return psi
end

function apply_second_order_trotter_adjoint_step!(
    psi::Vector{ComplexF64},
    terms::Vector{Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}},
    dt::Real,
    normalization::Real
)
    half_dt = -dt / (2 * normalization)

    for (coeff, P) in terms
        apply_pauli_rotation!(psi, P, coeff, half_dt)
    end

    for (coeff, P) in reverse(terms)
        apply_pauli_rotation!(psi, P, coeff, half_dt)
    end

    return psi
end

function apply_trotter_step!(
    psi::Vector{ComplexF64},
    terms::Vector{Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}},
    dt::Real,
    normalization::Real,
    order::Symbol
)
    if order == :first
        return apply_first_order_trotter_step!(psi, terms, dt, normalization)
    elseif order == :second
        return apply_second_order_trotter_step!(psi, terms, dt, normalization)
    else
        error("Unsupported Trotter order: $order")
    end
end

function apply_trotter_adjoint_step!(
    psi::Vector{ComplexF64},
    terms::Vector{Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}},
    dt::Real,
    normalization::Real,
    order::Symbol
)
    if order == :first
        return apply_first_order_trotter_adjoint_step!(psi, terms, dt, normalization)
    elseif order == :second
        return apply_second_order_trotter_adjoint_step!(psi, terms, dt, normalization)
    else
        error("Unsupported Trotter order: $order")
    end
end

function apply_trotter_unitary(
    x::AbstractVector{ComplexF64},
    terms::Vector{Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}},
    nsteps::Int,
    normalization::Real,
    order::Symbol;
    time::Real=pi
)
    psi = copy(Vector{ComplexF64}(x))
    dt = time / nsteps

    for _ in 1:nsteps
        apply_trotter_step!(psi, terms, dt, normalization, order)
    end

    return exp(-im * time / 2) * psi
end

function apply_trotter_unitary_adjoint(
    x::AbstractVector{ComplexF64},
    terms::Vector{Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}},
    nsteps::Int,
    normalization::Real,
    order::Symbol;
    time::Real=pi
)
    psi = copy(Vector{ComplexF64}(x))
    dt = time / nsteps

    for _ in 1:nsteps
        apply_trotter_adjoint_step!(psi, terms, dt, normalization, order)
    end

    return exp(im * time / 2) * psi
end

function reference_first_order_trotter_unitary(
    ham::Dict,
    normalization::Real,
    nqubits::Int;
    numsteps::Int=10,
    time::Real=pi
)
    U_s = sparse(I, 2^nqubits, 2^nqubits) .+ 0.0im
    Identity = sparse(I, 2^nqubits, 2^nqubits) .+ 0.0im
    id_key = "I"^nqubits
    H_terms = [(k, v) for (k, v) in ham if k != id_key]
    dt = time / (numsteps * normalization)

    for _ in 1:numsteps
        for (k, v) in H_terms
            theta = real(v) * dt
            U_s = (cos(theta) * Identity - im * sin(theta) * sparse(OP_from_string(k))) * U_s
        end
    end

    return exp(-im * time / 2) * U_s
end

function reference_trotter_unitary_by_order(
    ham::Dict{String,ComplexF64},
    normalization::Real,
    nqubits::Int,
    order::Symbol;
    numsteps::Int,
    time::Real=pi
)
    if order == :first
        return reference_first_order_trotter_unitary(
            ham,
            normalization,
            nqubits;
            numsteps=numsteps,
            time=time
        )
    elseif order == :second
        return reference_trotter_unitary(
            ham,
            normalization,
            nqubits;
            numsteps=numsteps,
            time=time
        )
    else
        error("Unsupported Trotter order: $order")
    end
end

function trotter_energy_arpack(
    meta::Dict{String,String},
    ham::Dict{String,ComplexF64},
    nsteps::Int;
    order::Symbol=:second,
    time::Real=pi,
    nev::Int=TROTTER_SPECTRAL_DEFAULT_NEV
)
    nqubits = parse(Int, meta["number of qubits"])
    norm_info = normalize_hamiltonian(meta, ham)
    U = reference_trotter_unitary_by_order(
        ham,
        norm_info.normalization,
        nqubits,
        order;
        numsteps=nsteps,
        time=time
    )

    C = sparse(U + U')
    n = size(C, 1)
    nvals = trotter_candidate_count(n; nev=nev)
    vals, _ = Arpack.eigs(C; nev=nvals, which=:LR)
    energies = cosine_value_to_energy.(real.(vals), norm_info.shift, norm_info.normalization, time)
    return minimum(real.(energies))
end

function trotter_energy_arnoldi(
    meta::Dict{String,String},
    ham::Dict{String,ComplexF64},
    nsteps::Int;
    order::Symbol=:second,
    time::Real=pi,
    nev::Int=TROTTER_SPECTRAL_DEFAULT_NEV
)
    nqubits = parse(Int, meta["number of qubits"])
    nelectrons = parse(Int, meta["number of active, occupied, single-occupancy orbitals"])
    norm_info = normalize_hamiltonian(meta, ham)
    terms = build_trotter_terms(ham, nqubits)
    n = 2^nqubits
    nvals = trotter_candidate_count(n; nev=nev)
    maxdim = trotter_krylov_dimension(n, nvals)
    mindim = min(maxdim, max(nvals, 2 * nvals))
    stateHF = construct_hf_state(nqubits, nelectrons)
    C = LinearMap{ComplexF64}(
        x -> apply_trotter_unitary(x, terms, nsteps, norm_info.normalization, order; time=time) +
             apply_trotter_unitary_adjoint(x, terms, nsteps, norm_info.normalization, order; time=time),
        n;
        ismutating=false
    )

    schur, _ = partialschur(
        C;
        v1=stateHF,
        nev=nvals,
        which=:LR,
        mindim=mindim,
        maxdim=maxdim
    )
    energies = cosine_value_to_energy.(real.(diag(schur.R)), norm_info.shift, norm_info.normalization, time)
    return minimum(real.(energies))
end

function trotter_energy(
    meta::Dict{String,String},
    ham::Dict{String,ComplexF64},
    nsteps::Int;
    method::Symbol=:arpack,
    order::Symbol=:second,
    time::Real=pi,
    nev::Int=TROTTER_SPECTRAL_DEFAULT_NEV
)
    if method == :arpack
        return trotter_energy_arpack(meta, ham, nsteps; order=order, time=time, nev=nev)
    elseif method == :arnoldi
        return trotter_energy_arnoldi(meta, ham, nsteps; order=order, time=time, nev=nev)
    else
        error("Unsupported Trotter eigensolver method: $method")
    end
end
