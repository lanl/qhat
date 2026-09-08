# trotter_spectral_energy.jl
#
# Shared helpers for extracting Trotter-derived ground energies from the
# largest eigenvalues of U + U'.

using ArnoldiMethod, Arpack, LinearAlgebra, LinearMaps, SparseArrays

include("hamiltonian_utils.jl")
include("error_bounds.jl")
include("statevector_simulators.jl")

const TROTTER_SPECTRAL_DEFAULT_NEV = 4
const TROTTER_ARNOLDI_BASE_TOL = 1e-12

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

function trotter_single_step_time(time::Real, nsteps::Int)
    nsteps >= 1 || throw(ArgumentError("nsteps must be positive"))
    return time / nsteps
end

function trotter_arnoldi_tolerance(nsteps::Int)
    nsteps >= 1 || throw(ArgumentError("nsteps must be positive"))
    # Eigenvalue separations in U + U' shrink quadratically with the
    # single-step evolution time. Tighten the residual tolerance by the same
    # factor, down to the useful limit of Float64 arithmetic.
    return max(eps(Float64), TROTTER_ARNOLDI_BASE_TOL / Float64(nsteps)^2)
end

function trotter_safety_parameters(
    ham::Dict{String,ComplexF64},
    normalization::Real,
    nqubits::Int,
    nsteps::Int,
    order::Symbol;
    term_ordering=:magnitude,
)
    ε = if order == :first
        only(first_order_commutator_error_bounds(
            ham,
            normalization,
            nqubits;
            nsteps_list=[nsteps],
            time=π,
            term_ordering=term_ordering,
        ))
    elseif order == :second
        bounds, _ = commutator_error_bounds(
            ham,
            normalization,
            nqubits;
            nsteps_list=[nsteps],
            time=π,
            term_ordering=term_ordering,
        )
        only(bounds)
    else
        return (error_bound=0.0, shift=0.0, scale=1.0)
    end

    0 <= ε < sqrt(2) || error(
        "The certified $order-order Trotter bound must be in [0, sqrt(2)); got ε=$ε",
    )
    b = 2asin(ε / 2)
    s = (π - 2b) / π
    return (error_bound=ε, shift=b, scale=s)
end

function safe_cosine_value_to_energy(
    cosine_eigenvalue::Real,
    shift::Real,
    normalization::Real,
    step_time::Real,
    safety_scale::Real,
    safety_phase_shift::Real,
)
    safe_phase = acos(clamp(cosine_eigenvalue / 2, -1.0, 1.0))
    original_phase = (safe_phase - safety_phase_shift) / safety_scale
    return shift + normalization * (original_phase / step_time - 0.5)
end

function format_trotter_energy_result(energy::Real, safety, return_details::Bool)
    return return_details ? (
        energy=Float64(real(energy)),
        safe_scaling_factor=Float64(safety.scale),
    ) : Float64(real(energy))
end

function build_trotter_terms(
    ham::Dict{String,ComplexF64},
    nqubits::Int;
    term_ordering=:magnitude,
)
    terms = Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}[]

    for (k, v) in ordered_hamiltonian_terms(
        ham, nqubits; term_ordering=term_ordering
    )
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
    time::Real=pi,
    term_ordering=:magnitude,
)
    U_s = sparse(I, 2^nqubits, 2^nqubits) .+ 0.0im
    Identity = sparse(I, 2^nqubits, 2^nqubits) .+ 0.0im
    H_terms = ordered_hamiltonian_terms(
        ham, nqubits; term_ordering=term_ordering
    )
    dt = time / (numsteps * normalization)

    for _ in 1:numsteps
        for (k, v) in H_terms
            theta = real(v) * dt
            U_s = (cos(theta) * Identity - im * sin(theta) * sparse(OP_from_string(k))) * U_s
        end
    end

    return exp(-im * time / 2) * U_s
end

# Apply one symmetric Trotter step of the given even `order` to U_s (left-multiply).
# order==2 is Strang; higher even orders use the Suzuki fractal recursion
#   S_{2k}(dt) = S_{2k-2}(p·dt)² S_{2k-2}((1-4p)·dt) S_{2k-2}(p·dt)²,  p = 1/(4 - 4^(1/(2k-1)))
# `ops` are precomputed (coeff, sparse_pauli) pairs; `dt` is already normalized.
function apply_symmetric_trotter_step(
    U_s,
    Identity,
    ops::Vector{Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}},
    dt::Real,
    order::Int
)
    if order == 2
        for (coeff, P) in ops
            θ = coeff * dt / 2
            U_s = (cos(θ) * Identity - im * sin(θ) * P) * U_s
        end
        for (coeff, P) in reverse(ops)
            θ = coeff * dt / 2
            U_s = (cos(θ) * Identity - im * sin(θ) * P) * U_s
        end
        return U_s
    else
        p = 1 / (4 - 4^(1 / (order - 1)))
        for δ in (p * dt, p * dt, (1 - 4p) * dt, p * dt, p * dt)
            U_s = apply_symmetric_trotter_step(U_s, Identity, ops, δ, order - 2)
        end
        return U_s
    end
end

# Reference symmetric (even-order) Trotter unitary built from the Suzuki recursion.
function reference_symmetric_trotter_unitary(
    ham::Dict,
    normalization::Real,
    nqubits::Int,
    order::Int;
    numsteps::Int,
    time::Real=pi,
    term_ordering=:magnitude,
)
    dim = 2^nqubits
    U_s = sparse(I, dim, dim) .+ 0.0im
    Identity = sparse(I, dim, dim) .+ 0.0im
    ops = Tuple{Float64,SparseMatrixCSC{ComplexF64,Int}}[
        (real(v), sparse(OP_from_string(k)))
        for (k, v) in ordered_hamiltonian_terms(
            ham, nqubits; term_ordering=term_ordering
        )
    ]
    dt = time / (numsteps * normalization)

    for _ in 1:numsteps
        U_s = apply_symmetric_trotter_step(U_s, Identity, ops, dt, order)
    end

    return exp(-im * time / 2) * U_s
end

function reference_trotter_unitary_by_order(
    ham::Dict{String,ComplexF64},
    normalization::Real,
    nqubits::Int,
    order::Symbol;
    numsteps::Int,
    time::Real=pi,
    term_ordering=:magnitude,
)
    if order == :first
        return reference_first_order_trotter_unitary(
            ham,
            normalization,
            nqubits;
            numsteps=numsteps,
            time=time,
            term_ordering=term_ordering,
        )
    elseif order == :second
        return reference_trotter_unitary(
            ham,
            normalization,
            nqubits;
            numsteps=numsteps,
            time=time,
            term_ordering=term_ordering,
        )
    elseif order == :fourth
        return reference_symmetric_trotter_unitary(
            ham,
            normalization,
            nqubits,
            4;
            numsteps=numsteps,
            time=time,
            term_ordering=term_ordering,
        )
    elseif order == :sixth
        return reference_symmetric_trotter_unitary(
            ham,
            normalization,
            nqubits,
            6;
            numsteps=numsteps,
            time=time,
            term_ordering=term_ordering,
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
    nev::Int=TROTTER_SPECTRAL_DEFAULT_NEV,
    term_ordering=:magnitude,
    safe_normalization::Bool=true,
    return_details::Bool=false,
)
    nqubits = parse(Int, meta["number of qubits"])
    norm_info = normalize_hamiltonian(meta, ham)
    step_time = trotter_single_step_time(time, nsteps)
    safety = safe_normalization ?
             trotter_safety_parameters(
                 ham,
                 norm_info.normalization,
                 nqubits,
                 nsteps,
                 order;
                 term_ordering=term_ordering,
             ) :
             (error_bound=0.0, shift=0.0, scale=1.0)
    safe_step_time = safety.scale * step_time
    safety_phase_shift = safety.shift * time / (π * nsteps)
    U = reference_trotter_unitary_by_order(
        ham,
        norm_info.normalization,
        nqubits,
        order;
        numsteps=1,
        time=safe_step_time,
        term_ordering=term_ordering,
    )
    U .*= exp(-im * safety_phase_shift)

    C = sparse(U + U')
    n = size(C, 1)
    nvals = trotter_candidate_count(n; nev=nev)
    # tol=0 requests machine precision from ARPACK. Its default was already
    # tighter than the scale-dependent ArnoldiMethod tolerance below.
    vals, _ = Arpack.eigs(C; nev=nvals, which=:LR, tol=0.0)
    energies = safe_cosine_value_to_energy.(
        real.(vals),
        norm_info.shift,
        norm_info.normalization,
        step_time,
        safety.scale,
        safety_phase_shift,
    )
    energy = minimum(real.(energies))
    return format_trotter_energy_result(energy, safety, return_details)
end

function trotter_energy_arnoldi(
    meta::Dict{String,String},
    ham::Dict{String,ComplexF64},
    nsteps::Int;
    order::Symbol=:second,
    time::Real=pi,
    nev::Int=TROTTER_SPECTRAL_DEFAULT_NEV,
    term_ordering=:magnitude,
    safe_normalization::Bool=true,
    return_details::Bool=false,
)
    nqubits = parse(Int, meta["number of qubits"])
    nelectrons = parse(Int, meta["number of active, occupied, single-occupancy orbitals"])
    norm_info = normalize_hamiltonian(meta, ham)
    step_time = trotter_single_step_time(time, nsteps)
    solver_tolerance = trotter_arnoldi_tolerance(nsteps)
    safety = safe_normalization ?
             trotter_safety_parameters(
                 ham,
                 norm_info.normalization,
                 nqubits,
                 nsteps,
                 order;
                 term_ordering=term_ordering,
             ) :
             (error_bound=0.0, shift=0.0, scale=1.0)
    safe_step_time = safety.scale * step_time
    safety_phase_shift = safety.shift * time / (π * nsteps)
    safe_phase = exp(-im * safety_phase_shift)
    terms = build_trotter_terms(ham, nqubits; term_ordering=term_ordering)
    n = 2^nqubits
    nvals = trotter_candidate_count(n; nev=nev)
    maxdim = trotter_krylov_dimension(n, nvals)
    mindim = min(maxdim, max(nvals, 2 * nvals))
    stateHF = construct_hf_state(nqubits, nelectrons)
    C = LinearMap{ComplexF64}(
        x -> safe_phase * apply_trotter_unitary(
                 x, terms, 1, norm_info.normalization, order; time=safe_step_time
             ) + conj(safe_phase) * apply_trotter_unitary_adjoint(
                 x, terms, 1, norm_info.normalization, order; time=safe_step_time
             ),
        n;
        ismutating=false
    )

    schur, _ = partialschur(
        C;
        v1=stateHF,
        nev=nvals,
        which=:LR,
        mindim=mindim,
        maxdim=maxdim,
        tol=solver_tolerance
    )
    energies = safe_cosine_value_to_energy.(
        real.(diag(schur.R)),
        norm_info.shift,
        norm_info.normalization,
        step_time,
        safety.scale,
        safety_phase_shift,
    )
    energy = minimum(real.(energies))
    return format_trotter_energy_result(energy, safety, return_details)
end

function trotter_energy(
    meta::Dict{String,String},
    ham::Dict{String,ComplexF64},
    nsteps::Int;
    method::Symbol=:arpack,
    order::Symbol=:second,
    time::Real=pi,
    nev::Int=TROTTER_SPECTRAL_DEFAULT_NEV,
    term_ordering=:magnitude,
    safe_normalization::Bool=true,
    return_details::Bool=false,
)
    if method == :arpack
        return trotter_energy_arpack(
            meta, ham, nsteps;
            order=order,
            time=time,
            nev=nev,
            term_ordering=term_ordering,
            safe_normalization=safe_normalization,
            return_details=return_details,
        )
    elseif method == :arnoldi
        return trotter_energy_arnoldi(
            meta, ham, nsteps;
            order=order,
            time=time,
            nev=nev,
            term_ordering=term_ordering,
            safe_normalization=safe_normalization,
            return_details=return_details,
        )
    else
        error("Unsupported Trotter eigensolver method: $method")
    end
end
