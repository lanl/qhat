# =============================================================================
# trotter_ground_energy.jl
#
# Read a QHAT Hamiltonian, build second-order Trotter product unitaries, and
# infer the ground state energy from the largest eigenvalues of U + U'.
#
# Usage:
#   julia --project=. trotter_ground_energy.jl [--arnoldi] [hamiltonian_file.dat]
# =============================================================================

using ArnoldiMethod, Arpack, LinearAlgebra, Printf, SparseArrays

include("src/parser.jl")
include("src/quantum_utils.jl")
include("src/hamiltonian_utils.jl")
include("src/error_bounds.jl")

const DEFAULT_FILE = "Li-Li_jw/Li-Li_2.90_hgbs-5_as-002-002_jw.dat"
const DEFAULT_NSTEPS = [1, 5, 10]
const TOTAL_TIME = π
const DEFAULT_NEV = 4

function parse_cli(args)
    use_arnoldi = false
    filepath = DEFAULT_FILE

    for arg in args
        if arg == "--arnoldi"
            use_arnoldi = true
        elseif startswith(arg, "--")
            error("Unknown option: $arg")
        else
            filepath = arg
        end
    end

    return (filepath=filepath, use_arnoldi=use_arnoldi)
end

function cosine_value_to_energy(cosine_eigenvalue::Real, shift::Real, normalization::Real, time::Real)
    theta = acos(clamp(cosine_eigenvalue / 2, -1.0, 1.0))
    return shift + normalization * (theta / time - 0.5)
end

function candidate_count(n::Int)
    return min(DEFAULT_NEV, n - 1)
end

function krylov_dimension(n::Int, nev::Int)
    return min(40, max(nev + 2, min(n - 1, 2 * nev + 1)))
end

function largest_real_eigenvalues_arpack(C, nev::Int)
    vals, _ = Arpack.eigs(C; nev=nev, which=:LR)
    return real.(vals)
end

function largest_real_eigenvalues_arnoldi(C, nev::Int, v1::AbstractVector{ComplexF64})
    n = size(C, 1)
    maxdim = krylov_dimension(n, nev)
    mindim = min(maxdim, max(nev, 2 * nev))
    schur, _ = partialschur(C; v1=v1, nev=nev, which=:LR, mindim=mindim, maxdim=maxdim)
    return real.(diag(schur.R))
end

function trotter_ground_energy(
    ham::Dict{String,ComplexF64},
    nqubits::Int,
    shift::Real,
    normalization::Real,
    nsteps::Int;
    time::Real=TOTAL_TIME,
    use_arnoldi::Bool=false,
    stateHF::Union{Nothing,Vector{ComplexF64}}=nothing
)
    U = reference_trotter_unitary(
        ham,
        normalization,
        nqubits;
        numsteps=nsteps,
        time=time
    )

    C = sparse(U + U')
    n = size(C, 1)
    nev = candidate_count(n)

    cosine_eigenvalues = if use_arnoldi
        stateHF === nothing && error("Arnoldi mode requires a Hartree-Fock starting state")
        largest_real_eigenvalues_arnoldi(C, nev, stateHF)
    else
        largest_real_eigenvalues_arpack(C, nev)
    end

    energies = cosine_value_to_energy.(cosine_eigenvalues, shift, normalization, time)
    return minimum(real.(energies))
end

function metadata_ground_energy(meta::Dict{String,String})
    return parse(Float64, meta["smallest eigenvalue"])
end

function main(filepath::String; use_arnoldi::Bool=false)
    meta, ham = parse_hamiltonian_file(filepath)

    nqubits = parse(Int, meta["number of qubits"])
    nelectrons = parse(Int, meta["number of active, occupied, single-occupancy orbitals"])
    norm_info = normalize_hamiltonian(meta, ham)
    exact_energy = metadata_ground_energy(meta)
    stateHF = use_arnoldi ? construct_hf_state(nqubits, nelectrons) : nothing
    method = use_arnoldi ? "arnoldi" : "arpack"

    println("file: ", basename(filepath))
    println(@sprintf("metadata_ground_energy: %.12f", exact_energy))
    println("nsteps,method,trotter_ground_energy,error,time_seconds")

    for nsteps in DEFAULT_NSTEPS
        energy_ref = Ref{Float64}()
        elapsed = @elapsed begin
            energy_ref[] = trotter_ground_energy(
                ham,
                nqubits,
                norm_info.shift,
                norm_info.normalization,
                nsteps;
                time=TOTAL_TIME,
                use_arnoldi=use_arnoldi,
                stateHF=stateHF
            )
        end
        energy = energy_ref[]
        println(@sprintf("%d,%s,%.12f,%.12e,%.6f",
                         nsteps, method, energy, energy - exact_energy, elapsed))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_cli(ARGS)
    main(opts.filepath; use_arnoldi=opts.use_arnoldi)
end
