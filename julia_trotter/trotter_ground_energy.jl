# =============================================================================
# trotter_ground_energy.jl
#
# Read a QHAT Hamiltonian, build second-order Trotter product unitaries, and
# infer the ground state energy from the largest eigenvalues of U + U'.
#
# Usage:
#   julia --project=. trotter_ground_energy.jl [--arnoldi] [hamiltonian_file.dat]
# =============================================================================

using Printf

include("src/parser.jl")
include("src/trotter_spectral_energy.jl")

const DEFAULT_FILE = "He-He/He-He_2.40_hgbs-5_as-004-004_jw.dat"
const DEFAULT_NSTEPS = [1, 5, 10]
const TOTAL_TIME = π

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

function metadata_ground_energy(meta::Dict{String,String})
    return parse(Float64, meta["smallest eigenvalue"])
end

function main(filepath::String; use_arnoldi::Bool=false)
    meta, ham = parse_hamiltonian_file(filepath)

    exact_energy = metadata_ground_energy(meta)
    method = use_arnoldi ? :arnoldi : :arpack

    println("file: ", basename(filepath))
    println(@sprintf("metadata_ground_energy: %.12f", exact_energy))
    println("nsteps,method,trotter_ground_energy,error,time_seconds")

    for nsteps in DEFAULT_NSTEPS
        energy_ref = Ref{Float64}()
        elapsed = @elapsed begin
            energy_ref[] = trotter_energy(
                meta,
                ham,
                nsteps;
                method=method,
                order=:second,
                time=TOTAL_TIME,
            )
        end
        energy = energy_ref[]
        println(@sprintf("%d,%s,%.12f,%.12e,%.6f",
                         nsteps, String(method), energy, energy - exact_energy, elapsed))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_cli(ARGS)
    main(opts.filepath; use_arnoldi=opts.use_arnoldi)
end
