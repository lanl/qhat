# =============================================================================
# bond_energy_curve.jl
#
# Scan Li-Li bond-length Hamiltonians for a fixed active space and encoding,
# print exact and Trotter ground-state energies, and save an overlay plot.
#
# Usage:
#   julia --project=. bond_energy_curve.jl [options]
#
# Options:
#   --basis BASIS        Basis directory to scan (default: hgbs-5)
#   --electrons N        Active occupied single-occupancy orbitals (default: 4)
#   --vacant N           Active vacant single-occupancy orbitals (default: 4)
#   --encoding jw|bk     Fermion-to-qubit encoding (default: jw)
#   --out PATH           Plot output path
#   --arnoldi            Use matrix-free Arnoldi instead of full-unitary Arpack
# =============================================================================

using Arpack, LinearAlgebra, Printf, Plots

include("src/parser.jl")
include("src/trotter_spectral_energy.jl")

const DEFAULT_ROOT = "Li-Li_bond"
const DEFAULT_BASIS = "hgbs-5"
const DEFAULT_ELECTRONS = 4
const DEFAULT_VACANT = 4
const DEFAULT_ENCODING = "jw"
const TROTTER_STEPS = [1, 5]
const TOTAL_TIME = π

function usage()
    println("""
Usage:
  julia --project=. bond_energy_curve.jl [options]

Options:
  --basis BASIS        Basis directory to scan (default: $(DEFAULT_BASIS))
  --electrons N        Active occupied single-occupancy orbitals (default: $(DEFAULT_ELECTRONS))
  --vacant N           Active vacant single-occupancy orbitals (default: $(DEFAULT_VACANT))
  --encoding jw|bk     Fermion-to-qubit encoding (default: $(DEFAULT_ENCODING))
  --out PATH           Plot output path
  --arnoldi            Use matrix-free Arnoldi instead of full-unitary Arpack
  --help               Show this message
""")
end

function require_value(args, i, flag)
    i < length(args) || error("Missing value after $flag")
    return args[i + 1]
end

function parse_cli(args)
    basis = DEFAULT_BASIS
    electrons = DEFAULT_ELECTRONS
    vacant = DEFAULT_VACANT
    encoding = DEFAULT_ENCODING
    outpath = nothing
    method = :arpack

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help"
            usage()
            exit(0)
        elseif arg == "--basis"
            basis = require_value(args, i, arg)
            i += 2
        elseif arg == "--electrons"
            electrons = parse(Int, require_value(args, i, arg))
            i += 2
        elseif arg == "--vacant"
            vacant = parse(Int, require_value(args, i, arg))
            i += 2
        elseif arg == "--encoding"
            encoding = lowercase(require_value(args, i, arg))
            i += 2
        elseif arg == "--out"
            outpath = require_value(args, i, arg)
            i += 2
        elseif arg == "--arnoldi"
            method = :arnoldi
            i += 1
        else
            error("Unknown option: $arg")
        end
    end

    encoding in ("jw", "bk") || error("Encoding must be 'jw' or 'bk', got '$encoding'")
    electrons >= 0 || error("--electrons must be nonnegative")
    vacant >= 0 || error("--vacant must be nonnegative")

    return (
        basis=basis,
        electrons=electrons,
        vacant=vacant,
        encoding=encoding,
        outpath=outpath,
        method=method
    )
end

function default_outpath(opts)
    return @sprintf(
        "bond_energy_curve_%s_as-%03d-%03d_%s.png",
        opts.basis,
        opts.electrons,
        opts.vacant,
        opts.encoding
    )
end

function expected_filename(bond::AbstractString, opts)
    return @sprintf(
        "Li-Li_%s_%s_as-%03d-%03d_%s.dat",
        bond,
        opts.basis,
        opts.electrons,
        opts.vacant,
        opts.encoding
    )
end

function metadata_energy(meta::Dict{String,String})
    for key in ("smallest eigenvalue", "smallest eigenvalue (Hartrees)")
        haskey(meta, key) && return parse(Float64, meta[key])
    end
    return nothing
end

function compute_ground_energy(meta::Dict{String,String}, ham::Dict{String,ComplexF64})
    nqubits = parse(Int, meta["number of qubits"])
    H = build_sparse_hamiltonian(ham, nqubits)
    vals, _ = Arpack.eigs(Hermitian(H); nev=1, which=:SR, maxiter=5000)
    return real(vals[1])
end

function trotter_order()
    return :second
end

function parse_energies(path::AbstractString, opts)
    meta, ham = parse_hamiltonian_file(path)

    exact_ref = Ref{Float64}()
    exact_seconds = @elapsed begin
        energy = metadata_energy(meta)
        exact_ref[] = energy === nothing ? compute_ground_energy(meta, ham) : energy
    end

    trotter = Dict{Int,Float64}()
    trotter_seconds = Dict{Int,Float64}()
    for nsteps in TROTTER_STEPS
        energy_ref = Ref{Float64}()
        elapsed = @elapsed begin
            energy_ref[] = trotter_energy(
                meta,
                ham,
                nsteps;
                method=opts.method,
                order=trotter_order(),
                time=TOTAL_TIME
            )
        end
        trotter[nsteps] = energy_ref[]
        trotter_seconds[nsteps] = elapsed
    end

    return (
        exact=exact_ref[],
        exact_seconds=exact_seconds,
        trotter=trotter,
        trotter_seconds=trotter_seconds
    )
end

function print_table_header()
    println("bond_length,exact_energy_hartree,trotter_1_step_hartree,trotter_5_step_hartree,exact_seconds,trotter_1_seconds,trotter_5_seconds,total_seconds,file")
    flush(stdout)
end

function print_table_row(point)
    println(@sprintf(
        "%.2f,%.12f,%.12f,%.12f,%.6f,%.6f,%.6f,%.6f,%s",
        point.bond_length,
        point.exact_energy,
        point.trotter_1,
        point.trotter_5,
        point.exact_seconds,
        point.trotter_1_seconds,
        point.trotter_5_seconds,
        point.total_seconds,
        point.file
    ))
    flush(stdout)
end

function collect_curve_points(opts; stream_rows::Bool=false)
    isdir(DEFAULT_ROOT) || error("Hamiltonian root not found: $(DEFAULT_ROOT)")

    points = NamedTuple[]
    for bond in sort(readdir(DEFAULT_ROOT), by=x -> parse(Float64, x))
        basis_dir = joinpath(DEFAULT_ROOT, bond, opts.basis)
        isdir(basis_dir) || continue

        filename = expected_filename(bond, opts)
        path = joinpath(basis_dir, filename)
        isfile(path) || continue

        energies = parse_energies(path, opts)

        push!(points, (
            bond_length=parse(Float64, bond),
            exact_energy=energies.exact,
            trotter_1=energies.trotter[1],
            trotter_5=energies.trotter[5],
            exact_seconds=energies.exact_seconds,
            trotter_1_seconds=energies.trotter_seconds[1],
            trotter_5_seconds=energies.trotter_seconds[5],
            total_seconds=energies.exact_seconds + energies.trotter_seconds[1] + energies.trotter_seconds[5],
            file=filename
        ))

        stream_rows && print_table_row(points[end])
    end

    isempty(points) && error(
        "No matching Hamiltonians found for basis=$(opts.basis), " *
        "as-$(lpad(opts.electrons, 3, '0'))-$(lpad(opts.vacant, 3, '0')), " *
        "encoding=$(opts.encoding)"
    )

    return points
end

function print_table(points)
    print_table_header()
    for point in points
        print_table_row(point)
    end
end

function save_plot(points, opts)
    outpath = opts.outpath === nothing ? default_outpath(opts) : opts.outpath
    xs = [point.bond_length for point in points]
    exact = [point.exact_energy for point in points]
    trotter_1 = [point.trotter_1 for point in points]
    trotter_5 = [point.trotter_5 for point in points]

    plot(
        xs,
        exact;
        label="exact",
        marker=:circle,
        linewidth=2,
        xlabel="Li-Li bond length",
        ylabel="Ground state energy (Ha)",
        title=@sprintf(
            "Li-Li %s as-%03d-%03d %s",
            opts.basis,
            opts.electrons,
            opts.vacant,
            uppercase(opts.encoding)
        ),
        legend=:best
    )
    plot!(xs, trotter_1; label="Trotter 1 step", marker=:square, linewidth=2)
    plot!(xs, trotter_5; label="Trotter 5 steps", marker=:diamond, linewidth=2)
    savefig(outpath)
    return outpath
end

function main(args)
    opts = parse_cli(args)
    print_table_header()
    points = collect_curve_points(opts; stream_rows=true)

    outpath = save_plot(points, opts)
    println("plot: ", outpath)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
