# =============================================================================
# trotter_ground_energy.jl
#
# Read a QHAT Hamiltonian, build second-order Trotter product unitaries, and
# infer the ground state energy from the largest eigenvalues of U + U'.
#
# Usage:
#   julia --project=. trotter_ground_energy.jl [options] [hamiltonian_file.dat]
#
# Options:
#
#   --arnoldi
#       Use the matrix-free Arnoldi eigensolver. The default is Arpack, which
#       constructs the full Trotter unitary. Arnoldi applies the unitary and
#       its adjoint directly to state vectors and is preferable when the full
#       matrix is too large.
#
#   --trotter-order=first|second
#       Select the product formula. First order applies one full rotation for
#       each Pauli term. Second order (the default) uses the symmetric Strang
#       sequence of forward and reverse half rotations. The matching first- or
#       second-order commutator bound is used by safe normalization.
#
#   --term-order=magnitude
#       Apply nonidentity Pauli terms in decreasing coefficient magnitude.
#       This is the default. Equal-magnitude terms are ordered by their Pauli
#       strings so the result is deterministic.
#
#   --term-order=increasing-magnitude
#       Apply nonidentity Pauli terms in increasing coefficient magnitude,
#       again breaking ties by Pauli string. When this reverses the decreasing
#       order, the two Trotter matrices can differ while their `U + U'`
#       spectra, and therefore their extracted energies, remain identical.
#       This occurs exactly for reversed second-order Strang products and can
#       also occur for first-order products of transpose-symmetric terms.
#
#   --term-order=lexicographic
#       Apply terms in lexicographic Pauli-string order, independently of
#       their coefficients.
#
#   --term-order=dict
#       Retain the iteration order of the Julia Hamiltonian Dict. This option
#       reproduces the old behavior, but the order should not be treated as a
#       file-order guarantee.
#
#   --term-order=PAULI1,PAULI2,...
#       Supply an explicit chronological application order. The list must
#       contain every nonidentity Pauli string exactly once. For example,
#       --term-order=XX,ZI,IZ applies XX first, then ZI, then IZ.
#
#   --safe-normalization=true|false
#       With true (the default), compute the certified commutator error ε and
#       use b=2asin(ε/2), s=(π-2b)/π to evolve the safely transformed
#       Hamiltonian H_safe=sH+bI. The transform is undone when converting the
#       measured phase to an energy in Hartrees. Each result row reports s as
#       `safe_scaling_factor`. With false, use b=0 and s=1, skip this
#       commutator-bound safety transform, and omit the scaling-factor column.
#
#   --no-safe-normalization
#       Shorthand for --safe-normalization=false.
#
#   --benchmark[=true|false]
#       Benchmark each complete Trotter diagonalization. This is false by
#       default. When enabled, the script measures each `trotter_energy` call,
#       adds `time_seconds` to the printed results, and returns the measured
#       time in each result record. The first measurement can include Julia
#       compilation time.
#
# If no Hamiltonian file is given, DEFAULT_FILE below is used. The script
# reports Trotter energies for the selected formula and the step counts in
# DEFAULT_NSTEPS.
# =============================================================================

using Printf

include("src/parser.jl")
include("src/trotter_spectral_energy.jl")

const DEFAULT_FILE = "He-He/He-He_2.40_hgbs-5_as-004-004_jw.dat"
const DEFAULT_NSTEPS = [1, 5, 10]
const TOTAL_TIME = π

function parse_boolean_option(value::AbstractString, option::AbstractString)
    normalized = lowercase(value)
    normalized == "true" && return true
    normalized == "false" && return false
    error("$option must be true or false; got '$value'")
end

function parse_term_ordering(value::AbstractString)
    normalized = lowercase(value)
    builtin = Dict(
        "magnitude" => :magnitude,
        "decreasing-magnitude" => :magnitude,
        "increasing-magnitude" => :increasing_magnitude,
        "lexicographic" => :lexicographic,
        "dict" => :dict,
    )
    haskey(builtin, normalized) && return builtin[normalized]

    pauli_strings = String.(split(value, ','))
    all(pauli -> !isempty(pauli) && all(in("IXYZ"), pauli), pauli_strings) &&
        return pauli_strings
    error("Unknown term order '$value'")
end

function parse_trotter_order(value::AbstractString)
    normalized = lowercase(value)
    normalized in ("first", "1") && return :first
    normalized in ("second", "2") && return :second
    error("--trotter-order must be first or second; got '$value'")
end

function parse_cli(args)
    use_arnoldi = false
    trotter_order = :second
    term_ordering = :magnitude
    safe_normalization = true
    benchmark = false
    filepath = DEFAULT_FILE

    for arg in args
        if arg == "--arnoldi"
            use_arnoldi = true
        elseif startswith(arg, "--trotter-order=")
            trotter_order = parse_trotter_order(split(arg, '='; limit=2)[2])
        elseif startswith(arg, "--term-order=")
            term_ordering = parse_term_ordering(split(arg, '='; limit=2)[2])
        elseif startswith(arg, "--safe-normalization=")
            safe_normalization = parse_boolean_option(
                split(arg, '='; limit=2)[2], "--safe-normalization"
            )
        elseif arg == "--no-safe-normalization"
            safe_normalization = false
        elseif arg == "--benchmark"
            benchmark = true
        elseif startswith(arg, "--benchmark=")
            benchmark = parse_boolean_option(
                split(arg, '='; limit=2)[2], "--benchmark"
            )
        elseif startswith(arg, "--")
            error("Unknown option: $arg")
        else
            filepath = arg
        end
    end

    return (
        filepath=filepath,
        use_arnoldi=use_arnoldi,
        trotter_order=trotter_order,
        term_ordering=term_ordering,
        safe_normalization=safe_normalization,
        benchmark=benchmark,
    )
end

function metadata_ground_energy(meta::Dict{String,String})
    return parse(Float64, meta["smallest eigenvalue"])
end

function main(
    filepath::String;
    use_arnoldi::Bool=false,
    trotter_order::Symbol=:second,
    term_ordering=:magnitude,
    safe_normalization::Bool=true,
    benchmark::Bool=false,
)
    meta, ham = parse_hamiltonian_file(filepath)

    exact_energy = metadata_ground_energy(meta)
    method = use_arnoldi ? :arnoldi : :arpack
    trotter_order in (:first, :second) || throw(ArgumentError(
        "trotter_order must be :first or :second",
    ))

    println("file: ", basename(filepath))
    println(@sprintf("metadata_ground_energy: %.12f", exact_energy))
    println("trotter_order: ", trotter_order)
    println("term_ordering: ", term_ordering)
    println("safe_normalization: ", safe_normalization)
    println("benchmark: ", benchmark)
    columns = ["nsteps", "method", "trotter_ground_energy", "error"]
    safe_normalization && push!(columns, "safe_scaling_factor")
    benchmark && push!(columns, "time_seconds")
    println(join(columns, ','))

    results = NamedTuple[]
    for nsteps in DEFAULT_NSTEPS
        compute_energy = () -> trotter_energy(
            meta,
            ham,
            nsteps;
            method=method,
            order=trotter_order,
            time=TOTAL_TIME,
            term_ordering=term_ordering,
            safe_normalization=safe_normalization,
            return_details=true,
        )

        elapsed = nothing
        if benchmark
            calculation_ref = Ref{Any}()
            elapsed = @elapsed calculation_ref[] = compute_energy()
            calculation = calculation_ref[]
        else
            calculation = compute_energy()
        end
        energy = calculation.energy
        safe_scaling_factor = calculation.safe_scaling_factor
        energy_error = energy - exact_energy

        row = @sprintf("%d,%s,%.12f,%.12e",
                       nsteps, String(method), energy, energy_error)
        safe_normalization && (row *= @sprintf(",%.12f", safe_scaling_factor))
        benchmark && (row *= @sprintf(",%.6f", elapsed))
        println(row)
        push!(results, (
            nsteps=nsteps,
            method=method,
            trotter_order=trotter_order,
            energy=energy,
            error=energy_error,
            safe_scaling_factor=safe_scaling_factor,
            time_seconds=elapsed,
        ))
    end

    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_cli(ARGS)
    main(
        opts.filepath;
        use_arnoldi=opts.use_arnoldi,
        trotter_order=opts.trotter_order,
        term_ordering=opts.term_ordering,
        safe_normalization=opts.safe_normalization,
        benchmark=opts.benchmark,
    )
end
