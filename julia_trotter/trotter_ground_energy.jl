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
#   --nsteps=N or --nsteps=N1,N2,...
#       Trotter step counts to evaluate. The default is a single step (1).
#       Give a comma-separated list of positive integers to sweep several step
#       counts; each produces one row in the printed table.
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
# reports Trotter energies for the selected formula and the step counts given
# by --nsteps (a single step by default).
# =============================================================================

using Printf

include("src/parser.jl")
include("src/trotter_spectral_energy.jl")

const DEFAULT_FILE = "He-He/He-He_2.40_hgbs-5_as-004-004_jw.dat"
const DEFAULT_NSTEPS = [1]
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

function parse_nsteps(value::AbstractString)
    entries = split(value, ',')
    steps = Int[]
    for entry in entries
        trimmed = strip(entry)
        isempty(trimmed) && continue
        n = tryparse(Int, trimmed)
        (n === nothing || n < 1) &&
            error("--nsteps must be positive integers; got '$entry'")
        push!(steps, n)
    end
    isempty(steps) && error("--nsteps requires at least one value")
    return steps
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
    nsteps = DEFAULT_NSTEPS
    filepath = DEFAULT_FILE

    for arg in args
        if arg == "--arnoldi"
            use_arnoldi = true
        elseif startswith(arg, "--nsteps=")
            nsteps = parse_nsteps(split(arg, '='; limit=2)[2])
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
        nsteps=nsteps,
    )
end

function metadata_ground_energy(meta::Dict{String,String})
    return parse(Float64, meta["smallest eigenvalue"])
end

"""
Print an aligned, comma-separated table. `columns` is a vector of header strings
and `rows` is a vector of string vectors (one per row, same length as
`columns`). Fields are comma-separated (so the output stays machine-parseable)
and space-padded to a common width per column; numeric-looking cells are
right-aligned, the rest left-aligned. Split lines on `,` and strip whitespace to
parse.
"""
function print_aligned_csv(columns::Vector{String}, rows::Vector{Vector{String}})
    ncols = length(columns)
    # A trailing comma follows every field except the last, so include it in the
    # width so the value columns still line up.
    cell(c, s) = c == ncols ? s : s * ","
    widths = [maximum(length, [cell(c, columns[c]); [cell(c, row[c]) for row in rows]]) for c in 1:ncols]
    right_align = [all(row -> occursin(r"^[-+]?[\d.eE]+$", row[c]), rows) for c in 1:ncols]

    pad(c, s) = right_align[c] ? lpad(cell(c, s), widths[c]) : rpad(cell(c, s), widths[c])
    line(cells) = rstrip(join([pad(c, cells[c]) for c in 1:ncols], " "))

    println(line(columns))
    for row in rows
        println(line(row))
    end
end

function main(
    filepath::String;
    use_arnoldi::Bool=false,
    trotter_order::Symbol=:second,
    term_ordering=:magnitude,
    safe_normalization::Bool=true,
    benchmark::Bool=false,
    nsteps=DEFAULT_NSTEPS,
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

    results = NamedTuple[]
    table_rows = Vector{String}[]
    for n in nsteps
        compute_energy = () -> trotter_energy(
            meta,
            ham,
            n;
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

        row = [
            @sprintf("%d", n),
            String(method),
            @sprintf("%.12f", energy),
            @sprintf("%.12e", energy_error),
        ]
        safe_normalization && push!(row, @sprintf("%.12f", safe_scaling_factor))
        benchmark && push!(row, @sprintf("%.6f", elapsed))
        push!(table_rows, row)
        push!(results, (
            nsteps=n,
            method=method,
            trotter_order=trotter_order,
            energy=energy,
            error=energy_error,
            safe_scaling_factor=safe_scaling_factor,
            time_seconds=elapsed,
        ))
    end

    print_aligned_csv(columns, table_rows)
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
        nsteps=opts.nsteps,
    )
end
