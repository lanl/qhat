# experiment_utils.jl
#
# Common infrastructure for parallel experiment execution.
# Provides templates for SLURM array jobs and result merging.

using JLD2, DataFrames, Printf

"""
    parallel_worker_template(process_fn, paths_file, outdir; error_fields=nothing)

Generic worker for SLURM array jobs.

Reads a file containing paths (one per line), processes the path corresponding
to SLURM_ARRAY_TASK_ID using the provided process function, and saves results.

# Arguments
- `process_fn`: Function that takes a filepath and returns a NamedTuple of results
- `paths_file`: Path to file with input paths (one per line)
- `outdir`: Directory to write partial results
- `error_fields`: Optional NamedTuple specifying default error record structure

# Environment Variables
- `SLURM_ARRAY_TASK_ID`: Task ID (1-indexed)

# Output
Writes `part-NNNNNN.jld2` containing a DataFrame with one row
"""
function parallel_worker_template(
    process_fn::Function,
    paths_file::String,
    outdir::String;
    error_fields::Union{Nothing, NamedTuple}=nothing
)
    mkpath(outdir)
    paths = filter(!isempty, strip.(readlines(paths_file)))

    task_id = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
    @assert 1 <= task_id <= length(paths) "SLURM_ARRAY_TASK_ID=$task_id out of bounds (1..$(length(paths)))"

    filepath = paths[task_id]
    start_time = time()

    row = try
        process_fn(filepath)
    catch err
        # Return error record
        if error_fields !== nothing
            merge(error_fields, (
                filepath=String(filepath),
                ok=false,
                error_message=sprint(showerror, err, catch_backtrace()),
            ))
        else
            (
                filepath=String(filepath),
                ok=false,
                error_message=sprint(showerror, err, catch_backtrace()),
            )
        end
    end

    elapsed_time = time() - start_time

    df = DataFrame([row])
    partfile = joinpath(outdir, @sprintf("part-%06d.jld2", task_id))
    JLD2.jldsave(partfile; db=df)

    ok_status = haskey(row, :ok) ? row.ok : true
    @info "Wrote part" partfile ok=ok_status filepath elapsed_time=@sprintf("%.2f s", elapsed_time)
end

"""
    merge_results(outdir, outfile)

Merge all part-*.jld2 files from a parallel run into single output file.

# Arguments
- `outdir`: Directory containing part-NNNNNN.jld2 files
- `outfile`: Output path for merged results

# Output
Writes JLD2 file with key "db" containing merged DataFrame
"""
function merge_results(outdir::String, outfile::String)
    files = sort(filter(f -> endswith(f, ".jld2"), readdir(outdir; join=true)))
    @assert !isempty(files) "No .jld2 part files found in $outdir"

    dfs = DataFrame[]
    for file in files
        push!(dfs, JLD2.load(file, "db"))
    end

    db = vcat(dfs...; cols=:union)
    JLD2.jldsave(outfile; db=db)
    @info "Merged" n_parts=length(files) n_rows=nrow(db) outfile
end

"""
    load_experiment_database(filename)

Load experiment results from JLD2 file.

# Arguments
- `filename`: Path to .jld2 file

# Returns
DataFrame containing results
"""
load_experiment_database(filename::String) = JLD2.load(filename, "db")

"""
    save_experiment_database(filename, df)

Save experiment results to JLD2 file.

# Arguments
- `filename`: Output path
- `df`: DataFrame to save
"""
function save_experiment_database(filename::String, df::DataFrame)
    JLD2.jldsave(filename; db=df)
end

"""
    input_exists(df, inputs)

Check if a row with given input values exists in DataFrame.

# Arguments
- `df`: DataFrame to search
- `inputs`: NamedTuple of (column => value) pairs to match

# Returns
`true` if a matching row exists, `false` otherwise
"""
function input_exists(df::DataFrame, inputs::NamedTuple)
    for row in eachrow(df)
        if all(getproperty(row, k) == v for (k, v) in pairs(inputs))
            return true
        end
    end
    return false
end

"""
    delete_data(df, inputs)

Remove rows matching given input values from DataFrame (in-place).

# Arguments
- `df`: DataFrame to modify
- `inputs`: NamedTuple of (column => value) pairs to match

# Returns
Modified DataFrame (same reference)
"""
function delete_data(df::DataFrame, inputs::NamedTuple)
    filter!(row -> !all(getproperty(row, k) == v for (k, v) in pairs(inputs)), df)
end
