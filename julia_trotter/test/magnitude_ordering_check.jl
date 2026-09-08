# Run with: julia --project=. test/magnitude_ordering_check.jl
# Compare every bundled 4- and 6-qubit Hamiltonian using full dense spectra.
# Safety-off comparisons assert the observed agreement for these datasets;
# this is not a claim about arbitrary Hamiltonians or term permutations.
using LinearAlgebra, Printf
include(joinpath(@__DIR__, "..", "src", "parser.jl"))
include(joinpath(@__DIR__, "..", "src", "trotter_spectral_energy.jl"))

function check_magnitude_orderings()
    root = dirname(@__DIR__)
    files = sort([joinpath(dir, file) for (dir, _, files) in walkdir(root)
                  for file in files if endswith(file, ".dat") &&
                  occursin(r"as-002-00[24]_", file)])
    println("file,qubits,formula,nsteps,safe_normalization,exact_reverse,tied_noncommuting_pairs,energy_decreasing,energy_increasing,error_decreasing,error_increasing,energy_difference,cosine_spectrum_difference,s_decreasing,s_increasing")
    for path in files
        meta, ham = parse_hamiltonian_file(path)
        n = parse(Int, meta["number of qubits"])
        info = normalize_hamiltonian(meta, ham)
        down = ordered_hamiltonian_terms(ham, n; term_ordering=:magnitude)
        up = ordered_hamiltonian_terms(ham, n; term_ordering=:increasing_magnitude)
        exact_reverse = first.(up) == reverse(first.(down))
        # Independent character-based check: an odd number of distinct,
        # nonidentity local Paulis means anticommutation.
        tied_noncommuting = count(
            abs(down[j][2]) == abs(down[k][2]) &&
            isodd(count(((a, b),) -> a != 'I' && b != 'I' && a != b,
                        zip(down[j][1], down[k][1])))
            for j in eachindex(down) for k in j+1:length(down)
        )
        reference_energy = eigmin(Hermitian(Matrix(build_sparse_hamiltonian(ham, n))))
        for order in (:first, :second), nsteps in (1, 5, 10), safe in (false, true)
            energies, spectra, scales = Float64[], Vector{Float64}[], Float64[]
            for ordering in (:magnitude, :increasing_magnitude)
                safety = safe ? trotter_safety_parameters(
                    ham, info.normalization, n, nsteps, order;
                    term_ordering=ordering,
                ) : (shift=0.0, scale=1.0)
                step_time = π / nsteps
                phase_shift = safety.shift / nsteps
                U = Matrix(reference_trotter_unitary_by_order(
                    ham, info.normalization, n, order;
                    numsteps=1, time=safety.scale * step_time,
                    term_ordering=ordering,
                )) * exp(-im * phase_shift)
                spectrum = eigvals(Hermitian(U + U'))
                push!(spectra, spectrum)
                push!(scales, safety.scale)
                push!(energies, safe_cosine_value_to_energy(
                    last(spectrum), info.shift, info.normalization,
                    step_time, safety.scale, phase_shift,
                ))
            end
            if !safe
                @assert maximum(abs, spectra[1] - spectra[2]) < 1e-12 "Cosine spectra differ: $path, $order, $nsteps"
                @assert abs(energies[1] - energies[2]) < 1e-10 "Energies differ: $path, $order, $nsteps"
            end
            @printf("%s,%d,%s,%d,%s,%s,%d,%.16g,%.16g,%.8e,%.8e,%.8e,%.8e,%.16g,%.16g\n",
                relpath(path, root), n, string(order), nsteps, string(safe),
                string(exact_reverse), tied_noncommuting, energies...,
                energies[1] - reference_energy, energies[2] - reference_energy,
                abs(energies[1] - energies[2]), maximum(abs, spectra[1] - spectra[2]),
                scales...)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    withenv("COMM_METHOD" => "exact") do
        check_magnitude_orderings()
    end
end
