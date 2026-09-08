module GroundEnergyScriptTests

using Test

include(joinpath(@__DIR__, "..", "trotter_ground_energy.jl"))

@testset "Ground-energy command-line options" begin
    defaults = parse_cli(String[])
    @test defaults.trotter_order == :second
    @test defaults.term_ordering == :magnitude
    @test defaults.safe_normalization
    @test !defaults.use_arnoldi
    @test !defaults.benchmark

    options = parse_cli([
        "--arnoldi",
        "--trotter-order=first",
        "--term-order=increasing-magnitude",
        "--safe-normalization=false",
        "--benchmark",
        "example.dat",
    ])
    @test options.use_arnoldi
    @test options.trotter_order == :first
    @test options.term_ordering == :increasing_magnitude
    @test !options.safe_normalization
    @test options.benchmark
    @test options.filepath == "example.dat"

    @test parse_cli(["--term-order=XX,ZI,IZ"]).term_ordering ==
          ["XX", "ZI", "IZ"]
    @test !parse_cli(["--no-safe-normalization"]).safe_normalization
    @test parse_cli(["--safe-normalization=TRUE"]).safe_normalization
    @test parse_cli(["--trotter-order=1"]).trotter_order == :first
    @test parse_cli(["--trotter-order=2"]).trotter_order == :second
    @test !parse_cli(["--benchmark=false"]).benchmark
    @test_throws ErrorException parse_cli(["--safe-normalization=maybe"])
    @test_throws ErrorException parse_cli(["--benchmark=maybe"])
    @test_throws ErrorException parse_cli(["--trotter-order=fourth"])
    @test_throws ErrorException parse_cli(["--term-order=unknown"])
end

@testset "Optional diagonalization benchmark" begin
    mktemp() do path, io
        write(io, """
        # number of qubits = 3
        # number of active, occupied, single-occupancy orbitals = 0
        # one-norm of sum of Pauli strings = 0.5
        # smallest eigenvalue = -0.5
        0.0 III
        0.5 ZII
        """)
        close(io)

        benchmarked = redirect_stdout(devnull) do
            main(
                path;
                trotter_order=:first,
                safe_normalization=false,
                benchmark=true,
            )
        end
        @test length(benchmarked) == length(DEFAULT_NSTEPS)
        @test all(result -> result.time_seconds isa Float64, benchmarked)
        @test all(result -> result.time_seconds >= 0, benchmarked)
        @test all(result -> result.trotter_order == :first, benchmarked)
        @test all(result -> result.safe_scaling_factor == 1.0, benchmarked)
        @test all(result -> isapprox(result.energy, -0.5; atol=1e-12), benchmarked)

        mktemp() do _, ordinary_output
            ordinary = redirect_stdout(ordinary_output) do
                main(path; safe_normalization=false, benchmark=false)
            end
            flush(ordinary_output)
            seekstart(ordinary_output)
            @test all(result -> isnothing(result.time_seconds), ordinary)
            @test !occursin("safe_scaling_factor", read(ordinary_output, String))
        end

        mktemp() do _, safe_output
            safe_results = redirect_stdout(safe_output) do
                main(path; safe_normalization=true, benchmark=false)
            end
            flush(safe_output)
            seekstart(safe_output)
            @test occursin("safe_scaling_factor", read(safe_output, String))
            @test all(result -> result.safe_scaling_factor == 1.0, safe_results)
        end
        @test_throws ArgumentError redirect_stdout(devnull) do
            main(path; trotter_order=:fourth)
        end
    end
end

end
