using Test

@testset "Julia Trotter" begin
    include("parser_hamiltonian_tests.jl")
    include("evolution_tests.jl")
    include("spectral_bound_tests.jl")
    include("ground_energy_script_tests.jl")
end
