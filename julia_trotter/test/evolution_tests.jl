module EvolutionTests

using LinearAlgebra
using SparseArrays
using Test

include(joinpath(@__DIR__, "..", "src", "statevector_simulators.jl"))

const X = ComplexF64[0 1; 1 0]
const Z = ComplexF64[1 0; 0 -1]

@testset "Pauli rotation and Hamiltonian action" begin
    psi = ComplexF64[1, 2im]
    psi /= norm(psi)
    rotated = copy(psi)
    apply_pauli_rotation!(rotated, sparse(X), 0.7, 0.3)
    @test rotated ≈ exp(-im * 0.7 * X * 0.3) * psi atol=1e-14
    @test norm(rotated) ≈ 1 atol=1e-14

    terms = [(0.7, sparse(X)), (-0.2, sparse(Z))]
    @test hamiltonian_matvec(psi, terms) ≈ (0.7 * X - 0.2 * Z) * psi
end

@testset "Product-formula convergence orders" begin
    terms = [(0.7, sparse(X)), (0.4, sparse(Z))]
    psi = ComplexF64[1, im]
    psi /= norm(psi)
    total_time = 0.8
    exact = exp(-im * (0.7 * X + 0.4 * Z) * total_time) * psi

    methods = (
        (first_order_trotter_statevec, 2.0),
        (second_order_trotter_statevec, 4.0),
        (fourth_order_trotter_statevec, 16.0),
    )
    for (method, expected_ratio) in methods
        errors = [norm(method(terms, psi, total_time / nsteps, nsteps) - exact)
                  for nsteps in (4, 8, 16)]
        @test errors[3] < errors[2] < errors[1]
        @test errors[1] / errors[2] ≈ expected_ratio rtol=0.02
        @test errors[2] / errors[3] ≈ expected_ratio rtol=0.02
        @test norm(method(terms, psi, total_time / 8, 8)) ≈ 1 atol=2e-14
    end

    one_term = [(0.7, sparse(X))]
    exact_one_term = exp(-im * 0.7 * X * total_time) * psi
    @test first_order_trotter_statevec(one_term, psi, total_time / 3, 3) ≈ exact_one_term
    @test second_order_trotter_statevec(one_term, psi, total_time / 3, 3) ≈ exact_one_term
    @test fourth_order_trotter_statevec(one_term, psi, total_time / 3, 3) ≈ exact_one_term
end

module ChebyshevTests

using LinearAlgebra
using SparseArrays
using Test

include(joinpath(@__DIR__, "..", "src", "densesimulators.jl"))

@testset "Chebyshev evolution" begin
    X = ComplexF64[0 1; 1 0]
    H = 0.4 * X
    psi = ComplexF64[1, 0]
    total_time = 0.7
    exact_state = exp(-im * H * total_time) * psi

    state, observables, error_info = chebyshev_simulation(
        psi,
        total_time,
        x -> H * x,
        1.0,
        [];
        α=4,
        order=12,
        returnstate=true,
    )
    @test state ≈ exact_state atol=1e-13
    @test isempty(observables)
    @test error_info.error < 1e-15

    unitary = chebyshev_exponentiation(sparse(H), 12, total_time, 1.0; α=4)
    @test Matrix(unitary) ≈ exp(-im * H * total_time) atol=1e-13
end

end

end
