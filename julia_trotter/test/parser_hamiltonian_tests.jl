module ParserHamiltonianTests

using LinearAlgebra
using SparseArrays
using Test

include(joinpath(@__DIR__, "..", "src", "parser.jl"))
include(joinpath(@__DIR__, "..", "src", "hamiltonian_utils.jl"))

const I2 = ComplexF64[1 0; 0 1]
const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const Z = ComplexF64[1 0; 0 -1]

@testset "Hamiltonian parser" begin
    mktemp() do path, io
        write(io, """
        # number of qubits = 2
        # one-norm of sum of Pauli strings = 4.5
        # free-form comment
        1.5 II
        0.5 XI
        0.25 XI
        -0.25 IZ
        0.1+0.2j YY
        """)
        close(io)

        metadata, ham = parse_hamiltonian_file(path)
        @test metadata["number of qubits"] == "2"
        @test metadata["free-form comment"] == ""
        @test ham["XI"] == 0.75 + 0.0im
        @test ham["YY"] == 0.1 + 0.2im
        @test length(ham) == 4
    end

    @test parse_coeff("-2.5e-3") == -2.5e-3 + 0.0im
    @test parse_coeff("1.25-0.5i") == 1.25 - 0.5im
    @test_throws ErrorException parse_coeff("not-a-number")
end

@testset "Pauli and Hamiltonian construction" begin
    @test Matrix(OP_from_string("I")) == I2
    @test Matrix(OP_from_string("X")) == X
    @test Matrix(OP_from_string("Y")) == Y
    @test Matrix(OP_from_string("Z")) == Z
    @test Matrix(OP_from_string("XY")) == kron(X, Y)

    metadata = Dict(
        "number of qubits" => "2",
        "one-norm of sum of Pauli strings" => "4.5",
    )
    ham = Dict(
        "II" => 1.5 + 0.0im,
        "XI" => 0.75 + 0.0im,
        "IZ" => -0.25 + 0.0im,
    )

    norm_info = normalize_hamiltonian(metadata, ham)
    @test norm_info.shift == 1.5
    @test norm_info.norm_bound == 3.0
    @test norm_info.normalization == 6.0

    @test first.(ordered_hamiltonian_terms(ham, 2)) == ["XI", "IZ"]
    @test first.(ordered_hamiltonian_terms(
        ham, 2; term_ordering=:lexicographic
    )) == ["IZ", "XI"]
    @test first.(ordered_hamiltonian_terms(
        ham, 2; term_ordering=:increasing_magnitude
    )) == ["IZ", "XI"]
    dict_order = [pauli for pauli in keys(ham) if pauli != "II"]
    @test first.(ordered_hamiltonian_terms(
        ham, 2; term_ordering=:dict
    )) == dict_order
    @test first.(ordered_hamiltonian_terms(
        ham, 2; term_ordering=["IZ", "XI"]
    )) == ["IZ", "XI"]

    tied_ham = Dict(
        "II" => 0.0 + 0.0im,
        "ZI" => -1.0 + 0.0im,
        "IX" => 1.0 + 0.0im,
        "XX" => 0.5 + 0.0im,
    )
    @test first.(ordered_hamiltonian_terms(tied_ham, 2)) == ["IX", "ZI", "XX"]
    @test first.(ordered_hamiltonian_terms(
        tied_ham, 2; term_ordering=:increasing_magnitude
    )) == ["XX", "IX", "ZI"]
    @test_throws ArgumentError ordered_hamiltonian_terms(
        ham, 2; term_ordering=:unknown
    )
    @test_throws ArgumentError ordered_hamiltonian_terms(
        ham, 2; term_ordering=["XI"]
    )
    @test_throws ArgumentError ordered_hamiltonian_terms(
        ham, 2; term_ordering=["XI", "XI"]
    )

    expected = 1.5 * kron(I2, I2) + 0.75 * kron(X, I2) - 0.25 * kron(I2, Z)
    @test Matrix(build_sparse_hamiltonian(ham, 2)) ≈ expected

    raw_terms = build_hamiltonian_terms(ham, metadata; normalize=false, scale_by_pi=false)
    @test raw_terms.shift == 1.5
    @test raw_terms.normalization == 1.0
    @test [coefficient for (coefficient, _) in raw_terms.H_terms] == [0.75, -0.25]
    @test sum(coefficient * Matrix(P) for (coefficient, P) in raw_terms.H_terms) ≈
          expected - 1.5 * kron(I2, I2)

    custom_terms = build_hamiltonian_terms(
        ham,
        metadata;
        normalize=false,
        scale_by_pi=false,
        term_ordering=["IZ", "XI"],
    )
    @test [coefficient for (coefficient, _) in custom_terms.H_terms] == [-0.25, 0.75]

    scaled_terms = build_hamiltonian_terms(ham, metadata; normalize=true, scale_by_pi=true)
    @test sum(coefficient * Matrix(P) for (coefficient, P) in scaled_terms.H_terms) ≈
          (π / 6) * (expected - 1.5 * kron(I2, I2))
end

@testset "Hartree–Fock basis state" begin
    state = construct_hf_state(4, 2)
    expected = zeros(ComplexF64, 16)
    expected[13] = 1
    @test state == expected
    @test construct_hf_state(3, 0)[1] == 1
    @test construct_hf_state(3, 3)[8] == 1
    @test norm(state) == 1
    @test_throws ArgumentError construct_hf_state(-1, 0)
    @test_throws ArgumentError construct_hf_state(3, 4)
end

end
