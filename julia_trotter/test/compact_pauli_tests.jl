module CompactPauliTests

using LinearAlgebra
using SparseArrays
using Test

include(joinpath(@__DIR__, "..", "src", "quantum_utils.jl"))
include(joinpath(@__DIR__, "..", "src", "statevector_simulators.jl"))

# OP_from_string may return a real-eltype sparse matrix (e.g. for "ZZI"); the
# legacy tuple methods dispatch on SparseMatrixCSC{ComplexF64,Int}, so convert.
sparse_op(s) = SparseMatrixCSC{ComplexF64,Int}(OP_from_string(s))

# Build a normalized random-ish (deterministic) complex state vector of length n.
function make_state(n)
    psi = ComplexF64[cos(0.3k) + im * sin(0.17k + 1) for k in 1:n]
    return psi / norm(psi)
end

@testset "Bit convention (leftmost_is_msb)" begin
    # OP_from_string("XI") = kron(X, I): the leftmost char acts on the MSB.
    # It should send |00⟩ (column 1) to |10⟩ (row 3).
    XI = OP_from_string("XI")
    e0 = ComplexF64[1, 0, 0, 0]
    @test XI * e0 ≈ ComplexF64[0, 0, 1, 0]

    term = compact_pauli_term("XI", 1.0)
    @test term.xmask == UInt64(0b10)   # X on the MSB -> bit 1
    @test term.zmask == UInt64(0)
end

@testset "apply_pauli_rotation! matches OP_from_string" begin
    strings = ["I", "X", "Y", "Z", "XI", "IX", "XY", "YX", "YZ", "ZY", "XYZ"]
    coeff = 0.63
    dt = 0.41
    theta = coeff * dt
    for p in strings
        nqubits = length(p)
        n = 2^nqubits
        psi = make_state(n)

        P = Matrix(OP_from_string(p))
        expected = cos(theta) * psi - im * sin(theta) * (P * psi)

        got = copy(psi)
        apply_pauli_rotation!(got, compact_pauli_term(p, coeff), dt)
        @test got ≈ expected atol = 1e-13

        # Also compare against the true matrix exponential.
        @test got ≈ exp(-im * coeff * P * dt) * psi atol = 1e-12
    end
end

@testset "hamiltonian_matvec matches sparse backend" begin
    strings = ["XYZ", "ZZI", "IXX", "YIZ"]
    coeffs = [0.7, -0.35, 0.22, 0.9]
    n = 2^3
    psi = make_state(n)

    sparse_terms = [(coeffs[i], sparse_op(strings[i])) for i in eachindex(strings)]
    compact_terms = [compact_pauli_term(strings[i], coeffs[i]) for i in eachindex(strings)]

    @test hamiltonian_matvec(psi, compact_terms) ≈ hamiltonian_matvec(psi, sparse_terms) atol = 1e-13
end

@testset "Compact Trotter matches sparse Trotter" begin
    strings = ["XYZ", "ZZI", "IXX", "YIZ"]
    coeffs = [0.7, -0.35, 0.22, 0.9]
    n = 2^3
    psi = make_state(n)
    dt = 0.15
    nsteps = 3

    sparse_terms = [(coeffs[i], sparse_op(strings[i])) for i in eachindex(strings)]
    compact_terms = [compact_pauli_term(strings[i], coeffs[i]) for i in eachindex(strings)]

    for method in (first_order_trotter_statevec,
                   second_order_trotter_statevec,
                   fourth_order_trotter_statevec)
        sparse_result = method(sparse_terms, psi, dt, nsteps)
        compact_result = method(compact_terms, psi, dt, nsteps)
        @test compact_result ≈ sparse_result atol = 1e-12
    end
end

@testset "Threaded path (>= 2^15 amplitudes)" begin
    # 15-qubit string exercises the threaded branch (n = 2^15 = PAULI_THREAD_THRESHOLD).
    p = "XYZIXZYIXZYIXZY"  # 15 chars
    @test length(p) == 15
    coeff = 0.53
    dt = 0.37
    n = 2^15
    psi = make_state(n)

    P = OP_from_string(p)
    expected = cos(coeff * dt) * psi - im * sin(coeff * dt) * (P * psi)

    got = copy(psi)
    apply_pauli_rotation!(got, compact_pauli_term(p, coeff), dt)
    @test got ≈ expected atol = 1e-11

    # Diagonal threaded branch.
    pz = "ZZZZZZZZZZZZZZZ"
    Pz = OP_from_string(pz)
    expected_z = cos(coeff * dt) * psi - im * sin(coeff * dt) * (Pz * psi)
    got_z = copy(psi)
    apply_pauli_rotation!(got_z, compact_pauli_term(pz, coeff), dt)
    @test got_z ≈ expected_z atol = 1e-11
end

end
