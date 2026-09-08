module SpectralBoundTests

using LinearAlgebra
using SparseArrays
using Test

include(joinpath(@__DIR__, "..", "src", "trotter_spectral_energy.jl"))

const I2 = ComplexF64[1 0; 0 1]
const X = ComplexF64[0 1; 1 0]
const Z = ComplexF64[1 0; 0 -1]

function dense_anticommuting_suffix_weights(pauli_strings, abs_coeffs)
    paulis = Matrix.(OP_from_string.(pauli_strings))
    weights = zeros(Float64, length(pauli_strings))
    for j in 1:length(paulis)-1, k in j+1:length(paulis)
        paulis[j] * paulis[k] == paulis[k] * paulis[j] ||
            (weights[j] += abs_coeffs[k])
    end
    return weights
end

function pauli_string(nqubits, entries::Pair{Int,Char}...)
    paulis = fill('I', nqubits)
    for (site, pauli) in entries
        paulis[site] = pauli
    end
    return String(paulis)
end

@testset "Hermitian operator norm on small matrices" begin
    for diagonal in ([3.0], [-2.0, 0.5], [-4.0, -1.0, 2.0, 3.0])
        A = sparse(ComplexF64.(diagm(0 => diagonal)))
        @test hermitian_opnorm(A) ≈ opnorm(Matrix(A)) rtol=1e-12 atol=1e-12
    end
end

@testset "Bit-packed Pauli commutativity" begin
    # Dense matrix multiplication supplies an independent commutativity test.
    # Exhausting all three-qubit strings exercises I, X, Y, and Z and every
    # possible odd/even parity of local anticommutations.
    paulis = ('I', 'X', 'Y', 'Z')
    pauli_strings = [string(a, b, c) for a in paulis for b in paulis for c in paulis]
    abs_coeffs = Float64.(eachindex(pauli_strings))
    dense_weights = dense_anticommuting_suffix_weights(
        pauli_strings, abs_coeffs
    )
    fast_weights = pauli_anticommuting_suffix_weights(
        pauli_strings, abs_coeffs
    )
    @test fast_weights == dense_weights

    # One local X/Z crossing anticommutes; two crossings commute.
    @test pauli_anticommuting_suffix_weights(
        ["XI", "ZI"], [1.0, 2.0]
    ) == [2.0, 0.0]
    @test pauli_anticommuting_suffix_weights(
        ["XX", "ZZ"], [1.0, 2.0]
    ) == [0.0, 0.0]

    # Put crossings in all three UInt64 chunks. The first pair has two
    # crossings and commutes; the second has one and anticommutes.
    p = pauli_string(130, 1 => 'X', 65 => 'X', 130 => 'X')
    two_crossings = pauli_string(130, 1 => 'Z', 65 => 'Z')
    one_crossing = pauli_string(130, 130 => 'Z')
    @test pauli_anticommuting_suffix_weights(
        [p, two_crossings, one_crossing], [1.0, 2.0, 3.0]
    ) == [3.0, 0.0, 0.0]

    @test_throws DimensionMismatch pauli_anticommuting_suffix_weights(
        ["X", "ZZ"], [1.0, 2.0]
    )
    @test_throws DimensionMismatch pauli_anticommuting_suffix_weights(
        ["X", "Z"], [1.0]
    )
    @test_throws ArgumentError pauli_anticommuting_suffix_weights(
        ["X", "A"], [1.0, 2.0]
    )
end

@testset "Analytic Strang commutator bound" begin
    term_mats = SparseMatrixCSC{ComplexF64,Int}[
        sparse(X),
        sparse(Z),
    ]
    # [Z,[Z,X]] = 4X and [X,[X,Z]] = 4Z, so Proposition 16
    # gives 4/12 + 4/24 = 1/2.
    @test commutator_bound_prefactor_exact(term_mats) ≈ 0.5 atol=1e-12

    ham = Dict("I" => 0.0 + 0.0im, "X" => 1.0 + 0.0im, "Z" => 1.0 + 0.0im)
    bounds, method = withenv("COMM_METHOD" => "exact") do
        commutator_error_bounds(
            ham,
            4.0,
            1;
            nsteps_list=[1, 2, 4],
            time=0.8,
        )
    end
    @test method == "exact"
    @test bounds ≈ [0.5 * (0.8 / (4n))^3 * n for n in (1, 2, 4)]
    @test bounds[1] / bounds[2] ≈ 4.0
    @test bounds[2] / bounds[3] ≈ 4.0

    fast_bounds, method = withenv("COMM_METHOD" => "fast") do
        commutator_error_bounds(
            ham,
            4.0,
            1;
            nsteps_list=[1],
            time=0.8,
        )
    end
    @test method == "fast"
    @test only(fast_bounds) ≈ bounds[1]

    @test commutator_bound_prefactor_fast(
        ["XI", "IZ"], [1.0, 2.0]; nqubits=2
    ) ≈ 0.0
    @test commutator_bound_prefactor_fast(
        ["XI", "ZI", "IZ"], ones(3); nqubits=2
    ) ≈ 5 / 6

    # The fast expression must remain a certified upper bound on the dense
    # nested-commutator calculation after commuting pairs are discarded.
    pauli_strings = ["XI", "ZI", "IZ", "XX", "YY"]
    abs_coeffs = [0.7, 1.1, 0.4, 0.9, 1.3]
    term_mats = SparseMatrixCSC{ComplexF64,Int}[
        coefficient * sparse(OP_from_string(pauli))
        for (pauli, coefficient) in zip(pauli_strings, abs_coeffs)
    ]
    exact_prefactor = commutator_bound_prefactor_exact(term_mats)
    for norm_mode in (:l1, :fro)
        fast_prefactor = commutator_bound_prefactor_fast(
            pauli_strings,
            abs_coeffs;
            nqubits=2,
            norm_mode=norm_mode,
        )
        @test fast_prefactor + 1e-12 >= exact_prefactor
    end

    commuting_strings = ["ZI", "IZ", "ZZ"]
    commuting_coeffs = [0.7, 1.1, 2.0]
    @test commutator_bound_prefactor_fast(
        commuting_strings, commuting_coeffs; nqubits=2
    ) == 0.0
end

@testset "Analytic first-order commutator bound" begin
    term_mats = SparseMatrixCSC{ComplexF64,Int}[
        sparse(X),
        sparse(Z),
    ]
    # ‖[Z,X]‖ = 2, so the first-order prefactor is (1/2)*2 = 1.
    @test first_order_commutator_bound_prefactor_exact(term_mats) ≈ 1.0 atol=1e-12
    @test first_order_commutator_bound_prefactor_fast(
        ["X", "Z"], [1.0, 1.0]
    ) ≈ 1.0
    @test first_order_commutator_bound_prefactor_fast(
        ["XI", "IX", "XX"], ones(3)
    ) ≈ 0.0
    @test pauli_anticommuting_suffix_weights(
        ["XI", "ZI", "IZ"], [1.0, 2.0, 4.0]
    ) == [2.0, 0.0, 0.0]
    @test pauli_anticommuting_suffix_weights(
        ["I"^64 * "X", "I"^64 * "Z"], [1.0, 3.0]
    ) == [3.0, 0.0]

    # For first order, the fast result is exactly the coefficient-weighted
    # sum over anticommuting pairs. It must also bound the dense expression,
    # whose suffix commutators can exhibit cancellation.
    pauli_strings = ["XI", "ZI", "IZ", "XX", "YY"]
    abs_coeffs = [0.7, 1.1, 0.4, 0.9, 1.3]
    dense_weights = dense_anticommuting_suffix_weights(
        pauli_strings, abs_coeffs
    )
    fast_prefactor = first_order_commutator_bound_prefactor_fast(
        pauli_strings, abs_coeffs
    )
    @test fast_prefactor == sum(abs_coeffs .* dense_weights)
    term_mats = SparseMatrixCSC{ComplexF64,Int}[
        coefficient * sparse(OP_from_string(pauli))
        for (pauli, coefficient) in zip(pauli_strings, abs_coeffs)
    ]
    @test fast_prefactor + 1e-12 >=
          first_order_commutator_bound_prefactor_exact(term_mats)

    ham = Dict("I" => 0.0 + 0.0im, "X" => 1.0 + 0.0im, "Z" => 1.0 + 0.0im)
    normalization = 4.0
    total_time = 0.8
    nsteps_list = [1, 2, 4]
    bounds = first_order_commutator_error_bounds(
        ham,
        normalization,
        1;
        nsteps_list=nsteps_list,
        time=total_time,
    )
    @test bounds ≈ [total_time^2 / (normalization^2 * n) for n in nsteps_list]
    @test bounds[1] / bounds[2] ≈ 2.0
    @test bounds[2] / bounds[3] ≈ 2.0

    exact = exp(-im * total_time / 2) *
            exp(-im * total_time * (X + Z) / normalization)
    for (nsteps, bound) in zip(nsteps_list, bounds)
        approximation = reference_first_order_trotter_unitary(
            ham,
            normalization,
            1;
            numsteps=nsteps,
            time=total_time,
        )
        @test opnorm(Matrix(approximation) - exact) <= bound + 1e-14
    end

    commuting_terms = SparseMatrixCSC{ComplexF64,Int}[sparse(Z), sparse(2Z)]
    @test first_order_commutator_bound_prefactor_exact(commuting_terms) ≈ 0.0 atol=1e-14
    @test_throws ArgumentError first_order_commutator_error_bounds(
        ham, normalization, 1; nsteps_list=[0], time=total_time
    )

    first_safety = trotter_safety_parameters(ham, normalization, 1, 2, :first)
    first_error_bound = π^2 / (2 * normalization^2)
    @test first_safety.error_bound ≈ first_error_bound
    @test first_safety.shift ≈ 2asin(first_error_bound / 2)
    @test first_safety.scale ≈ (π - 2first_safety.shift) / π
    @test first_safety.scale * 0 + first_safety.shift ≈ first_safety.shift
    @test first_safety.scale * π + first_safety.shift ≈ π - first_safety.shift

    exact_energy = -sqrt(2)
    step_time = π / 2
    safe_phase_shift = first_safety.shift / 2
    safe_phase = first_safety.scale * step_time *
                 (0.5 + exact_energy / normalization) + safe_phase_shift
    @test safe_cosine_value_to_energy(
        2cos(safe_phase),
        0.0,
        normalization,
        step_time,
        first_safety.scale,
        safe_phase_shift,
    ) ≈ exact_energy

    second_safety = trotter_safety_parameters(ham, normalization, 1, 2, :second)
    second_error_bound = π^3 / (8 * normalization^3)
    @test second_safety.error_bound ≈ second_error_bound
    @test second_safety.shift ≈ 2asin(second_error_bound / 2)
    @test second_safety.scale ≈ (π - 2second_safety.shift) / π
end

@testset "Trotter term ordering" begin
    ham = Dict(
        "I" => 0.0 + 0.0im,
        "X" => 0.2 + 0.0im,
        "Y" => -0.8 + 0.0im,
        "Z" => 0.5 + 0.0im,
    )
    @test first.(ordered_hamiltonian_terms(ham, 1)) == ["Y", "Z", "X"]
    @test [coefficient for (coefficient, _) in build_trotter_terms(ham, 1)] ==
          [-0.8, 0.5, 0.2]
    @test first.(ordered_hamiltonian_terms(
        ham, 1; term_ordering=:increasing_magnitude
    )) == ["X", "Z", "Y"]
    @test [
        coefficient for (coefficient, _) in build_trotter_terms(
            ham, 1; term_ordering=:increasing_magnitude
        )
    ] == [0.2, 0.5, -0.8]

    explicit_order = ["X", "Y", "Z"]
    normalization = 2.0
    time = 0.3
    rotations = Dict(
        pauli => exp(-im * real(ham[pauli]) * Matrix(OP_from_string(pauli)) *
                     time / normalization)
        for pauli in explicit_order
    )
    expected_first = exp(-im * time / 2) *
                     rotations["Z"] * rotations["Y"] * rotations["X"]
    actual_first = reference_trotter_unitary_by_order(
        ham,
        normalization,
        1,
        :first;
        numsteps=1,
        time=time,
        term_ordering=explicit_order,
    )
    @test Matrix(actual_first) ≈ expected_first atol=2e-14

    half_rotations = Dict(
        pauli => exp(-im * real(ham[pauli]) * Matrix(OP_from_string(pauli)) *
                     time / (2normalization))
        for pauli in explicit_order
    )
    expected_second = exp(-im * time / 2) *
                      half_rotations["X"] * half_rotations["Y"] *
                      rotations["Z"] * half_rotations["Y"] *
                      half_rotations["X"]
    actual_second = reference_trotter_unitary_by_order(
        ham,
        normalization,
        1,
        :second;
        numsteps=1,
        time=time,
        term_ordering=explicit_order,
    )
    @test Matrix(actual_second) ≈ expected_second atol=2e-14

    # The two magnitude strategies must reach the product-formula loops in
    # opposite coefficient order and produce different unitaries when the
    # Pauli terms do not commute.
    for order in (:first, :second)
        decreasing_unitary = reference_trotter_unitary_by_order(
            ham,
            normalization,
            1,
            order;
            numsteps=1,
            time=time,
            term_ordering=:magnitude,
        )
        increasing_unitary = reference_trotter_unitary_by_order(
            ham,
            normalization,
            1,
            order;
            numsteps=1,
            time=time,
            term_ordering=:increasing_magnitude,
        )
        @test opnorm(Matrix(decreasing_unitary - increasing_unitary)) > 1e-8
        if order == :second
            decreasing_spectrum = eigvals(Hermitian(Matrix(
                decreasing_unitary + decreasing_unitary',
            )))
            increasing_spectrum = eigvals(Hermitian(Matrix(
                increasing_unitary + increasing_unitary',
            )))
            @test decreasing_spectrum ≈ increasing_spectrum atol=2e-14
        end
    end

    # With transpose-symmetric Pauli factors, reversing a first-order product
    # transposes the unitary. U + U' then has the same eigenvalues even though
    # the two unitary matrices are different.
    real_ham = Dict(
        "II" => 0.0 + 0.0im,
        "XI" => 0.9 + 0.0im,
        "ZI" => 0.6 + 0.0im,
        "XZ" => 0.2 + 0.0im,
    )
    decreasing_unitary = reference_trotter_unitary_by_order(
        real_ham,
        normalization,
        2,
        :first;
        numsteps=1,
        time=time,
        term_ordering=:magnitude,
    )
    increasing_unitary = reference_trotter_unitary_by_order(
        real_ham,
        normalization,
        2,
        :first;
        numsteps=1,
        time=time,
        term_ordering=:increasing_magnitude,
    )
    @test opnorm(Matrix(decreasing_unitary - increasing_unitary)) > 1e-8
    @test eigvals(Hermitian(Matrix(decreasing_unitary + decreasing_unitary'))) ≈
          eigvals(Hermitian(Matrix(increasing_unitary + increasing_unitary'))) atol=2e-14

    first_bound = only(first_order_commutator_error_bounds(
        ham,
        normalization,
        1;
        nsteps_list=[2],
        term_ordering=explicit_order,
    ))
    safety = trotter_safety_parameters(
        ham,
        normalization,
        1,
        2,
        :first;
        term_ordering=explicit_order,
    )
    @test safety.error_bound ≈ first_bound

    second_bound, _ = withenv("COMM_METHOD" => "exact") do
        commutator_error_bounds(
            ham,
            normalization,
            1;
            nsteps_list=[2],
            term_ordering=explicit_order,
        )
    end
    second_safety = withenv("COMM_METHOD" => "exact") do
        trotter_safety_parameters(
            ham,
            normalization,
            1,
            2,
            :second;
            term_ordering=explicit_order,
        )
    end
    @test second_safety.error_bound ≈ only(second_bound)

    metadata = Dict(
        "number of qubits" => "2",
        "number of active, occupied, single-occupancy orbitals" => "0",
        "one-norm of sum of Pauli strings" => "1.5",
    )
    two_qubit_ham = Dict(
        "II" => 0.0 + 0.0im,
        "XI" => 0.2 + 0.0im,
        "YI" => -0.8 + 0.0im,
        "ZI" => 0.5 + 0.0im,
    )
    details = trotter_energy(
        metadata,
        two_qubit_ham,
        2;
        method=:arpack,
        order=:first,
        nev=2,
        return_details=true,
    )
    expected_safety = trotter_safety_parameters(
        two_qubit_ham,
        normalize_hamiltonian(metadata, two_qubit_ham).normalization,
        2,
        2,
        :first,
    )
    @test details.safe_scaling_factor ≈ expected_safety.scale
    unsafe_details = trotter_energy(
        metadata,
        two_qubit_ham,
        2;
        method=:arpack,
        order=:first,
        nev=2,
        safe_normalization=false,
        return_details=true,
    )
    @test unsafe_details.safe_scaling_factor == 1.0
end

@testset "Trotter unitary actions" begin
    ham = Dict("II" => 0.0 + 0.0im, "XI" => 0.7 + 0.0im, "ZI" => -0.4 + 0.0im)
    terms = build_trotter_terms(ham, 2)
    normalization = 2.5
    total_time = 0.6
    nsteps = 3
    psi = ComplexF64[1, 2im, -1, 0.5]
    psi /= norm(psi)

    for order in (:first, :second)
        U = reference_trotter_unitary_by_order(
            ham,
            normalization,
            2,
            order;
            numsteps=nsteps,
            time=total_time,
        )
        identity4 = Matrix{ComplexF64}(I, 4, 4)
        rotation((coefficient, P), step_time) =
            exp(-im * coefficient * Matrix(P) * step_time / normalization)
        sequence = order == :first ?
                   [(term, total_time / nsteps) for term in terms] :
                   vcat(
                       [(term, total_time / (2nsteps)) for term in terms],
                       [(term, total_time / (2nsteps)) for term in reverse(terms)],
                   )
        expected_step = foldl(
            (product, item) -> rotation(item[1], item[2]) * product,
            sequence;
            init=identity4,
        )
        expected = exp(-im * total_time / 2) * expected_step^nsteps

        @test Matrix(U) ≈ expected atol=2e-14
        @test Matrix(U)' * Matrix(U) ≈ Matrix{ComplexF64}(I, 4, 4) atol=2e-14
        @test apply_trotter_unitary(
            psi, terms, nsteps, normalization, order; time=total_time
        ) ≈ U * psi atol=2e-14
        @test apply_trotter_unitary_adjoint(
            psi, terms, nsteps, normalization, order; time=total_time
        ) ≈ U' * psi atol=2e-14
    end

    # With a single Pauli term every product formula equals the exact exponential.
    single_term_ham = Dict("II" => 0.0 + 0.0im, "XI" => 0.7 + 0.0im)
    expected = exp(-im * total_time / 2) *
               exp(-im * (0.7 / normalization) * kron(X, I2) * total_time)
    for order in (:first, :second, :fourth, :sixth)
        U = reference_trotter_unitary_by_order(
            single_term_ham,
            normalization,
            2,
            order;
            numsteps=2,
            time=total_time,
        )
        @test Matrix(U) ≈ expected atol=3e-14
    end
end

@testset "Spectral ground energy" begin
    metadata = Dict(
        "number of qubits" => "2",
        "number of active, occupied, single-occupancy orbitals" => "1",
        "one-norm of sum of Pauli strings" => "1.9",
    )
    ham = Dict("II" => 1.2 + 0.0im, "ZI" => 0.7 + 0.0im)
    exact_ground_energy = 1.2 - 0.7

    normalization = normalize_hamiltonian(metadata, ham).normalization
    cosine_eigenvalue = 2cos(0.8 * (0.5 - 0.7 / normalization))
    @test cosine_value_to_energy(cosine_eigenvalue, 1.2, normalization, 0.8) ≈
          exact_ground_energy

    @test trotter_single_step_time(π, 5) ≈ π / 5
    @test_throws ArgumentError trotter_single_step_time(π, 0)
    @test trotter_arnoldi_tolerance(1) ≈ 1e-12
    @test trotter_arnoldi_tolerance(10) ≈ 1e-14
    @test trotter_arnoldi_tolerance(100) == eps(Float64)
    @test_throws ArgumentError trotter_arnoldi_tolerance(0)

    for safe_normalization in (true, false),
        method in (:arpack, :arnoldi),
        order in (:first, :second)
        energy = trotter_energy(
            metadata,
            ham,
            3;
            method=method,
            order=order,
            time=0.8,
            nev=2,
            safe_normalization=safe_normalization,
        )
        @test energy ≈ exact_ground_energy atol=2e-12
    end
end

@testset "Short-time spectral signal" begin
    metadata = Dict(
        "number of qubits" => "2",
        "number of active, occupied, single-occupancy orbitals" => "1",
        "one-norm of sum of Pauli strings" => "2.0",
    )
    ham = Dict("II" => 0.0 + 0.0im, "XI" => 1.0 + 0.0im, "ZI" => 1.0 + 0.0im)
    normalization = normalize_hamiltonian(metadata, ham).normalization
    nsteps = 25
    step_time = π / nsteps

    for order in (:first, :second)
        safety = trotter_safety_parameters(ham, normalization, 2, nsteps, order)
        safe_step_time = safety.scale * step_time
        safety_phase_shift = safety.shift / nsteps
        one_step_unitary = reference_trotter_unitary_by_order(
            ham,
            normalization,
            2,
            order;
            numsteps=1,
            time=safe_step_time,
        )
        one_step_unitary .*= exp(-im * safety_phase_shift)
        cosine_eigenvalues = eigvals(Hermitian(Matrix(
            one_step_unitary + one_step_unitary',
        )))
        expected = minimum(safe_cosine_value_to_energy.(
            cosine_eigenvalues,
            0.0,
            normalization,
            step_time,
            safety.scale,
            safety_phase_shift,
        ))

        for method in (:arpack, :arnoldi)
            energy = trotter_energy(
                metadata,
                ham,
                nsteps;
                method=method,
                order=order,
                time=π,
                nev=2,
            )
            @test energy ≈ expected atol=2e-10
        end
    end
end

end
