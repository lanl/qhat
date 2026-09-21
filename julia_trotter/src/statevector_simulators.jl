# statevector_simulators.jl
#
# State-vector evolution methods using Pauli rotation formula.
# These methods operate on dense state vectors (not operator matrices).
#
# For Pauli strings P with P² = I, we have:
#   exp(-i*c*P*dt)|ψ⟩ = cos(c*dt)|ψ⟩ - i*sin(c*dt)*P|ψ⟩

using LinearAlgebra, SparseArrays

"""
    apply_pauli_rotation!(psi, P_sparse, coeff, dt)

Apply exp(-i*coeff*P*dt) to state vector using Pauli involution:
    exp(-i*c*P*dt)|ψ⟩ = cos(c*dt)|ψ⟩ - i*sin(c*dt)*P|ψ⟩
where P² = I (Pauli string property).

In-place operation on psi.

# Arguments
- `psi`: State vector (modified in-place)
- `P_sparse`: Sparse Pauli operator matrix
- `coeff`: Coefficient multiplying P
- `dt`: Time step

# Returns
Modified `psi` (same reference)
"""
function apply_pauli_rotation!(
    psi::Vector{ComplexF64},
    P_sparse::SparseMatrixCSC{ComplexF64, Int},
    coeff::Real,
    dt::Real
)
    theta = coeff * dt
    c, s = cos(theta), sin(theta)
    Ppsi = P_sparse * psi
    @. psi = c * psi - im * s * Ppsi
    return psi
end

"""
    first_order_trotter_statevec(H_terms, psi0, dt, nsteps)

First-order Trotter evolution:
    U(T) ≈ [exp(-iHₘdt) ⋯ exp(-iH₂dt) exp(-iH₁dt)]^nsteps

The vector order is chronological: `H₁` is applied to the state first.

# Arguments
- `H_terms`: Vector of `(coeff, P_sparse)` tuples where H = Σⱼ coeff_j * P_j
- `psi0`: Initial state vector
- `dt`: Time step
- `nsteps`: Number of Trotter steps

# Returns
Final state vector after evolution

# Error
O(dt²) per step, O(T*dt) = O(T²/nsteps) total
"""
function first_order_trotter_statevec(
    H_terms::Vector{Tuple{T, SparseMatrixCSC{ComplexF64, Int}}},
    psi0::Vector{ComplexF64},
    dt::Real,
    nsteps::Int
) where T<:Real
    psi = copy(psi0)
    for _ in 1:nsteps
        for (coeff, P) in H_terms
            apply_pauli_rotation!(psi, P, coeff, dt)
        end
    end
    return psi
end

"""
    second_order_trotter_statevec(H_terms, psi0, dt, nsteps)

Second-order Trotter evolution (Strang splitting):
    S2(dt) = exp(-iH₁dt/2) ⋯ exp(-iHₘ₋₁dt/2) exp(-iHₘdt) exp(-iHₘ₋₁dt/2) ⋯ exp(-iH₁dt/2)

Applied for `nsteps` to evolve total time T = dt * nsteps.

# Arguments
- `H_terms`: Vector of `(coeff, P_sparse)` tuples
- `psi0`: Initial state vector
- `dt`: Time step
- `nsteps`: Number of Trotter steps

# Returns
Final state vector after evolution

# Error
O(dt³) per step, O(T*dt²) = O(T³/nsteps²) total
"""
function second_order_trotter_statevec(
    H_terms::Vector{Tuple{T, SparseMatrixCSC{ComplexF64, Int}}},
    psi0::Vector{ComplexF64},
    dt::Real,
    nsteps::Int
) where T<:Real
    psi = copy(psi0)
    nterms = length(H_terms)

    for _ in 1:nsteps
        # Forward half-steps (dt/2) for all but last term
        for idx in 1:nterms-1
            coeff, P = H_terms[idx]
            apply_pauli_rotation!(psi, P, coeff, dt / 2)
        end

        # Full step (dt) for last term
        coeff_end, P_end = H_terms[end]
        apply_pauli_rotation!(psi, P_end, coeff_end, dt)

        # Backward half-steps (dt/2) in reverse order
        for idx in nterms-1:-1:1
            coeff, P = H_terms[idx]
            apply_pauli_rotation!(psi, P, coeff, dt / 2)
        end
    end
    return psi
end

"""
    fourth_order_trotter_statevec(H_terms, psi0, dt, nsteps)

Fourth-order Trotter evolution using recursive formula:
    S4(dt) = S2(p*dt)² S2((1-4p)*dt) S2(p*dt)²
where p = 1/(4 - 4^(1/3)).

# Arguments
- `H_terms`: Vector of `(coeff, P_sparse)` tuples
- `psi0`: Initial state vector
- `dt`: Time step
- `nsteps`: Number of Trotter steps

# Returns
Final state vector after evolution

# Error
O(dt⁵) per step, O(T*dt⁴) = O(T⁵/nsteps⁴) total

# References
Suzuki, "General theory of fractal path integrals", J. Math. Phys. 32, 400 (1991)
"""
function fourth_order_trotter_statevec(
    H_terms::Vector{Tuple{T, SparseMatrixCSC{ComplexF64, Int}}},
    psi0::Vector{ComplexF64},
    dt::Real,
    nsteps::Int
) where T<:Real
    p = 1 / (4 - 4^(1/3))
    psi = copy(psi0)
    nterms = length(H_terms)

    # Helper: apply one S2 step with given timestep
    function apply_s2_step!(ψ, δt)
        # Forward
        for idx in 1:nterms-1
            coeff, P = H_terms[idx]
            apply_pauli_rotation!(ψ, P, coeff, δt / 2)
        end
        # Middle
        coeff_end, P_end = H_terms[end]
        apply_pauli_rotation!(ψ, P_end, coeff_end, δt)
        # Backward
        for idx in nterms-1:-1:1
            coeff, P = H_terms[idx]
            apply_pauli_rotation!(ψ, P, coeff, δt / 2)
        end
        return ψ
    end

    for _ in 1:nsteps
        # S2(p*dt)²
        apply_s2_step!(psi, p * dt)
        apply_s2_step!(psi, p * dt)

        # S2((1-4p)*dt)
        apply_s2_step!(psi, (1 - 4p) * dt)

        # S2(p*dt)²
        apply_s2_step!(psi, p * dt)
        apply_s2_step!(psi, p * dt)
    end

    return psi
end

"""
    hamiltonian_matvec(psi, H_terms)

Compute H|ψ⟩ = Σⱼ cⱼPⱼ|ψ⟩ without forming the full Hamiltonian matrix.

Useful for Chebyshev expansion or other iterative methods.

# Arguments
- `psi`: State vector
- `H_terms`: Vector of `(coeff, P_sparse)` tuples

# Returns
Result vector H|ψ⟩
"""
function hamiltonian_matvec(
    psi::Vector{ComplexF64},
    H_terms::Vector{Tuple{T, SparseMatrixCSC{ComplexF64, Int}}}
) where T<:Real
    result = zeros(ComplexF64, size(psi))
    for (coeff, P) in H_terms
        result .+= coeff .* (P * psi)
    end
    return result
end


# =============================================================================
# Compact Pauli backend
#
# Bit-mask representation of a Pauli string that applies exp(-i c P dt) directly
# to amplitude pairs of a dense state vector, without ever forming a 2^n x 2^n
# sparse matrix.
#
# The kernels here are single-threaded: a single rotation touches only ~n
# amplitudes, too little work to amortize a per-call thread fork/join. Threading
# is better applied at a coarser grain by the caller (e.g. evolving U and U'
# concurrently in the Arnoldi matvec).
#
# This backend coexists with the sparse-matrix routines above; the methods are
# selected by dispatch on `CompactPauliTerm`.
# =============================================================================

"""
Compact representation of `coeff * P`, where P = ⊗_j σ_{p_j}.

- `xmask`: bit set for X or Y
- `zmask`: bit set for Z or Y
- `yphase`: i^(number of Y operators)

For computational-basis state |j⟩:

    P|j⟩ = yphase * (-1)^popcount(zmask & j) * |j ⊻ xmask⟩
"""
struct CompactPauliTerm
    coeff::Float64
    xmask::UInt64
    zmask::UInt64
    yphase::ComplexF64
end

"""
    compact_pauli_term(pauli, coeff; leftmost_is_msb=true)

Convert a Pauli string to a `CompactPauliTerm`.

`leftmost_is_msb=true` means the first character acts on the most-significant
computational-basis bit. For this repo `OP_from_string("XI") = kron(X, I)`, so
the leftmost character is the MSB; the default matches that convention.
"""
function compact_pauli_term(
    pauli::AbstractString,
    coeff::Real;
    leftmost_is_msb::Bool=true,
)
    nqubits = length(pauli)
    if nqubits > 63
        error("CompactPauliTerm supports at most 63 qubits; received $nqubits")
    end

    xmask = UInt64(0)
    zmask = UInt64(0)
    ny = 0

    for (pos, op) in enumerate(pauli)
        bit = leftmost_is_msb ? (nqubits - pos) : (pos - 1)
        bitmask = UInt64(1) << bit
        if op == 'I'
            # nothing
        elseif op == 'X'
            xmask |= bitmask
        elseif op == 'Y'
            xmask |= bitmask
            zmask |= bitmask
            ny += 1
        elseif op == 'Z'
            zmask |= bitmask
        else
            error("Unsupported Pauli character '$op' in string '$pauli'")
        end
    end

    yphase = if ny % 4 == 0
        ComplexF64(1.0, 0.0)
    elseif ny % 4 == 1
        ComplexF64(0.0, 1.0)
    elseif ny % 4 == 2
        ComplexF64(-1.0, 0.0)
    else
        ComplexF64(0.0, -1.0)
    end

    return CompactPauliTerm(Float64(coeff), xmask, zmask, yphase)
end

"""
Scalar `s` such that `P|basis⟩ = s * |basis ⊻ xmask⟩`.

Each Z or Y qubit sitting on a set bit of `basis` contributes a factor `-1`, so
the sign is `(-1)^(number of set bits in zmask & basis)`; the Y operators also
contribute the overall `yphase = i^(#Y)`.
"""
@inline function compact_pauli_phase(term::CompactPauliTerm, basis::UInt64)::ComplexF64
    parity = isodd(count_ones(term.zmask & basis))
    return parity ? -term.yphase : term.yphase
end

"""
    apply_pauli_rotation!(psi, term::CompactPauliTerm, dt)

Apply exp(-i * term.coeff * P * dt) to `psi` in place, without forming P.

Diagonal (I/Z only) strings scale each amplitude independently. Off-diagonal
strings pair each basis index j with k = j ⊻ xmask; each pair is touched exactly
once (so the update is safe to parallelize at a coarser grain if ever needed).
"""
function apply_pauli_rotation!(
    psi::AbstractVector{ComplexF64},
    term::CompactPauliTerm,
    dt::Real,
)
    theta = term.coeff * dt
    c = cos(theta)
    alpha = -im * sin(theta)
    n = length(psi)

    # Diagonal case: only I and Z, so each amplitude scales independently.
    if term.xmask == 0
        @inbounds for idx in eachindex(psi)
            basis = UInt64(idx - 1)
            phase = compact_pauli_phase(term, basis)
            psi[idx] = (c + alpha * phase) * psi[idx]
        end
        return psi
    end

    # Off-diagonal case: xmask ≠ 0, so P sends each basis state |j⟩ to a
    # *different* state |k⟩ = |j ⊻ xmask⟩. We update the two amplitudes of each
    # pair {j, k} together, and must visit every pair exactly once.
    #
    # Enumeration trick: take the lowest set bit of xmask as a "pivot". Exactly
    # one member of each pair has that bit = 0, so iterating over the n/2 indices
    # whose pivot bit is 0 hits every pair once. We turn a counter 0..n/2-1 into
    # such an index by inserting a 0 bit at the pivot position — bits below the
    # pivot stay, bits at/above it shift up by one (so the pivot bit is always 0).
    pivot = trailing_zeros(term.xmask)
    below_pivot = pivot == 0 ? UInt64(0) : (UInt64(1) << pivot) - UInt64(1)
    npairs = n >> 1

    # For the pair (j, k): P|j⟩ = phase_j·|k⟩, and because Pauli strings are
    # Hermitian the reverse element ⟨j|P|k⟩ is conj(phase_j). With the rotation
    # exp(-i c P dt) = cos·I - i sin·P and alpha = -i sin(cθ):
    #   psi[j] ← cos·psi[j] + alpha·conj(phase_j)·psi[k]
    #   psi[k] ← cos·psi[k] + alpha·phase_j·psi[j]
    @inbounds for pair_index in 0:(npairs - 1)
        counter = UInt64(pair_index)
        j = (counter & below_pivot) | ((counter & ~below_pivot) << 1)
        k = j ⊻ term.xmask
        phase_j = compact_pauli_phase(term, j)
        ji = Int(j) + 1
        ki = Int(k) + 1
        a = psi[ji]
        b = psi[ki]
        psi[ji] = c * a + alpha * conj(phase_j) * b
        psi[ki] = c * b + alpha * phase_j * a
    end

    return psi
end

"""
    first_order_trotter_statevec(H_terms::Vector{CompactPauliTerm}, psi0, dt, nsteps)

First-order Trotter evolution using the compact backend.
"""
function first_order_trotter_statevec(
    H_terms::Vector{CompactPauliTerm},
    psi0::AbstractVector{ComplexF64},
    dt::Real,
    nsteps::Int,
)
    psi = copy(psi0)
    for _ in 1:nsteps
        for term in H_terms
            apply_pauli_rotation!(psi, term, dt)
        end
    end
    return psi
end

"""
    second_order_trotter_statevec(H_terms::Vector{CompactPauliTerm}, psi0, dt, nsteps)

Second-order (Strang) Trotter evolution using the compact backend.
"""
function second_order_trotter_statevec(
    H_terms::Vector{CompactPauliTerm},
    psi0::AbstractVector{ComplexF64},
    dt::Real,
    nsteps::Int,
)
    psi = copy(psi0)
    for _ in 1:nsteps
        for term in H_terms
            apply_pauli_rotation!(psi, term, dt / 2)
        end
        for term in Iterators.reverse(H_terms)
            apply_pauli_rotation!(psi, term, dt / 2)
        end
    end
    return psi
end

"""
    fourth_order_trotter_statevec(H_terms::Vector{CompactPauliTerm}, psi0, dt, nsteps)

Fourth-order Trotter evolution using the compact backend.
"""
function fourth_order_trotter_statevec(
    H_terms::Vector{CompactPauliTerm},
    psi0::AbstractVector{ComplexF64},
    dt::Real,
    nsteps::Int,
)
    p = 1 / (4 - 4^(1 / 3))
    psi = copy(psi0)

    function apply_s2_step!(state, delta_t)
        for term in H_terms
            apply_pauli_rotation!(state, term, delta_t / 2)
        end
        for term in Iterators.reverse(H_terms)
            apply_pauli_rotation!(state, term, delta_t / 2)
        end
        return state
    end

    for _ in 1:nsteps
        apply_s2_step!(psi, p * dt)
        apply_s2_step!(psi, p * dt)
        apply_s2_step!(psi, (1 - 4p) * dt)
        apply_s2_step!(psi, p * dt)
        apply_s2_step!(psi, p * dt)
    end
    return psi
end

"""
Add `coeff * P * psi` to `result` for a compact Pauli term.

Each source amplitude `psi[basis]` contributes to `result[basis ⊻ xmask]`.
`basis ↦ basis ⊻ xmask` is a bijection, so distinct iterations write distinct
`result` slots (safe to parallelize at a coarser grain if ever needed).
"""
function add_compact_pauli_action!(
    result::AbstractVector{ComplexF64},
    psi::AbstractVector{ComplexF64},
    term::CompactPauliTerm,
)
    @inbounds for idx in eachindex(psi)
        basis = UInt64(idx - 1)
        destination = basis ⊻ term.xmask
        phase = compact_pauli_phase(term, basis)
        result[Int(destination) + 1] += term.coeff * phase * psi[idx]
    end
    return result
end

"""
    hamiltonian_matvec(psi, H_terms::Vector{CompactPauliTerm})

Compute H|ψ⟩ = Σⱼ cⱼPⱼ|ψ⟩ using the compact backend.
"""
function hamiltonian_matvec(
    psi::AbstractVector{ComplexF64},
    H_terms::Vector{CompactPauliTerm},
)
    result = zeros(ComplexF64, length(psi))
    for term in H_terms
        add_compact_pauli_action!(result, psi, term)
    end
    return result
end
