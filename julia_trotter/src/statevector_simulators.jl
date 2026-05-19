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
    U(T) ≈ [exp(-iH₁dt) exp(-iH₂dt) ⋯ exp(-iHₘdt)]^nsteps

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
