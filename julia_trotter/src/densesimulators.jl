include("tensorsimulators.jl")
using SparseArrays, LinearAlgebra, SpecialFunctions
dense_matrix_exp(H,t)  = exp(-im*Matrix(H)*t)

"""
These functions when called return the full U amtrix after expenentiation/trotterization
"""

function chebyshev_exponentiation(H, order, t, energy_bound; α::Int = 4)

        # materialize Hermitian/Symmetric views to full sparse CSC to avoid CHOLMOD path
        # For Hermitian matrices, convert to sparse directly (symmetry is already handled)
    Hs = H isa Hermitian ? sparse(H) : H
    n  = size(Hs, 1)

    scale = ceil(Int, t * energy_bound)
    dt    = t / scale

    # make an explicit sparse identity with matching eltype
    state = spdiagm(0 => ones(eltype(Hs), n))

    coeffs = [2*(im^n)*besselj(n, 1.0/α ) for n in 1:order]
    coeffs = [besselj(0, 1.0/α), coeffs ...]
    for i in 1:(Int(scale))*α
        cheby_list = [state, -Hs*dt*state]
        for  k =2:order
            push!(cheby_list, (-2.0*dt)*Hs*cheby_list[end] - cheby_list[end-1])
        end
        state = sum(coeffs[i]*cheby_list[i] for i in 1:order+1)
        # push!(Obs_vals, [ state'*X*state for X in Obs])
    end

    error1 =  sum(2*abs(besselj(n, 1.0/α )) for n in order+1:1000)*scale*α
    error2 =  sum(2*abs(besselj(n, 1.0/α )) for n in order+1:4000)*scale*α
    @show error1, error2
    @assert abs(error1) < 1E-8
    @info "In those error values are close, then the method has converged to that error bound"
    return state
end


function dense_first_order_trotter(H_list, dt::Float64)
    return mapreduce(x -> dense_matrix_exp(x, -dt), *, H_list)
end


function dense_second_order_trotter(H_list, dt::Float64)
    U_temp1 = mapreduce(x->dense_matrix_exp(x, -dt/2), *, H_list[1:end-1] )
    U_temp2 = mapreduce(x->dense_matrix_exp(x, -dt/2), *, reverse(H_list[1:end-1]))
    return U_temp1 * dense_matrix_exp(H_list[end],-dt) * U_temp2
end


function dense_fourth_order_trotter(H_list, dt::Float64)
    u_4 = 1/(4 - 4^(1/3))
    U1 = dense_second_order_trotter(H_list, dt*u_4) 
    U2 = dense_second_order_trotter(H_list, dt*(1- 4*u_4))
    U1sq = U1*U1
    return U1sq*U2*U1sq
end




