using ITensors
using Plots
using SpecialFunctions
include("quantum_utils.jl")

compute_observable(state::ITensor, Op) = Array(conj(state)*noprime(Op*state))[1]
    ##Just a function for ease of use



function chebyshev_simulation(ψ, T::Real, Ham_mult::Function, energy_bound::Real, Observables; α=3, order = 5, returnstate=false, normalize = false)
    """
    Implements exp(-iHT)|ψ⟩ based on the Jacobi-Angers fromula and the following property of chebyshev polynomials (T_n(cos(x)) = cos(nx)).
energy_bound:an upper bound on the spectral norm (largest abs eigvenval) of the Hamiltonian.
α: an aditional factor that can increased to improve the accuracy of the expansion. This will cause an O(α) increase in the runtime.
order: the maximum order in the polynomial expansion

    """

    @assert α  >= 1
    scale = ceil(T*energy_bound)
    dt = T/scale
    coeffs = [2*(im^n)*besselj(n, 1.0/α ) for n in 1:order]
    coeffs = [besselj(0, 1.0/α), coeffs ...]
    global statetemp = ψ
    Obs_dict = Dict(k => compute_observable(statetemp,V) for (k,V) in Observables)

    for i in 1:(Int(scale))*α
        cheby_list = [statetemp , Ham_mult(-dt*statetemp)]
        for  k =2:order
            ##Computing chebyshev polynomials by recusrion
            push!(cheby_list, (-2.0*dt)*Ham_mult(cheby_list[end]) - cheby_list[end-1])
        end
        global statetemp = sum(coeffs[i]*cheby_list[i] for i in 1:order+1)

        if normalize == true
          global statetemp = statetemp/norm(statetemp)
        end

        for k in keys(Obs_dict)
          Obs_dict[k] =  compute_observable(statetemp, Observables[k])
        end

    end

    error1 =  sum(2*abs(besselj(n, 1.0/α )) for n in order+1:4000)*scale*α
    error2 =  sum(2*abs(besselj(n, 1.0/α )) for n in order+1:8000)*scale*α
    @show error1, error2
    @assert error1 < 1e-7
    @assert abs(error1 -error2) < 1e-7
    @info "In those error values are close, then the method has converged to that error bound"

    if returnstate == false
        return Obs_dict, (error = error2, dt = dt)
    else
        return statetemp, Obs_dict, (error = error2, dt = dt)
    end
end

 


function chebyshev_simulation(ψ::ITensor, T::Real, H_op_list::Vector{ITensor}, energy_bound::Real, Observables; α=3, order = 5, returnstate=false, normalize = false)
    """
    Implements exp(-iHT)|ψ⟩ based on the Jacobi-Angers fromula and the following property of chebyshev polynomials (T_n(cos(x)) = cos(nx)).
energy_bound:an upper bound on the spectral norm (largest abs eigvenval) of the Hamiltonian.
α: an aditional factor that can increased to improve the accuracy of the expansion. This will cause an O(α) increase in the runtime.
order: the maximum order in the polynomial expansion

    """

    @assert α  >= 1
    Ham_mult(x) = apply_H(x, H_op_list)
    scale = ceil(T*energy_bound)
    dt = T/scale
    coeffs = [2*(im^n)*besselj(n, 1.0/α ) for n in 1:order]
    coeffs = [besselj(0, 1.0/α), coeffs ...]
    global statetemp = ψ
    Obs_dict = Dict(k => compute_observable(statetemp,V) for (k,V) in Observables)

    for i in 1:(Int(scale))*α
        cheby_list = [statetemp , Ham_mult(-dt*statetemp)]
        for  k =2:order
            ##Computing chebyshev polynomials by recusrion
            push!(cheby_list, (-2.0*dt)*Ham_mult(cheby_list[end]) - cheby_list[end-1])
        end
        global statetemp = sum(coeffs[i]*cheby_list[i] for i in 1:order+1)

        if normalize == true
          global statetemp = statetemp/norm(statetemp)
        end

        for k in keys(Obs_dict)
          Obs_dict[k] =  compute_observable(statetemp, Observables[k])
        end

    end

    error1 =  sum(2*abs(besselj(n, 1.0/α )) for n in order+1:4000)*scale*α
    error2 =  sum(2*abs(besselj(n, 1.0/α )) for n in order+1:8000)*scale*α
    @show error1, error2
    @assert error1 < 1e-7
    @assert abs(error1 -error2) < 1e-7
    @info "In those error values are close, then the method has converged to that error bound"

    if returnstate == false
        return Obs_dict, (error = error2, dt = dt)
    else
        return statetemp, Obs_dict, (error = error2, dt = dt)
    end
end

function apply_H(x::ITensor, H_op_list)
    """
    Compute H |x⟩ whenre H is given as a list of local operators.
    """
    state = noprime(x)
    state = noprime(H_op_list[1]*state)
    for O  in H_op_list[2:end]
        state = state + noprime(O*x)
    end
    return state
  end


function apply_trotter(ψ::ITensor, U_list)
    """Applies one step of second  order trotter formula.
U_list is an array
    """
    for U in U_list
        for R in U
            ψ = noprime(R*ψ)
        end
    end
    return ψ
end

function second_order_trotter( ψ::ITensor,T, nsteps, H_blocks_list, Observables; returnstate = false)
    
    dt = T/nsteps
    R_ops = [[exp(-im*O*dt/2) for O in X] for X in H_blocks_list[1:end-1]]
    R_ops = vcat(R_ops...)
    R_ops_end = [exp(-im*O*dt) for O in H_blocks_list[end]]
    global state = ψ
    Obs_dict = Dict(k => compute_observable(state,V) for (k,V) in Observables)


    for tind in 1:nsteps
        state = apply_trotter(state, [R_ops, R_ops_end, reverse(R_ops)])
        global current_time = tind*dt
    end
    @assert current_time == T

    for k in keys(Obs_dict)
      Obs_dict[k] = compute_observable(state, Observables[k])
    end


    if returnstate == false
        return Obs_dict, (dt = dt,)
    else
        return state, Obs_dict, (dt = dt,)
    end
end

function fourth_order_trotter(ψ::ITensor,T, nsteps, H_blocks_list, Observables; returnstate = false)
    """
Computed from the recursive fromula given in Child's paper on the theory of commutator error scaling
"""
    dt = T/nsteps
    u_4 = 1/(4 - 4^(1/3))

    R_ops_1 = [[exp(-im*O*dt*u_4/2) for O in X] for X in H_blocks_list[1:end-1]]
    R_ops_1 = vcat(R_ops_1...)
    R_ops_end_1 = [exp(-im*O*dt*u_4) for O in H_blocks_list[end]]

    R_ops_2 = [[exp(-im*O*dt*(1-4*u_4)/2) for O in X] for X in H_blocks_list[1:end-1]]
    R_ops_2 = vcat(R_ops_2...)
    R_ops_end_2 = [exp(-im*O*dt*(1-4*u_4)) for O in H_blocks_list[end]]

    global state = ψ
    Obs_dict = Dict(k => compute_observable(state,V) for (k,V) in Observables)
    for tind in 1:nsteps
        ##XXX This is suboptimal as the some of the circuit elements can be combined.

        state = apply_trotter(state, [R_ops_1, R_ops_end_1, reverse(R_ops_1)])
        state = apply_trotter(state, [R_ops_1, R_ops_end_1, reverse(R_ops_1)])
        state = apply_trotter(state, [R_ops_2, R_ops_end_2, reverse(R_ops_2)])
        state = apply_trotter(state, [R_ops_1, R_ops_end_1, reverse(R_ops_1)])
        state = apply_trotter(state, [R_ops_1, R_ops_end_1, reverse(R_ops_1)])
   #     global current_time = tind*dt
    end
   # @assert current_time ≈  T
    for k in keys(Obs_dict)
        Obs_dict[k] = compute_observable(state, Observables[k])
    end



    if returnstate == false
        return Obs_dict, (dt = dt,)
    else
        return state, Obs_dict, (dt = dt,)
    end
end




function fourth_order_trotter_tfim(ψ::ITensor,T, nsteps, ZZ_ops, X_ops, Observables; returnstate = false)
    """
Computed from the recursive fromula given in Child's paper on the theory of commutator error scaling
"""
    dt = T/nsteps
    u_4 = 1/(4 - 4^(1/3))
    RZZ_ops1 = [exp(-im*T*dt*u_4/2) for T in ZZ_ops]
    RZZ_ops2 = [exp(-im*T*dt*(1-4*u_4)/2) for T in ZZ_ops]
    RX_ops1 = [exp(-im*T*dt*u_4) for T in X_ops]
    RX_ops2 = [exp(-im*T*dt*(1-4*u_4)) for T in X_ops]
    
    Obs_list = []
    global state = ψ 
    for tind in 1:nsteps
        ##XXX This is suboptimal as the ZZ roations can be combined.
        ###Below I am blindly applying the recusiuon
        ##TODO Faster implementation
        state = apply_trotter(state, [RZZ_ops1, RX_ops1, RZZ_ops1])
        state = apply_trotter(state, [RZZ_ops1, RX_ops1, RZZ_ops1])
        state = apply_trotter(state, [RZZ_ops2, RX_ops2, RZZ_ops2])
        state = apply_trotter(state, [RZZ_ops1, RX_ops1, RZZ_ops1])
        state = apply_trotter(state, [RZZ_ops1, RX_ops1, RZZ_ops1])
        current_time = tind*dt
        push!(Obs_list, [conj(state)*noprime(X*state) for X in Observables])
    end

    Obs_list = [[Array(x)[1] for x in y] for y in Obs_list]
    Obs_array = hcat(Obs_list...)'

    if returnstate == false
        return Obs_array, (dt = dt,)
    else
        return state, Obs_array, (dt = dt,)
    end
end









function second_order_trotter_tfim(ψ::ITensor,T, nsteps, ZZ_ops, X_ops, Observables; returnstate = false)
    dt = T/nsteps
    RZZ_ops = [exp(-im*T*dt/2) for T in ZZ_ops]
    RX_ops = [exp(-im*T*dt) for T in X_ops]
    Obs_list = []
    global state = ψ 

    for tind in 1:nsteps
        state = apply_trotter(state, [RZZ_ops, RX_ops, RZZ_ops])
        current_time = tind*dt
        push!(Obs_list, [conj(state)*noprime(X*state) for X in Observables])
    end

    Obs_list = [[Array(x)[1] for x in y] for y in Obs_list]
    Obs_array = hcat(Obs_list...)'

    if returnstate == false
        return Obs_array, (dt = dt,)
    else
        return state, Obs_array, (dt = dt,)
    end
end




function apply_H(x, H_op_list)
    state = noprime(x)
    for O  in H_op_list
        state = state + noprime(O*x)
    end
    return state
end

function apply_trotter(ψ::ITensor, U_list)
    """Applies one step of seonc order trotter formula.
U_list is an array
    """
    for U in U_list
        for R in U
            ψ = noprime(R*ψ)
        end
    end
    return ψ
end

