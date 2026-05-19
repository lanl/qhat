##The aim here is to define a trotterization framework that onbly cares about how things mulitply.
##This way input datastarcutes can be changed as long as mulitplicaitons are correctly defined

function second_order_trotter_step(ψ::ITensor,dt,Layers::Vector{F}) where {F <: Function}
    """
Layer[i] is a function that takes state, timestep and applies everything in that layer by time dt
    """
    depth = length(Layers)

    for ii in 1:depth-1
        ψ = Layers[ii](ψ,dt/2)
    end
    ψ = Layers[end](ψ,dt)

    for ii in (depth-1):-1:1
        ψ = Layers[ii](ψ,dt/2)
    end
    return ψ 
end


function second_order_trotter(ψ::ITensor,T::Real,nsteps::Int,Layers::Vector{F}) where {F <: Function}
    dt = T/nsteps
    for t in 1:nsteps
        ψ= second_order_trotter_step(ψ, dt, Layers)
    end

    return ψ 
    end

function fourth_order_trotter(ψ::ITensor,T::Real,nsteps::Int,Layers::Vector{F}) where {F <: Function}
    dt = T/nsteps
    for t in 1:nsteps
        p = inv(4 - (4^(1/3)))
        ψ= second_order_trotter_step(ψ, p*dt, Layers)
        ψ= second_order_trotter_step(ψ, p*dt, Layers)
        ψ= second_order_trotter_step(ψ, (1-4*p)*dt, Layers)
        ψ= second_order_trotter_step(ψ, p*dt, Layers)
        ψ= second_order_trotter_step(ψ, p*dt, Layers)
    end

    return ψ 
end



