# =============================================================================
# first_order_bond_energy_curve.jl
#
# First-order Trotter version of bond_energy_curve.jl.
#
# Usage:
#   julia --project=. first_order_bond_energy_curve.jl [options]
# =============================================================================

include("bond_energy_curve.jl")

function default_outpath(opts)
    return @sprintf(
        "first_order_bond_energy_curve_%s_as-%03d-%03d_%s.png",
        opts.basis,
        opts.electrons,
        opts.vacant,
        opts.encoding
    )
end

function trotter_order()
    return :first
end

function save_plot(points, opts)
    outpath = opts.outpath === nothing ? default_outpath(opts) : opts.outpath
    xs = [point.bond_length for point in points]
    exact = [point.exact_energy for point in points]
    trotter_1 = [point.trotter_1 for point in points]
    trotter_5 = [point.trotter_5 for point in points]

    plot(
        xs,
        exact;
        label="exact",
        marker=:circle,
        linewidth=2,
        xlabel="Li-Li bond length",
        ylabel="Ground state energy (Ha)",
        title=@sprintf(
            "Li-Li %s as-%03d-%03d %s first-order",
            opts.basis,
            opts.electrons,
            opts.vacant,
            uppercase(opts.encoding)
        ),
        legend=:best
    )
    plot!(xs, trotter_1; label="First-order 1 step", marker=:square, linewidth=2)
    plot!(xs, trotter_5; label="First-order 5 steps", marker=:diamond, linewidth=2)
    savefig(outpath)
    return outpath
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
