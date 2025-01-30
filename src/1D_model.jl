#import Pkg; Pkg.add("GLMakie")

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!, ImpenetrableBoundaryCondition
using Printf
using GLMakie
grid = RectilinearGrid(size=128, z=(-10, 10), topology=(Flat, Flat, Bounded))
sinking = AdvectiveForcing(u=0,v=0, w=-0.09)
rising = AdvectiveForcing(u=0,v=0, w=+0.09)

b_to_a( z, t, a, b) = + a * b*0.2
a_to_b( z, t, a, b) = - a * b*0.1

a_reaction = Forcing(a_to_b, field_dependencies=(:a, :b))
b_reaction = Forcing(b_to_a, field_dependencies=(:a, :b))


model = NonhydrostaticModel(; grid,
                                    velocities = nothing,
                                    tracers = (:a, :b),
                                    buoyancy = nothing,
                                    closure = ScalarDiffusivity(κ=1e-3),
                                    forcing = (;a=(a_reaction,sinking),b = (b_reaction,rising)))

println(model.forcing.a)
println(model.forcing.b)

#aᵢ(z) = exp(-((z - 3)^2) / 2)  # Wider spread (dividing by 4)
#bᵢ(z) = exp(-((z + 6)^2) / 2)
aᵢ(z) = sin(π * z)
bᵢ(z) = cos(π * z)
set!(model, a=aᵢ, b=bᵢ)

simulation = Simulation(model; Δt=1e-2, stop_iteration=1500)
progress(sim) = @info @sprintf("Iter: %d, time: %.2e", iteration(sim), time(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

z = znodes(grid, Center())
a = interior(model.tracers.a, 1, 1, :)
b = interior(model.tracers.b, 1, 1, :)

fig = Figure()
ax = Axis(fig[1, 1],
          xlabel = "Reactant concentration",
          ylabel = "z",
          title = "Tracer reactions at t=0")
xlims!(ax, -1, 4)

ℓa = lines!(ax, a, z, label="a")
ℓb = lines!(ax, b, z, label="b")

axislegend(ax, position=:rb)

display(fig)

function update_plot!(sim)
    ℓa[1] = interior(sim.model.tracers.a, 1, 1, :)
    ℓb[1] = interior(sim.model.tracers.b, 1, 1, :)
    ax.title[] = @sprintf("Tracer reactions at t=%.2e", time(sim))
end

record(fig, "tracer_reactions.mp4", 1:400, framerate=24) do nn
    simulation.stop_iteration += 10
    run!(simulation)
    update_plot!(simulation)
end