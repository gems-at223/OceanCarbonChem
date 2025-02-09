#import Pkg; Pkg.add("GLMakie")

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!, ImpenetrableBoundaryCondition
using Printf
using GLMakie
grid = RectilinearGrid(size=128, z=(-10, 10), topology=(Flat, Flat, Bounded))
sinking = AdvectiveForcing(u=0,v=0, w=-0.2)
rising = AdvectiveForcing(u=0,v=0, w=+0.2)

k_1_plus=0.0036
k_1_neg=3.3*10^4
k_2_plus=0
k_2_neg=0
k_3_plus=0
k_3_neg=0
k_4_plus=8500
k_4_neg=3.7*10^-4
k_5_plus=10^10
k_5_neg=9
k_6_neg=2.8*10^10
k_6_plus=1.3*10^-3
k_7_neg=10^10
k_7_plus=20


CO2_ev( z, t, H,HCO3,OH,CO2) = ((k_1_neg*H+k_4_neg)*HCO3) - (k_1_plus -(k_4_plus*OH))*CO2

HCO3_ev( z, t, CO2,H,HCO3,OH,Carbonate) = (k_1_plus*CO2) -(k_1_neg*H*HCO3) +(k_4_plus*CO2*OH) - (k_4_neg*HCO3) +(k_5_plus*H*Carbonate) -(k_5_neg*HCO3)

Carbonate_ev(z,t,HCO3,H,Carbonate)=(k_5_neg*HCO3) - (k_5_plus*H*Carbonate)

H_ev(z,t,H,HCO3,CO2,Carbonate,OH,BOH3,BOH4)=((k_5_neg-(k_1_neg*H))*HCO3) + (k_1_plus*CO2) - (k_5_plus*H*Carbonate) - k_6_plus -(k_6_neg*H*OH) +(k_7_plus*BOH3) - (k_7_neg*H*BOH4)

OH_ev(z,t,HCO3,CO2,OH,H)=(k_4_neg*HCO3) - (k_4_plus*CO2*OH)+k_6_plus-(k_6_neg*H*OH)

BOH3_ev(z,t,BOH3,BOH4,H)= -(k_7_plus*BOH3)+ (k_7_neg*H*BOH4)

BOH4_ev(z,t,BOH3,BOH4,H)= (k_7_plus*BOH3)- (k_7_neg*H*BOH4)

CO2_reaction = Forcing(CO2_ev, field_dependencies=(:H, :HCO3,:OH,:CO2))
HCO3_reaction = Forcing(HCO3_ev, field_dependencies=(:CO2, :H,:HCO3,:OH,:Carbonate))
Carbonate_reaction=Forcing(Carbonate_ev,field_dependencies=(:HCO3,:H,:Carbonate))
H_reaction=Forcing(H_ev,field_dependencies=(:H,:HCO3,:CO2,:Carbonate,:OH,:BOH3,:BOH4))
OH_reaction=Forcing(OH_ev,field_dependencies=(:HCO3,:CO2,:OH,:H))
BOH3_reaction=Forcing(BOH3_ev,field_dependencies=(:BOH3,:BOH4,:H))
BOH4_reaction=Forcing(BOH4_ev,field_dependencies=(:BOH3,:BOH4,:H))


model = NonhydrostaticModel(; grid,
                                    velocities = nothing,
                                    tracers = (:CO2, :HCO3,:Carbonate,:H,:OH,:BOH3,:BOH4),
                                    buoyancy = nothing,
                                    closure = ScalarDiffusivity(κ=1e-4),
                                    forcing = (;CO2=(CO2_reaction),HCO3 = (HCO3_reaction),Carbonate=(Carbonate_reaction),
                                    H=(H_reaction),OH=(OH_reaction),BOH3=(BOH3_reaction),BOH4=(BOH4_reaction)
                                    ))


CO2ᵢ(z) = exp(-((z - 3)^2) / 2)  # Wider spread (dividing by 4)
HCO3ᵢ(z) = exp(-((z + 3)^2) / 2)
Carbonateᵢ(z) = exp(-((z - 3)^2) / 2)  # Wider spread (dividing by 4)
Hᵢ(z) = exp(-((z + 6)^2) / 2)
OHᵢ(z) = exp(-((z - 3)^2) / 2)  # Wider spread (dividing by 4)
BOH3ᵢ(z) = exp(-((z + 6)^2) / 2)
BOH4ᵢ(z) = exp(-((z + 6)^2) / 2)

#aᵢ(z) = sin(π * z)
#bᵢ(z) = cos(π * z)
set!(model, CO2=CO2ᵢ, HCO3=HCO3ᵢ,Carbonate=Carbonateᵢ,H=Hᵢ,OH=OHᵢ,BOH3=BOH3ᵢ,BOH4=BOH4ᵢ)

simulation = Simulation(model; Δt=1e-9, stop_iteration=100)
progress(sim) = @info @sprintf("Iter: %d, time: %.2e", iteration(sim), time(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

z = znodes(grid, Center())
CO2 = interior(model.tracers.CO2, 1, 1, :)
HCO3 = interior(model.tracers.HCO3, 1, 1, :)

fig = Figure()
ax = Axis(fig[1, 1],
          xlabel = "Reactant concentration",
          ylabel = "z",
          title = "Tracer reactions at t=0")
xlims!(ax, -1, 4)

ℓa = lines!(ax, CO2, z, label="a")
ℓb = lines!(ax, HCO3, z, label="b")

axislegend(ax, position=:rb)

display(fig)

function update_plot!(sim)
    ℓa[1] = interior(sim.model.tracers.CO2, 1, 1, :)
    ℓb[1] = interior(sim.model.tracers.HCO3, 1, 1, :)
    ax.title[] = @sprintf("Tracer reactions at t=%.2e", time(sim))
end

record(fig, "tracer_reactions.mp4", 1:400, framerate=24) do nn
    simulation.stop_iteration += 5
    run!(simulation)
    update_plot!(simulation)
end