
include("QSS.jl")
using OceanBioME, Oceananigans
using Oceananigans.Units

grid = RectilinearGrid(CPU(), size = (128, 32), extent = (100meters, 500meters), topology = (Bounded, Flat, Bounded))

biogeochemistry = DynamicCarbonateChemistry()


forcing=AdvectiveForcing(w=1)
model = NonhydrostaticModel(; grid, biogeochemistry,
                              advection = WENO(; grid),
                              tracers = (:T, :c₁, :c₂, :c₃, :c₅, :c₆, :c₇),  
                              closure = ScalarDiffusivity(κ=1e-3),
                              buoyancy=SeawaterBuoyancy(constant_salinity = true))
                             # forcing=(;c₃=(forcing),c₁ = (forcing)))

@inline front(x, z, μ, δ) = μ + δ * tanh((x - 1000 + 4 * z) / 500)

Tᵢ(x, z) = front(x, z, 9, 0.05)

c1ᵢ(x, z) = ifelse(z > -100, 7.57*10^-6, 3.57*10^-6)
c2iᵢ(x, z) = ifelse( z > -100, 1.67*10^-3, 1.67*10^-2)
c3ᵢ(x, z) = ifelse( z > -100, 3.15*10^-4, 1.15*10^-3)
#c4ᵢ(x, z) = front(x, z, 0.0000001, 0.0000001)
c5ᵢ(x, z) = ifelse( z > -100, 9.60*10^-6, 9.60*10^-6)
c6ᵢ(x, z) = ifelse( z > -100, 2.97e-4, 1.97e-4)
c7ᵢ(x, z) = ifelse(z > -100, 1.19e-4, 1.19e-4)


set!(model,T=Tᵢ, c₁=c1ᵢ, c₂=c2iᵢ, c₃=c3ᵢ, c₅=c5ᵢ, c₆=c6ᵢ, c₇=c7ᵢ)

simulation = Simulation(model; Δt = 0.003seconds, stop_time = 0.5seconds)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "buoyancy_front.jld2",
                                                       schedule = TimeInterval(0.01seconds),
                                                       overwrite_existing = true)

run!(simulation)

T = FieldTimeSeries("buoyancy_front.jld2", "c₁")
N = FieldTimeSeries("buoyancy_front.jld2", "c₂")
P = FieldTimeSeries("buoyancy_front.jld2", "c₃")

xc, yc, zc = nodes(T)

times = T.times

using CairoMakie

n = Observable(1)

T_lims = (6.57e-6, 7.5699999999999995e-6)
N_lims = (0, 0.02)
P_lims = ( 0.000315,0.00115)

Tₙ = @lift interior(T[$n], :, 1, :)
Nₙ = @lift interior(N[$n], :, 1, :)
Pₙ = @lift interior(P[$n], :, 1, :)

fig = Figure(size = (1000, 520), fontsize = 20)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title)

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 770, yticks = [-400, -200, 0])
ax1 = Axis(fig[1, 1]; title = "Temperature (°C)", axis_kwargs...)
ax2 = Axis(fig[2, 1]; title = "Nutrients concentration (mmol N / m³)",axis_kwargs...)
ax3 = Axis(fig[3, 1]; title = "Phytoplankton concentration (mmol N / m³)", axis_kwargs...)

hm1 = heatmap!(ax1, xc, zc, Tₙ, colorrange = T_lims, colormap = Reverse(:lajolla), interpolate = true)
hm2 = heatmap!(ax2, xc, zc, Nₙ, colorrange = N_lims, colormap = Reverse(:bamako), interpolate = true)
hm3 = heatmap!(ax3, xc, zc, Pₙ, colorrange = P_lims, colormap = Reverse(:bamako), interpolate = true)

Colorbar(fig[1, 2], hm1, ticks = [6.57e-6, 7.5699999999999995e-6])
#Colorbar(fig[2, 2], hm2, ticks = [0, 2, 4])
Colorbar(fig[3, 2], hm3, ticks = [0.000315,0.00115])

rowgap!(fig.layout, 0)

record(fig, "buoyancy_front.gif", 1:length(times)) do i
    n[] = i
end
P_min = minimum(interior(P[1], :, 1, :))
P_max = maximum(interior(P[1], :, 1, :))
println("P concentration range: $P_min to $P_max")

T_min = minimum(interior(T[1], :, 1, :))
T_max = maximum(interior(T[1], :, 1, :))
println("P concentration range: $T_min to $T_max")

N_min = minimum(interior(N[1], :, 1, :))
N_max = maximum(interior(N[1], :, 1, :))
println("P concentration range: $N_min to $N_max")
