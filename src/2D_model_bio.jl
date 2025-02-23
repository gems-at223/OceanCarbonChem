
include("DynamicCarbonateChemistry/QSS.jl")
using OceanBioME, Oceananigans
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.BoundaryConditions

grid = RectilinearGrid(CPU(), size = (24,24),x = (-0.025meters, 0.025meters), z=(-0.025meters, 0.025meters), topology = (Bounded, Flat, Bounded))
biogeochemistry = DynamicCarbonateChemistry()

horizontal_closure = HorizontalScalarDiffusivity(ν=1e-3, κ=2e-9)
vertical_closure = VerticalScalarDiffusivity(ν=1e-3, κ=2e-9)

CO₂_flux = 1e-5  # Example flux value
R=0.007



@inline function inside_cylinder(x, z)
    return ((x)^2 + (z)^2) <= R^2
end


underlying_grid=grid
immersed_grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(inside_cylinder))

boundary_conditions = (c₃ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₂ = FieldBoundaryConditions(top = FluxBoundaryCondition(CO₂_flux), bottom = GradientBoundaryCondition(0.0)),
                       c₁ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₄ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₅ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₆ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₇ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)))
#forcing=AdvectiveForcing(w=-0.0000001)
#grid=ImmersedBoundaryGrid(bgrid, boundary)

model = NonhydrostaticModel(; grid=immersed_grid, biogeochemistry,
                             #advection=WENO(;grid),
                              tracers = ( :c₁, :c₂, :c₃, :c₅, :c₆, :c₇),
                              boundary_conditions = boundary_conditions,
                              #closure = (horizontal_closure, vertical_closure)
                              closure = ScalarDiffusivity(ν=1e-3, κ=5e-6))

@inline front(x, z, μ, δ) = μ + δ * tanh((x - 0.0001 + 0.00004 * z) / 0.00005)

#Tᵢ(x, z) = front(x, z, 9, 0.05)

c1ᵢ(x, z) = front(x, z, 7.57*10^-6, 3.57*10^-6)#7.57*10^-6  
c2i(x, z) = front(x,z,1.67*10^-3, (1.67*10^-3)/2) # 1.67*10^-3  
c3ᵢ(x, z) =  front(x,z,3.15*10^-4, (3.15*10^-4)/2) #3.15*10^-4
#C4ᵢ(x, z) = 6.31*10^-9
c5ᵢ(x, z) =  9.60*10^-6
c6ᵢ(x, z) =  2.97e-4
c7ᵢ(x, z) =  1.19e-4


set!(model, c₁=c1ᵢ, c₂=c2i, c₃=c3ᵢ, c₅=c5ᵢ, c₆=c6ᵢ, c₇=c7ᵢ)

simulation = Simulation(model; Δt = 0.0005seconds, stop_time = 1seconds)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "buoyancy_front.jld2",
                                                       schedule = TimeInterval(0.0005*50seconds),
                                                       overwrite_existing = true)

run!(simulation)

T = FieldTimeSeries("buoyancy_front.jld2", "c₁")
N = FieldTimeSeries("buoyancy_front.jld2", "c₂")
P = FieldTimeSeries("buoyancy_front.jld2", "c₃")

xc, yc, zc = nodes(T)

times = T.times

using CairoMakie

n = Observable(1)

T_lims = (0.0 , 1.1139999999999998e-5)
N_lims = ( 0.0 ,0.0025050000000000003)
P_lims = ( 0.0 , 0.00047250000000000005)

Tₙ = @lift interior(T[$n], :, 1, :)
Nₙ = @lift interior(N[$n], :, 1, :)
Pₙ = @lift interior(P[$n], :, 1, :)

fig = Figure(size = (1000, 1000), fontsize = 20)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title)

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 500, yticks = [-0.025,-0.01,0,0.01,0.025])
ax1 = Axis(fig[1, 1]; title = "CO2 mol/kg", axis_kwargs...)
ax2 = Axis(fig[2, 1]; title = "HCO3- mol/kg",axis_kwargs...)
ax3 = Axis(fig[3, 1]; title = "CO3-2 mol/kg", axis_kwargs...)

hm1 = heatmap!(ax1, xc, zc, Tₙ, colorrange = T_lims, colormap = Reverse(:lajolla), interpolate = true)
hm2 = heatmap!(ax2, xc, zc, Nₙ, colorrange = N_lims, colormap = Reverse(:bamako), interpolate = true)
hm3 = heatmap!(ax3, xc, zc, Pₙ, colorrange = P_lims, colormap = Reverse(:bamako), interpolate = true)

Colorbar(fig[1, 2], hm1, ticks = [4.144952723679106e-6 , 7.347076171977047e-6])
Colorbar(fig[2, 2], hm2, ticks = [0.0008689035081994547 , 0.0016178595528293655])
Colorbar(fig[3, 2], hm3, ticks = [ 0.00016389497310348997 , 0.00030516512523428153])

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
