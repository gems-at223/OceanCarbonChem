include("DynamicCarbonateChemistry/QSS.jl")


using OceanBioME, Oceananigans
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.BoundaryConditions

grid = RectilinearGrid(CPU(), size = (16,16),x = (-0.001meters, 0.001meters), z=(-0.001meters, 0.001meters), topology = (Bounded, Flat, Bounded))
biogeochemistry = DynamicCarbonateChemistry()

horizontal_closure = HorizontalScalarDiffusivity(ν=1e-3, κ=2e-9)
vertical_closure = VerticalScalarDiffusivity(ν=1e-3, κ=2e-9)

CO₂_flux = 4e-6  # Example flux value
R=0.0002



@inline function inside_cylinder(x, z)
    return ((x)^2 + (z)^2) <= R^2
end
# Calculate the distance from the center
@inline function distance_from_center(x, z)
    return sqrt(x^2 + z^2)
end

# Forcing function that only applies outside the sphere
@inline function CO2_forcing(x, z,c₁, max_distance=R)
    # Only apply forcing if outside the sphere (cylinder)
    if !inside_cylinder(x, z)  # Forcing is only applied outside the sphere
        # Calculate the distance from the center of the sphere
        dist = distance_from_center(x, z)
        return - ((3.53*10^-12/(0.0003*2*π*dist))*(c₁/((5*10^-6)+c₁)))
    else
        return 0.0  # No forcing inside the sphere
    end
end

@inline function HCO3_forcing(x, z,c₁, max_distance=R)
    # Only apply forcing if outside the sphere (cylinder)
    if !inside_cylinder(x, z)  # Forcing is only applied outside the sphere
        # Calculate the distance from the center of the sphere
        dist = distance_from_center(x, z)
        #println(- ((3.53*10^-12/(0.002*4*π*dist^2))*(1-(c₁/((5*10^-6)+c₁)))))
        return - ((3.53*10^-12/(0.0003*2*π*dist))*(1-(c₁/((5*10^-6)+c₁))))
    else
        return 0.0  # No forcing inside the sphere
    end
end

@inline function OH_forcing(x, z,c₁, max_distance=R)
    # Only apply forcing if outside the sphere (cylinder)
    if !inside_cylinder(x, z)  # Forcing is only applied outside the sphere
        # Calculate the distance from the center of the sphere
        dist = distance_from_center(x, z)
        return  ((3.53*10^-12/(0.0009*4*π*dist^2))*(1-(c₁/((5*10^-6)+c₁))))
    else
        return 0.0  # No forcing inside the sphere
    end
end


 forcing = Forcing(CO2_forcing, field_dependencies=(:c₁))
 forcinghco3=Forcing(HCO3_forcing, field_dependencies=(:c₁))
 OHforcing=Forcing(OH_forcing, field_dependencies=(:c₁))

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
                             # boundary_conditions = boundary_conditions,
                              #closure = (horizontal_closure, vertical_closure)
                              closure = ScalarDiffusivity(ν=1e-3, κ=2e-9),
                              forcing = (;c₁=(forcing),c₂=(forcinghco3),c₅=(OHforcing)))


@inline front(x, z, μ, δ) = μ + δ * tanh((x - 0.01 + 0.00004 * z) / 0.005)

#Tᵢ(x, z) = front(x, z, 9, 0.05)

c1ᵢ(x, z) = 7.57*10^-6  #front(x, z, 7.57*10^-6, 3.57*10^-6) #
c2i(x, z) =2*1.67*10^-3 #front(x,z,1.67*10^-4, (1.67*10^-3))   
c3ᵢ(x, z) = 3.15*10^-4 #front(x,z,3.15*10^-4, (3.15*10^-4)/2)    
#C4ᵢ(x, z) = 6.31*10^-9
c5ᵢ(x, z) =  9.60*10^-6
c6ᵢ(x, z) =  2.97e-4
c7ᵢ(x, z) =  1.19e-4


set!(model, c₁=c1ᵢ, c₂=c2i, c₃=c3ᵢ, c₅=c5ᵢ, c₆=c6ᵢ, c₇=c7ᵢ)

simulation = Simulation(model; Δt = 0.001seconds, stop_time = 2seconds)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "buoyancy_front.jld2",
                                                       schedule = TimeInterval(0.001*10seconds),
                                                       overwrite_existing = true)

run!(simulation)

CO2 = FieldTimeSeries("buoyancy_front.jld2", "c₁")
HCO3 = FieldTimeSeries("buoyancy_front.jld2", "c₂")
CO3 = FieldTimeSeries("buoyancy_front.jld2", "c₃")
OH= FieldTimeSeries("buoyancy_front.jld2", "c₅")

xc, yc, zc = nodes(CO2)

times = CO2.times

using CairoMakie

n = Observable(1)

CO2_lims = (0.0 , 7.5699999999999995e-6)
HCO3_lims = ( 0.0 , 0.00334)
CO3_lims = ( 0.0 , 0.00315)
OH_lims = ( 0.0 , 9.6e-6)

CO2ₙ = @lift interior(CO2[$n], :, 1, :)
HCO3ₙ = @lift interior(HCO3[$n], :, 1, :)
CO3ₙ = @lift interior(CO3[$n], :, 1, :)
OHₙ = @lift interior(OH[$n], :, 1, :)

fig = Figure(size = (1000, 1000), fontsize = 20)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title)

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 500, yticks = [-0.025,-0.01,0,0.01,0.025])
ax1 = Axis(fig[1, 1]; title = "CO2 mol/kg", axis_kwargs...)
ax2 = Axis(fig[2, 1]; title = "HCO3- mol/kg",axis_kwargs...)
ax3 = Axis(fig[3, 1]; title = "CO3-2 mol/kg", axis_kwargs...)
ax4 = Axis(fig[4, 1]; title = "OH- mol/kg", axis_kwargs...)

hm1 = heatmap!(ax1, xc, zc, CO2ₙ, colorrange = CO2_lims, colormap = Reverse(:lajolla), interpolate = true)
hm2 = heatmap!(ax2, xc, zc, HCO3ₙ, colorrange = HCO3_lims, colormap = Reverse(:bamako), interpolate = true)
hm3 = heatmap!(ax3, xc, zc, CO3ₙ, colorrange = CO3_lims, colormap = Reverse(:bamako), interpolate = true)
hm4 = heatmap!(ax4, xc, zc, OHₙ, colorrange = OH_lims, colormap = Reverse(:bamako), interpolate = true)

Colorbar(fig[1, 2], hm1, ticks = [0.0 , 7.5699999999999995e-6])
Colorbar(fig[2, 2], hm2, ticks = [0.0 , 0.00334])
Colorbar(fig[3, 2], hm3, ticks = [ 0.0 , 0.00315])
Colorbar(fig[4, 2], hm4, ticks = [ 0.0 ,  9.6e-6])

rowgap!(fig.layout, 0)

record(fig, "images/2D_results/organism_results.gif", 1:length(times)) do i
    n[] = i
end
CO2_min = minimum(interior(CO2[1], :, 1, :))
CO2_max = maximum(interior(CO2[1], :, 1, :))
println("CO2 concentration range: $CO2_min to $CO2_max")

HCO3_min = minimum(interior(HCO3[1], :, 1, :))
HCO3_max = maximum(interior(HCO3[1], :, 1, :))
println("HCO3 concentration range: $HCO3_min to $HCO3_max")

CO3_min = minimum(interior(CO3[1], :, 1, :))
CO3_max = maximum(interior(CO3[1], :, 1, :))
println("CO3 concentration range: $CO3_min to $CO3_max")

OH_min = minimum(interior(OH[1], :, 1, :))
OH_max = maximum(interior(OH[1], :, 1, :))
println("OH concentration range: $OH_min to $OH_max")