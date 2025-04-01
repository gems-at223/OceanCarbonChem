include("../DynamicCarbonateChemistry/QSS.jl")
include("../utils/param_init.jl")

using OceanBioME, Oceananigans
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.BoundaryConditions
using Oceananigans.Grids: xnodes, ynodes, znodes

grid = RectilinearGrid(CPU(), size = (128,128),
x = (-0.002meters, 0.002meters),y=1, z=(-0.002meters, 0.002meters), topology = (Bounded, Flat, Bounded))
println(grid)
biogeochemistry = DynamicCarbonateChemistry2()

params=calculate_params(8.1, 420)
c₁=params[1]
c₂=params[2]
c₃=params[3]
c₅=params[4]
c₆=params[5]
c₇=params[6]

R=0.0002



@inline function inside_cylinder(x, z)
    return ((x)^2 + (z)^2) <= R^2
end

@inline function distance_from_center(x, z)
    return sqrt(x^2 + z^2)
end


function compute_area_around_boundary(grid, R, distance)
    total_area = 0.0
    
    # Get the coordinates of the grid points in the x and z directions
    x_coords = xnodes(grid,Center())
    z_coords = znodes(grid,Center())
    
    # Define grid dimensions based on the coordinate arrays
    Nx = length(x_coords)
    Nz = length(z_coords)
    
    # Calculate cell size (dx and dz are usually available from the grid)
    dx = 0.004 / 128
    dz = 0.004 / 128
    
    # Loop through all grid points
    for i in 1:Nx
        for k in 1:Nz
            # Get point coordinates directly from the arrays
            x = x_coords[i]
            z = z_coords[k]
            
            # Calculate distance from center of cylinder
            dist = sqrt(x^2 + z^2)
            
            # Check if this point is within the specified distance from the boundary
            if !inside_cylinder(x, z) && (abs(dist - R) <= distance)
               # println("Point at x=$x, z=$z is within the specified distance from the boundary")
                # Add the cell's area to the total
                # For more accuracy, you might want to calculate the actual area
                # of the cell associated with this point
                total_area += dx * dz
            end
        end
    end
    
    return total_area
end

underlying_grid=grid
immersed_grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(inside_cylinder))
distance=0.00002
area=compute_area_around_boundary(immersed_grid,R,0.00002)
println("Area: $area")


@inline function CO2_forcing(x,z,t,c₁)
    dist = distance_from_center(x, z)
    if !inside_cylinder(x, z)  && (abs(dist - R) <= distance)
        return - (((12.7e-9/3600)/(area))*(c₁/((5*10^-6)+c₁)))
    else
        return 0.0  
    end
end

@inline function respiration(x,z,t,c₁)
    dist = distance_from_center(x, z)
    if !inside_cylinder(x, z)  && (abs(dist - R) <= distance)
        return ((3e-9)/3600)/(area)
    else
        return 0.0  
    end
end


@inline function calcification(x,z,t,c₃)
    dist = distance_from_center(x, z)
    if !inside_cylinder(x, z)  && (abs(dist - R) <= distance)

        return -((3.25e-9)/3600)/(area)
    else
        return 0.0  
    end
end

@inline function HCO3_forcing(x,z,t,c₁)
    dist = distance_from_center(x, z)
    if !inside_cylinder(x, z)  &&  (abs(dist - R) <= distance)
        return - (((12.7e-9/3600)/(area))*(1-(c₁/((5*10^-6)+c₁))))
    else
        return 0.0  
    end
end

@inline function OH_forcing(x,z,t,c₁)
    dist = distance_from_center(x, z)
    if !inside_cylinder(x, z) && (abs(dist - R) <= distance)
        return  (((12.7e-9/3600)/(area))*(1-(c₁/((5*10^-6)+c₁))))
    else
        return 0.0  
    end
end

boundary_conditions = (
    c₁ = FieldBoundaryConditions(top = ValueBoundaryCondition(10.5*10^-6), 
                                 bottom = ValueBoundaryCondition(10.5*10^-6),
                                 west = ValueBoundaryCondition(10.5*10^-6),
                                 east = ValueBoundaryCondition(10.5*10^-6)),
    
   c₂ = FieldBoundaryConditions(top = ValueBoundaryCondition(1.698*10^-3), 
                                 bottom = ValueBoundaryCondition(1.698*10^-3),
                                 west = ValueBoundaryCondition(1.698*10^-3), 
                                 east = ValueBoundaryCondition(1.698*10^-3)),
    
    c₃ = FieldBoundaryConditions(top = ValueBoundaryCondition(3.15*10^-4), 
                                 bottom = ValueBoundaryCondition(3.15*10^-4),
                                 west = ValueBoundaryCondition(3.15*10^-4),
                                 east = ValueBoundaryCondition(3.15*10^-4)),

    c₅ = FieldBoundaryConditions(top = ValueBoundaryCondition(10.1*10^-6), 
                                 bottom = ValueBoundaryCondition(10.1*10^-6),
                                 west = ValueBoundaryCondition(10.1*10^-6),
                                 east = ValueBoundaryCondition(10.1*10^-6)),
                                 
    c₆ = FieldBoundaryConditions(top = ValueBoundaryCondition(2.97e-4), 
                                 bottom = ValueBoundaryCondition(2.97e-4),
                                 west = ValueBoundaryCondition(2.97e-4),
                                 east = ValueBoundaryCondition(2.97e-4)),
                                 
    c₇ = FieldBoundaryConditions(top = ValueBoundaryCondition( 1.19e-4), 
                                 bottom = ValueBoundaryCondition( 1.19e-4),
                                 west = ValueBoundaryCondition( 1.19e-4),
                                 east = ValueBoundaryCondition( 1.19e-4))
)
 forcing = Forcing(CO2_forcing, field_dependencies=(:c₁))
 forcinghco3=Forcing(HCO3_forcing, field_dependencies=(:c₁))
 OHforcing=Forcing(OH_forcing, field_dependencies=(:c₁))
 resp=Forcing(respiration, field_dependencies=(:c₁))
calc=Forcing(calcification, field_dependencies=(:c₃))





model = NonhydrostaticModel(; grid=immersed_grid,
                            biogeochemistry = biogeochemistry,
                              #advection=WENO(;grid),
                              tracers = ( :c₁, :c₂, :c₃,:c₅,:c₆, :c₇),
                              boundary_conditions = boundary_conditions,
                              closure = ScalarDiffusivity(ν=1e-6,κ=(c₁ = 2e-9,
                               c₂ = 1.18e-9, c₃ = 0.955e-9, c₅ = 5.27e-9, c₆ = 1.11e-9, c₇ = 0.98e-9)),
                               #forcing=(;c₁=(forcing),c₂=(forcinghco3),c₅=(OHforcing)),)
                              #forcing = (;c₁=(resp)))
                              forcing=(;c₃=(calc)))
                              

#,c₂=(forcinghco3),c₅ = (OHforcing)
@inline front(x, z, μ, δ) = μ + δ * tanh((x - 0.01 + 0.00004 * z) / 0.005)

#Tᵢ(x, z) = front(x, z, 9, 0.05)

set!(model, 
    c₁ = 10.370*10^-6,   # CO₂ (scaled up from 1.5e-5)
    c₂ = 1.698*10^-3,    # HCO₃⁻ (scaled up from 1.9e-3)
    c₃ = 3.15*10^-4 ,    # CO₃²⁻ (scaled up from 2.5e-4)
    c₅ = 10.10*10^-6,   # OH⁻ (scaled up from 3.16e-7)
    c₆ = 2.97e-4,   # B(OH)₃ (scaled up from 3.75e-4)
    c₇ = 1.19e-4    # B(OH)₄⁻ (scaled up from 1.25e-4)
    )



#set!(model, c₁=c1ᵢ, c₂=c2i, c₃=c3ᵢ, c₆=c6ᵢ, c₇=c7ᵢ)

simulation = Simulation(model; Δt = 0.01seconds, stop_time = 300seconds)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "buoyancy_front.jld2",
                                                       schedule = TimeInterval(10seconds),
                                                       overwrite_existing = true)

run!(simulation)

CO2 = FieldTimeSeries("buoyancy_front.jld2", "c₁")
HCO3 = FieldTimeSeries("buoyancy_front.jld2", "c₂")
CO3 = FieldTimeSeries("buoyancy_front.jld2", "c₃")
OH= FieldTimeSeries("buoyancy_front.jld2", "c₅")
BOH3= FieldTimeSeries("buoyancy_front.jld2", "c₆")
BOH4= FieldTimeSeries("buoyancy_front.jld2", "c₇")

xc, yc, zc = nodes(CO2)

times = CO2.times

using CairoMakie

n = Observable(1)

CO2_lims = (10.37e-6 , 20.55e-6)
HCO3_lims = ( (1.6*10^-3) ,  1.7*10^-3)
CO3_lims = (1.7*10^-4 , 3.19*10^-4)
OH_lims = ( 8e-6
,  9.5e-6)

CO2ₙ = @lift interior(CO2[$n], :, 1, :)
HCO3ₙ = @lift interior(HCO3[$n], :, 1, :)
CO3ₙ = @lift interior(CO3[$n], :, 1, :)
OHₙ = @lift interior(OH[$n], :, 1, :)

fig = Figure(size = (1000, 1000), fontsize = 20)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title)

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 500, yticks = [-0.01,-0.005,0,0.005,0.01])
ax1 = Axis(fig[1, 1]; title = "CO2 mol/kg", axis_kwargs...)
ax2 = Axis(fig[2, 1]; title = "HCO3- mol/kg",axis_kwargs...)
ax3 = Axis(fig[3, 1]; title = "CO3-2 mol/kg", axis_kwargs...)
ax4 = Axis(fig[4, 1]; title = "OH- mol/kg", axis_kwargs...)

hm1 = heatmap!(ax1, xc, zc, CO2ₙ, colorrange = CO2_lims, colormap = Reverse(:lajolla), interpolate = true)
hm2 = heatmap!(ax2, xc, zc, HCO3ₙ, colorrange = HCO3_lims, colormap = Reverse(:bamako), interpolate = true)
hm3 = heatmap!(ax3, xc, zc, CO3ₙ, colorrange = CO3_lims, colormap = Reverse(:bamako), interpolate = true)
hm4 = heatmap!(ax4, xc, zc, OHₙ, colorrange = OH_lims, colormap = Reverse(:bamako), interpolate = true)

Colorbar(fig[1, 2], hm1, ticks = [10.37e-6 , 20.55e-6])
Colorbar(fig[2, 2], hm2, ticks = [1.7*10^-3 ,  1.8*10^-3])
Colorbar(fig[3, 2], hm3, ticks = [ 1.7*10^-4 , 3.19*10^-4])
Colorbar(fig[4, 2], hm4, ticks = [ 8e-6
 ,  10.5e-6])

rowgap!(fig.layout, 0)

record(fig, "images/2D_results/organism_calc_300_sec_128xz.gif", 1:length(times)) do i
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



tolerance = 0.01

z_center_index = argmin(abs.(zc))
if z_center_index === nothing
    error("Center z-plane not found in zc array.")
end

x_center_index = argmin(abs.(xc))
if x_center_index === nothing
    error("Center x-plane not found in xc array.")
end

println("Center x-plane index: $x_center_index")

println("Center z-plane index: $z_center_index")
println(zc)
# Define the radius of the immersed boundary

# Extract concentration data along the center z-plane, avoiding the immersed boundary
CO2_profile = [interior(CO2[end], i, 1, z_center_index)[] for i in 1:length(xc) if xc[i] > R]
HCO3_profile = [interior(HCO3[end], i, 1, z_center_index)[] for i in 1:length(xc) if xc[i] > R]
CO3_profile = [interior(CO3[end], i, 1, z_center_index)[] for i in 1:length(xc) if xc[i] > R]
OH_profile = [interior(OH[end], i, 1, z_center_index)[] for i in 1:length(xc) if xc[i] > R]
BOH3_profile=[interior(BOH3[end], i, 1, z_center_index)[] for i in 1:length(xc) if xc[i] > R]
BOH4_profile=[interior(BOH4[end], i, 1, z_center_index)[] for i in 1:length(xc) if xc[i] > R]

# Extract the corresponding x-coordinates, avoiding the immersed boundary
x_profile = [xc[i] for i in 1:length(xc) if xc[i] > R]

println(x_profile)
println(xc)

pOH=-log10.(OH_profile)
pH=14 .- pOH .-0.783

# Plot concentration vs. x-axis
fig_line = Figure()
ax_line = Axis(fig_line[1, 1], xlabel = "x (m)", ylabel = "CO2 (mol/kg)")

#ylims!(ax_line, 0, 60e-6)

HCO3_line=Axis(fig_line[1, 2], xlabel = "x (m)", ylabel = "HCO3⁻ (mol/kg)")
CO3_line=Axis(fig_line[2, 1], xlabel = "x (m)", ylabel = "CO3²⁻ (mol/kg)")
#BOH3_line=Axis(fig_line[2, 2], xlabel = "Length (m)", ylabel = "BOH3 (mol / kg)")
OH_line=Axis(fig_line[3, 1], xlabel = "x (m)", ylabel = "OH⁻ (mol/kg)")
BOH4_line=Axis(fig_line[3, 2], xlabel = "Length (m)", ylabel = "BOH4 (mol / kg)")
H_line=Axis(fig_line[2, 2], xlabel = "x (m)", ylabel = "pH")

#ylims!(ax_line, 0, 80)



lines!(ax_line, x_profile*10^6, CO2_profile*10^6, label = "CO2")
lines!(HCO3_line, x_profile*10^6, HCO3_profile*10^6, label = "HCO3")
lines!(CO3_line, x_profile*10^6, CO3_profile*10^6, label = "CO3")
lines!(OH_line, x_profile*10^6, OH_profile*10^6, label = "OH")
lines!(H_line, x_profile*10^6, pH, label = "BOH3")
lines!(BOH4_line, x_profile*10^6, BOH4_profile*10^6, label = "BOH4")

#axislegend(ax_line)

save("images/2D_results/concentration_vs_x_calc_300_sec_128xz.png", fig_line)
display(fig_line)

println(CO3_profile)



