include("../DynamicCarbonateChemistry/evolutionequations.jl")

using Oceananigans, OceanBioME
using Oceananigans.Units
using Oceananigans.BoundaryConditions

using CairoMakie




grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (16, ), x = 1, y = 1, z = (-0.01meters, 0meters))


biogeochemistry = Biogeochemistry(DynamicCarbonateChemistry())
                               
CO₂_sink_rate = 1e-4  # Example sink rate value

#function CO₂_sink(z,c₁)
    #if z > -0.1 && z < -0.05  # Apply sink only in the middle of the domain
        #return 0
    #else
        #return (((-0.05+0.1)*(z-0.05)^2))*(c₁/(5*10^6+c₁))
    #end
#end
#function HCO₃_sink(z,c₁)
    #if z > -0.1 && z < -0.05  # Apply sink only in the middle of the domain
     #   return 0
    #else
     #   return (((-0.05+0.1)*(z-0.05)^2))*(1-(c₁/(5*10^6+c₁)))
    #end
#end
#function HCO₃_sink(z,c₁)
    #if z > -0.1125 && z < -0.0875  # Apply sink only in the middle of the domain
    #    return 0
    #else
    #    return ((((-0.0875+0.1125)*(z+0.1)^2)))*(1-(c₁/(5*10^6+c₁)))
    #end
#end
function HCO₃_sink(z,c₁)
    return ((((-0.02)*(z)^2)))*(1-(c₁/(5*10^6+c₁)))
end

function CO₂_sink(z,c₁)
    return ((((-0.02)*(z)^2)))*((c₁/(5*10^6+c₁)))

end


Co2_forcing=Forcing(CO₂_sink, field_dependencies=(:c₁))
HCO₃_forcing=Forcing(HCO₃_sink, field_dependencies=(:c₁))
function custom_boundary_condition(z)
    if z > -0.125 && z < -0.075  # No flux between -0.1 and -0.05
        return 0.0
    else
        return c
    end
end

#FieldBoundaryConditions(top = FluxBoundaryCondition(CO₂_flux), bottom = GradientBoundaryCondition(0.0)),

CO₂_flux = 1e-5  # Example flux value
boundary_conditions = (c₃ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₂ = FieldBoundaryConditions(top = FluxBoundaryCondition(CO₂_flux), bottom = GradientBoundaryCondition(0.0)),
                       c₁ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₄ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₅ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₆ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₇ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)))
biogeochemistry = Biogeochemistry(DynamicCarbonateChemistry())

model = NonhydrostaticModel(; grid,
                              biogeochemistry,
                              closure = ScalarDiffusivity(κ=2e-4),
                              boundary_conditions = boundary_conditions,
                              tracers = (:c₁, :c₂, :c₃, :c₄, :c₅, :c₆, :c₇),
                              #forcing=(;c₁=(CO₂_sink),c₂=(HCO₃_sink))
                             )
        

set!(model, 
    c₁ = 7.57*10^-6,    # CO₂ (scaled up from 1.5e-5)
    c₂ = 1.67*10^-4,    # HCO₃⁻ (scaled up from 1.9e-3)
    c₃ = 3.15*10^-4,    # CO₃²⁻ (scaled up from 2.5e-4)
    c₄ = 6.31*10^-9,   # H⁺ (scaled up from 3.16e-8)
    c₅ = 9.60*10^-6,   # OH⁻ (scaled up from 3.16e-7)
    c₆ = 2.97e-4,   # B(OH)₃ (scaled up from 3.75e-4)
    c₇ = 1.19e-4    # B(OH)₄⁻ (scaled up from 1.25e-4)
    )
# run

simulation = Simulation(model, Δt = 0.0000005seconds, stop_time = 0.01seconds)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "column_np.jld2",
                                                       schedule = TimeInterval(0.0005seconds),
                                                       overwrite_existing = true)


run!(simulation)

CO2 = FieldTimeSeries("column_np.jld2", "c₁")
HCO3 = FieldTimeSeries("column_np.jld2", "c₂")
CO3= FieldTimeSeries("column_np.jld2", "c₃")
H= FieldTimeSeries("column_np.jld2", "c₄")
OH= FieldTimeSeries("column_np.jld2", "c₅")
BOH= FieldTimeSeries("column_np.jld2", "c₇")

#sed = FieldTimeSeries("column_np_sediment.jld2", "N_storage")

fig = Figure()

axCO2 = Axis(fig[1, 1], ylabel = "z (m)")
axHCO3 = Axis(fig[2, 1], ylabel = "z (m)")
axCO3 = Axis(fig[3, 1], ylabel = "z (m)")
axH = Axis(fig[4, 1], ylabel = "z (m)")
axOH = Axis(fig[5, 1], ylabel = "z (m)")
#axBOH = Axis(fig[6, 1], ylabel = "z (m)")

#axSed = Axis(fig[3, 1:2], ylabel = "Sediment (mmol N / m²)", xlabel = "Time (years)")

_, _, zc = nodes(grid, Center(), Center(), Center())
times = N.times

hmCO2 = heatmap!(axCO2, times , zc, CO2[1, 1, 1:grid.Nz, 1:end]',
               interpolate = true, colormap = Reverse(:batlow))

hmHCO3 = heatmap!(axHCO3, times, zc, HCO3[1, 1, 1:grid.Nz, 1:end]',
               interpolate = true, colormap = Reverse(:batlow))

hmCO3 = heatmap!(axCO3, times, zc, CO3[1, 1, 1:grid.Nz, 1:end]',
               interpolate = true, colormap = Reverse(:batlow))

hmH = heatmap!(axH, times, zc, H[1, 1, 1:grid.Nz, 1:end]',
                interpolate = true, colormap = Reverse(:batlow))

hmOH = heatmap!(axOH, times, zc, OH[1, 1, 1:grid.Nz, 1:end]',interpolate = true, colormap = Reverse(:batlow))

#hmBOH = heatmap!(axBOH, times, zc, BOH[1, 1, 1:grid.Nz, 1:end]',
 #                interpolate = true, colormap = Reverse(:batlow))




#lines!(axSed, times ./ year, sed[1, 1, 1, :])

Colorbar(fig[1, 2], hmCO2, label = "CO₂ (mol/ kg)")
Colorbar(fig[2, 2], hmHCO3, label = "HCO₃ (mol / kg)")
Colorbar(fig[3, 2], hmCO3, label = "HCO3- (mol / kg)")
Colorbar(fig[4, 2], hmH, label = "H⁺ (mol / kg)")
#Colorbar(fig[5, 2], hmB, label = "B(OH)₃ (mol / kg)")
#Colorbar(fig[6, 2], hmBOH, label = "B(OH)₄⁻ (mol / kg)")


save("images/1D_results/1D_model.png",fig)

display(fig)

fig_line = Figure()
ax_line = Axis(fig_line[1, 1], xlabel = "Length (m)", ylabel = "CO2 (mol / kg)")
HCO3_line=Axis(fig_line[1, 2], xlabel = "Length (m)", ylabel = "HCO3⁻ (mol / kg)")
CO3_line=Axis(fig_line[2, 1], xlabel = "Length (m)", ylabel = "CO3²⁻ (mol / kg)")
H_line=Axis(fig_line[2, 2], xlabel = "Length (m)", ylabel = "H⁺ (mol / kg)")
OH_line=Axis(fig_line[3, 1], xlabel = "Length (m)", ylabel = "OH⁻ (mol / kg)")


final_time_index = size(HCO3, 4)
CO2_data_final = CO2[:, :, :, final_time_index]
HCO3_data_final = HCO3[:, :, :, final_time_index]
CO3_data_final = CO3[:, :, :, final_time_index]
H_data_final = H[:, :, :, final_time_index]
OH_data_final = OH[:, :, :, final_time_index]


# Extract OH data across the depth of the domain
CO2_profile = CO2_data_final[1, 1, 1:16]
HCO3_profile = HCO3_data_final[1, 1, 1:16]
CO3_profile = CO3_data_final[1, 1, 1:16]
H_profile = H_data_final[1, 1, 1:16]
OH_profile = OH_data_final[1, 1, 1:16]



lines!(ax_line, zc,CO2_profile, label = "CO2 at final time step")
lines!(HCO3_line, zc,HCO3_profile, label = "HCO3 at final time step")
lines!(CO3_line, zc,CO3_profile, label = "CO3 at final time step")
lines!(H_line, zc,H_profile, label = "H at final time step")
lines!(OH_line, zc,OH_profile, label = "OH at final time step")


save("images/1D_results/line_plot.png", fig_line)

display(fig_line)

sleep(10)