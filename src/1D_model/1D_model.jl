include("../DynamicCarbonateChemistry/QSS.jl")

using Oceananigans, OceanBioME
using Oceananigans.Units
using Oceananigans.BoundaryConditions

using CairoMakie




grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (64, ), x = 1, y = 1, z = (-0.01meters, 0meters))


                               
CO₂_sink_rate = 1e-9/3600 # Example sink rate value

function CO₂_sink(z,t,c₁)
    if z<-0.00
        return -(CO₂_sink_rate/(4*pi*(z+0.00)^2))*(c₁/((5*10^-6)+c₁))
    else
        return 0
    end
end

function HCO₃_sink(z,t,c₁)
    if z<-0.00
        return -(CO₂_sink_rate/(4*pi*(z+0.00)^2))*(1-(c₁/((5*10^-6)+c₁)))
    else
        return 0
    end
end

# Define the total uptake of HCO₃⁻ by integrating the HCO₃⁻ sink over the entire domain
function total_HCO₃_uptake( z,t,c₁)
    total_uptake = 0.0
    z_coords = znodes(grid,Center())
    Nz = length(z_coords)

    for zp in 1:Nz
        if zp<0.00
            total_uptake += HCO₃_sink(zp, t, c₁)* (0.01/64)
            total_uptake += 0
        end
    end
    return total_uptake
end

# Define the OH⁻ flux equivalent to the total uptake of HCO₃⁻ at z=0
function OH_forcing_func(z,t, c₁)
    if z>-0.0002
        return abs(total_HCO₃_uptake(z,t, c₁))
    else
        return 0
    end
end
function OH_sink(z,t,c₁)
    if z<0.00
        return (CO₂_sink_rate/(4*pi*(z+0.00)^2))*(1-(c₁/((5*10^-6)+c₁)))
    else
        return 0
    end
end

function CO2_source(z,t,c₁)
    if z==0
        return 3e-9/3600
    else
        return 0
    end
end


Co2_forcing=Forcing(CO₂_sink, field_dependencies=(:c₁))
HCO₃_forcing=Forcing(HCO₃_sink, field_dependencies=(:c₁))
OH_forcing=Forcing(OH_forcing_func, field_dependencies=(:c₁))

#FieldBoundaryConditions(top = FluxBoundaryCondition(CO₂_flux), bottom = GradientBoundaryCondition(0.0)),

CO₂_flux = -12e-9/3600  # Example flux value
CO₃_flux = 3e-9/3600
boundary_conditions = (c₃ = FieldBoundaryConditions(top = FluxBoundaryCondition(CO₃_flux), 
bottom = GradientBoundaryCondition(0.0)),
                       c₂ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₁ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₄ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₅ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₆ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)),
                       c₇ = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0), bottom = GradientBoundaryCondition(0.0)))
biogeochemistry = Biogeochemistry(DynamicCarbonateChemistry2())
model = NonhydrostaticModel(; grid,
                              biogeochemistry,
                              closure = ScalarDiffusivity(ν=1e-3,κ=(c₁ = 2e-9, c₂ = 1.18e-9, c₃ = 0.955e-9, c₅ = 5.27e-9, c₆ = 1.11e-9, c₇ = 0.98e-9)),
                              
                              boundary_conditions = boundary_conditions,
                              tracers = (:c₁, :c₂, :c₃, :c₅, :c₆, :c₇),
                              forcing=(;c₁=Co2_forcing,c₂=HCO₃_forcing,c₅=OH_forcing),
                              )
        

#set!(model, 
 #   c₁ = 7.57*10^-6,    # CO₂ (scaled up from 1.5e-5)
  #  c₂ = 1.67*10^-3,    # HCO₃⁻ (scaled up from 1.9e-3)
   # c₃ = 3.15*10^-4,    # CO₃²⁻ (scaled up from 2.5e-4)
   # c₄ = 6.31*10^-9,   # H⁺ (scaled up from 3.16e-8)
   # c₅ = 9.60*10^-6,   # OH⁻ (scaled up from 3.16e-7)
   # c₆ = 2.97e-4,   # B(OH)₃ (scaled up from 3.75e-4)
   # c₇ = 1.19e-4    # B(OH)₄⁻ (scaled up from 1.25e-4)
   # )

set!(model, 
    c₁ = 24.061856*10^-6,    # CO₂ (scaled up from 1.5e-5)
    c₂ = 3588.183336*10^-6,    # HCO₃⁻ (scaled up from 1.9e-3)
    c₃ = 422.754808*10^-6,    # CO₃²⁻ (scaled up from 2.5e-4)
    #c₄ = 8.31*10^-9,   # H⁺ (scaled up from 3.16e-8)
    c₅ = 1.16*10^−6,   # OH⁻ (scaled up from 3.16e-7)
    c₆ = 1626.474134*10^-6,   # B(OH)₃ (scaled up from 3.75e-4)
    c₇ = 428.611580*10^-6    # B(OH)₄⁻ (scaled up from 1.25e-4)
    )

#simulation = Simulation(model, Δt = 0.001seconds, stop_time = 20seconds)
simulation = Simulation(model, Δt = 0.002seconds, stop_time = 10seconds)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "column_np.jld2",
                                                       schedule = TimeInterval(2seconds),
                                                       overwrite_existing = true)


run!(simulation)

CO2 = FieldTimeSeries("column_np.jld2", "c₁")
HCO3 = FieldTimeSeries("column_np.jld2", "c₂")
CO3= FieldTimeSeries("column_np.jld2", "c₃")
#H= FieldTimeSeries("column_np.jld2", "c₄")
OH= FieldTimeSeries("column_np.jld2", "c₅")
BOH3= FieldTimeSeries("column_np.jld2", "c₆")
BOH4= FieldTimeSeries("column_np.jld2", "c₇")

#sed = FieldTimeSeries("column_np_sediment.jld2", "N_storage")

fig = Figure()

axCO2 = Axis(fig[1, 1], ylabel = "z (m)")
axHCO3 = Axis(fig[2, 1], ylabel = "z (m)")
axCO3 = Axis(fig[3, 1], ylabel = "z (m)")
#axH = Axis(fig[4, 1], ylabel = "z (m)")
axOH = Axis(fig[5, 1], ylabel = "z (m)")
#axBOH3 = Axis(fig[6, 1], ylabel = "z (m)")

#axSed = Axis(fig[3, 1:2], ylabel = "Sediment (mmol N / m²)", xlabel = "Time (years)")

_, _, zc = nodes(grid, Center(), Center(), Center())
times = CO2.times

hmCO2 = heatmap!(axCO2, times , zc, CO2[1, 1, 1:grid.Nz, 1:end]',
               interpolate = true, colormap = Reverse(:batlow))

hmHCO3 = heatmap!(axHCO3, times, zc, HCO3[1, 1, 1:grid.Nz, 1:end]',
               interpolate = true, colormap = Reverse(:batlow))

hmCO3 = heatmap!(axCO3, times, zc, CO3[1, 1, 1:grid.Nz, 1:end]',
               interpolate = true, colormap = Reverse(:batlow))

#hmH = heatmap!(axH, times, zc, H[1, 1, 1:grid.Nz, 1:end]',
 #               interpolate = true, colormap = Reverse(:batlow))

hmOH = heatmap!(axOH, times, zc, OH[1, 1, 1:grid.Nz, 1:end]',interpolate = true, colormap = Reverse(:batlow))

#hmBOH = heatmap!(axBOH3, times, zc, BOH[1, 1, 1:grid.Nz, 1:end]',
 #                interpolate = true, colormap = Reverse(:batlow))




#lines!(axSed, times ./ year, sed[1, 1, 1, :])

Colorbar(fig[1, 2], hmCO2, label = "CO₂ (mol/ kg)")
Colorbar(fig[2, 2], hmHCO3, label = "HCO₃ (mol / kg)")
Colorbar(fig[3, 2], hmCO3, label = "CO3- (mol / kg)")
#Colorbar(fig[4, 2], hmH, label = "H⁺ (mol / kg)")
Colorbar(fig[5, 2], hmOH, label = "OH⁻ (mol / kg)")
#Colorbar(fig[5, 2], hmB, label = "B(OH)₃ (mol / kg)")
#Colorbar(fig[6, 2], hmBOH, label = "B(OH)₄⁻ (mol / kg)")


save("1D_model_f.png",fig)

display(fig)

fig_line = Figure()
ax_line = Axis(fig_line[1, 1], xlabel = "r(m)", ylabel = "CO2 (mol / kg)")
HCO3_line=Axis(fig_line[1, 2], xlabel = "r (m)", ylabel = "HCO3⁻ (mol / kg)")
CO3_line=Axis(fig_line[2, 1], xlabel = "r (m)", ylabel = "CO3²⁻ (mol / kg)")
#H_line=Axis(fig_line[2, 2], xlabel = "Length (m)", ylabel = "H⁺ (mol / kg)")
OH_line=Axis(fig_line[3, 1], xlabel = "r (m)", ylabel = "OH⁻ (mol / kg)")
BOH3_line=Axis(fig_line[3, 2], xlabel = "r (m)", ylabel = "B(OH)₃ (mol / kg)")
BOH4_line=Axis(fig_line[2, 2], xlabel = "r (m)", ylabel = "B(OH)₄⁻ (mol / kg)")


final_time_index = size(HCO3, 4)
CO2_data_final = CO2[:, :, :, final_time_index]
HCO3_data_final = HCO3[:, :, :, final_time_index]
CO3_data_final = CO3[:, :, :, final_time_index]
#H_data_final = H[:, :, :, final_time_index]
OH_data_final = OH[:, :, :, final_time_index]
BOH3_data_final = BOH3[:, :, :, final_time_index]
BOH4_data_final = BOH4[:, :, :, final_time_index]


# Extract OH data across the depth of the domain
CO2_profile = CO2_data_final[1, 1, 1:64]*10^6
HCO3_profile = HCO3_data_final[1, 1, 1:64]*10^6
CO3_profile = CO3_data_final[1, 1, 1:64]*10^6
#H_profile = H_data_final[1, 1, 1:16]
OH_profile = OH_data_final[1, 1, 1:64]*10^6
BOH3_profile = BOH3_data_final[1, 1, 1:64]*10^6
BOH4_profile = BOH4_data_final[1, 1, 1:64]*10^6


#xticks = range(minimum(zc) * -1 * 10^6, stop = maximum(zc) * -1 * 10^6, length = 10)

lines!(ax_line, zc*-1*10^6,CO2_profile, label = "CO2 at final time step")
lines!(HCO3_line, zc*-1*10^6,HCO3_profile, label = "HCO3 at final time step")
lines!(CO3_line, zc*-1*10^6,CO3_profile, label = "CO3 at final time step")
#lines!(H_line, zc*-1,H_profile, label = "H at final time step")
lines!(OH_line, zc*-1,OH_profile, label = "OH at final time step")
lines!(BOH3_line, zc*-1*10^6,BOH3_profile, label = "BOH3 at final time step")
lines!(BOH4_line, zc*-1*10^6,BOH4_profile, label = "BOH4 at final time step")

#ax_line.xticks = (xticks, string.(xticks))


save("line_plot_f.png", fig_line)

display(fig_line)

sleep(10)