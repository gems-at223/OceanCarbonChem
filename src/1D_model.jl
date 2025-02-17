include("QSS.jl")
using Oceananigans, OceanBioME
using Oceananigans.Units
using CairoMakie




grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (32, ), x = 1, y = 1, z = (-100, 0))


biogeochemistry = Biogeochemistry(DynamicCarbonateChemistry())
                               


# put the model together

model = NonhydrostaticModel(; grid,
                              biogeochemistry,
                              closure = ScalarDiffusivity(κ=1e-1))

set!(model, 
    c₁ = 7.57*10^-6,    # CO₂ (scaled up from 1.5e-5)
    c₂ = 1.67*10^-2,    # HCO₃⁻ (scaled up from 1.9e-3)
    c₃ = 3.15*10^-4,    # CO₃²⁻ (scaled up from 2.5e-4)
    #c₄ = 6.31*10^-9,   # H⁺ (scaled up from 3.16e-8)
    c₅ = 9.60*10^-6,   # OH⁻ (scaled up from 3.16e-7)
    c₆ = 2.97e-4,   # B(OH)₃ (scaled up from 3.75e-4)
    c₇ = 1.19e-4    # B(OH)₄⁻ (scaled up from 1.25e-4)
    )
# run

simulation = Simulation(model, Δt = 0.0001seconds, stop_time = 1seconds)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "column_np.jld2",
                                                       schedule = TimeInterval(0.005seconds),
                                                       overwrite_existing = true)


run!(simulation)

N = FieldTimeSeries("column_np.jld2", "c₁")
P = FieldTimeSeries("column_np.jld2", "c₃")
#sed = FieldTimeSeries("column_np_sediment.jld2", "N_storage")

fig = Figure()

axN = Axis(fig[1, 1], ylabel = "z (m)")
axP = Axis(fig[2, 1], ylabel = "z (m)")
#axSed = Axis(fig[3, 1:2], ylabel = "Sediment (mmol N / m²)", xlabel = "Time (years)")

_, _, zc = nodes(grid, Center(), Center(), Center())
times = N.times

hmN = heatmap!(axN, times , zc, N[1, 1, 1:grid.Nz, 1:end]',
               interpolate = true, colormap = Reverse(:batlow))

hmP = heatmap!(axP, times, zc, P[1, 1, 1:grid.Nz, 1:end]',
               interpolate = true, colormap = Reverse(:batlow))

#lines!(axSed, times ./ year, sed[1, 1, 1, :])

Colorbar(fig[1, 2], hmN, label = "CO₂ (mol/ kg)")
Colorbar(fig[2, 2], hmP, label = "CO₃ (mol / kg)")

save("1D_model.png",fig)

display(fig)

sleep(10)