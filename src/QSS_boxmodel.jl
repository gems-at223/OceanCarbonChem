#include("evolutionequations.jl")
include("visualization.jl")
include("QSS.jl")

using OceanBioME, Oceananigans
using Oceananigans.Units
using CairoMakie




function setup_model()
    biogeochemistry = Biogeochemistry(DynamicCarbonateChemistry())
    clock = Clock(; time = 0.0)
    @inline temp(t) = 2.4 * cos(t * 2π / year + 50days) + 26

    
    model = BoxModel(; 
        biogeochemistry,
        prescribed_tracers = (; T = temp),
        clock
    )
    
    # Initial conditions
    #set!(model, 
        #c₁ = 1.5e-1,    # CO₂
        #c₂ = 1.9e-3,    # HCO₃⁻
        #c₃ = 2.5e-4,    # CO₃²⁻
        #c₄ = 3.16e-5,   # H⁺
        #c₅ = 3.16e-5,   # OH⁻
        #c₆ = 3.75e-8,   # B(OH)₃
        #c₇ = 1.25e-6    # B(OH)₄⁻
#)
    set!(model, 
    c₁ = 7.57*10^-6,    # CO₂ (scaled up from 1.5e-5)
    c₂ = 1.67*10^-2,    # HCO₃⁻ (scaled up from 1.9e-3)
    c₃ = 3.15*10^-4,    # CO₃²⁻ (scaled up from 2.5e-4)
    #c₄ = 6.31*10^-9,   # H⁺ (scaled up from 3.16e-8)
    c₅ = 9.60*10^-6,   # OH⁻ (scaled up from 3.16e-7)
    c₆ = 2.97e-4,   # B(OH)₃ (scaled up from 3.75e-4)
    c₇ = 1.19e-4    # B(OH)₄⁻ (scaled up from 1.25e-4)
    )

    
    return model
end

function run_simulation(model)
    simulation = Simulation(model; 
        Δt = 0.0000009seconds, 
        stop_time = 1seconds
    )
    
    simulation.output_writers[:fields] = JLD2OutputWriter(
        model, 
        model.fields;
        filename = "box_np.jld2",
        schedule = TimeInterval(0.0001seconds),
        overwrite_existing = true
    )
    
    run!(simulation)
    return simulation
end

function main()
    model = setup_model()
    simulation = run_simulation(model)
    
    # Load results
    c1 = FieldTimeSeries("box_np.jld2", "c₁")
    c2 = FieldTimeSeries("box_np.jld2", "c₂")
    c3 = FieldTimeSeries("box_np.jld2", "c₃")
    #c4 = FieldTimeSeries("box_np.jld2", "c₄")
    c5 = FieldTimeSeries("box_np.jld2", "c₅")
    c6 = FieldTimeSeries("box_np.jld2", "c₆")
    c7 = FieldTimeSeries("box_np.jld2", "c₇")
    
    # Plot results
    fig = plot_results(c1, c2, c3, c3, c5, c6, c7, c1.times)
    save("perturbed_model_results.png",fig)
    display(fig)
end

main()