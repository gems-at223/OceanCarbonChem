include("evolutionequations.jl")
include("visualization.jl")

using OceanBioME, Oceananigans
using Oceananigans.Units



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
    set!(model, 
        c₁ = 1.5e-1,    # CO₂
        c₂ = 1.9e-3,    # HCO₃⁻
        c₃ = 2.5e-4,    # CO₃²⁻
        c₄ = 3.16e-5,   # H⁺
        c₅ = 3.16e-5,   # OH⁻
        c₆ = 3.75e-8,   # B(OH)₃
        c₇ = 1.25e-6    # B(OH)₄⁻
    )
    
    return model
end

function run_simulation(model)
    simulation = Simulation(model; 
        Δt = 0.0000007seconds, 
        stop_time = 0.04seconds
    )
    
    simulation.output_writers[:fields] = JLD2OutputWriter(
        model, 
        model.fields;
        filename = "box_np.jld2",
        schedule = TimeInterval(0.00001seconds),
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
    c4 = FieldTimeSeries("box_np.jld2", "c₄")
    c5 = FieldTimeSeries("box_np.jld2", "c₅")
    c6 = FieldTimeSeries("box_np.jld2", "c₆")
    c7 = FieldTimeSeries("box_np.jld2", "c₇")
    
    # Plot results
    fig = plot_results(c1, c2, c3, c4, c5, c6, c7, c1.times)
    display(fig)
end

main()