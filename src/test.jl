using OceanBioME, Oceananigans
using Oceananigans.Units
using CairoMakie
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Units
include("visualization.jl")
import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity

@kwdef struct CarbonChemistry{FT} <: AbstractContinuousFormBiogeochemistry
     α1 :: FT = 0.037 # s^(-1)
     β1 :: FT = 2.66*10^4 # kg/mol*s
     α2 :: FT = 4.05*10^3 # kg/mol*s
     β2 :: FT = 1.76*10^-4 # s^(-1)
     α3 :: FT = 5.0*10^10 # kg/mol*s
     β3 :: FT = 59.4 # s^(-1)
     α4 :: FT = 6.0*10^9 # kg/mol*s
     β4 :: FT = 3.06*10^5
     α5 :: FT = 1.40*10^-3 # kg/mol*s
     β5 :: FT = 2.31*10^-10 # kg/mol*s
     α6 :: FT = 1.04*10^7 # kg/mol*s
     β6 :: FT = 249.0 # s^(-1)
     α7 :: FT = 6.92*10^6 # kg/mol*s
     β7 :: FT = 3.26*10^6 # kg/mol*s

end

required_biogeochemical_tracers(::CarbonChemistry) = (:c1, :c2, :c3, :c5, :c6, :c7)

@inline function (bgc::CarbonChemistry)(::Val{:c1}, x, y, z,t, c1, c2, c3, c5, c6, c7)
    c4 = (bgc.α1 * c1 + bgc.β3 * c2 + bgc.α5) / (bgc.β1 * c2 + bgc.α3 * c3 + bgc.β5 * c5)
    dc1 = -(bgc.α1 + bgc.α2 * c5) * c1 + (bgc.β1 * c4 + bgc.β2) * c2
    return dc1
end

@inline function (bgc::CarbonChemistry)(::Val{:c2}, x, y, z,t, c1, c2, c3, c5, c6, c7)
    c4 = (bgc.α1 * c1 + bgc.β3 * c2 + bgc.α5) / (bgc.β1 * c2 + bgc.α3 * c3 + bgc.β5 * c5)
    dc2 = (bgc.α1 + bgc.α2 * c5) * c1 - (bgc.β1 * c4 + bgc.β2 + bgc.β3 + bgc.α4 * c5 + bgc.β7 * c7) * c2 + (bgc.α3 * c4 + bgc.β4 + bgc.α7 * c6) * c3
    println(dc2)
    return dc2
end

@inline function (bgc::CarbonChemistry)(::Val{:c3}, x, y, z,t, c1, c2, c3, c5, c6, c7)
    c4 = (bgc.α1 * c1 + bgc.β3 * c2 + bgc.α5) / (bgc.β1 * c2 + bgc.α3 * c3 + bgc.β5 * c5)
    dc3 = (bgc.β3 + bgc.α4 * c5 + bgc.β7 * c7) * c2 - (bgc.α3 * c4 + bgc.β4 + bgc.α7 * c6) * c3
    return dc3
end

@inline function (bgc::CarbonChemistry)(::Val{:c5}, x, y, z,t, c1, c2, c3, c5, c6, c7)
    c4 = (bgc.α1 * c1 + bgc.β3 * c2 + bgc.α5) / (bgc.β1 * c2 + bgc.α3 * c3 + bgc.β5 * c5)
    dc5 = -1*(bgc.α2 * c5 * c1) + (bgc.β2 - bgc.α4 * c5) * c2 + (bgc.β4 * c3) + (bgc.α5 - bgc.β5 * c4 * c5) - (bgc.α6 * c5 * c6 - bgc.β6 * c7)
    return dc5
end

@inline function (bgc::CarbonChemistry)(::Val{:c6}, x, y, z,t, c1, c2, c3, c5, c6, c7)
    dc6 = (bgc.β7 * c7 * c2) - (bgc.α7 * c6 * c3) - (bgc.α6 * c5 * c6 - bgc.β6 * c7)
return dc6
end

@inline function (bgc::CarbonChemistry)(::Val{:c7}, x, y, z,t, c1, c2, c3, c5, c6, c7)
    dc7 = -1*(bgc.β7 * c7 * c2) + (bgc.α7 * c6 * c3) + (bgc.α6 * c5 * c6 - bgc.β6 * c7)
    return dc7
end


function setup_model()
    biogeochemistry = Biogeochemistry(CarbonChemistry())
    clock = Clock(; time = 0.0)
    @inline temp(t) = 2.4 * cos(t * 2π / year + 50days) + 26

    
    model = BoxModel(; 
        biogeochemistry,
        prescribed_tracers = (; T = temp),
        clock
    )
    set!(model, 
    c1 = 7.57*10^-6,    # CO₂ (scaled up from 1.5e-5)
    c2 = 5*1.67*10^-3,    # HCO₃⁻ (scaled up from 1.9e-3)
    c3 = 3.15*10^-4,    # CO₃²⁻ (scaled up from 2.5e-4)
    #c4 = 6.31*10^-9,   # H⁺ (scaled up from 3.16e-8)
    c5 = 9.60*10^-6,   # OH⁻ (scaled up from 3.16e-7)
    c6 = 2.97e-4,   # B(OH)₃ (scaled up from 3.75e-4)
    c7 = 1.19e-4    # B(OH)₄⁻ (scaled up from 1.25e-4)
    )
    return model
end

function run_simulation(model)
    simulation = Simulation(model; 
        Δt = 0.00000005seconds, 
        stop_time = 0.003seconds
    )
    
    simulation.output_writers[:fields] = JLD2OutputWriter(
        model, 
        model.fields;
        filename = "box_np.jld2",
        schedule = TimeInterval(0.0005seconds),
        overwrite_existing = true
    )
    
    run!(simulation)
    return simulation
end

function main()
    model = setup_model()
    simulation = run_simulation(model)
    
    # Load results
    c1 = FieldTimeSeries("box_np.jld2", "c1")
    c2 = FieldTimeSeries("box_np.jld2", "c2")
    c3 = FieldTimeSeries("box_np.jld2", "c3")
    c4 = FieldTimeSeries("box_np.jld2", "c3")
    c5 = FieldTimeSeries("box_np.jld2", "c5")
    c6 = FieldTimeSeries("box_np.jld2", "c6")
    c7 = FieldTimeSeries("box_np.jld2", "c7")
    
    # Plot results
    fig = plot_results(c1, c2, c3, c4, c5, c6, c7, c1.times)
    save("perturbed_full_model_results2.png",fig)
    display(fig)
end

main()