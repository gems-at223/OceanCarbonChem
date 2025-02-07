using OceanBioME, Oceananigans
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Units

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity

#@kwdef struct DynamicCarbonateChemistry{FT} <: AbstractContinuousFormBiogeochemistry
    #α₁ :: FT = 0.0037        ## s-1
    #β₁ :: FT = 2.66*10^4      # kg mol-1 s-1
    #α₄ :: FT = 4.05e3       # kg mol-1 s-1
    #β₄ :: FT = 1.76*10^-4     # s-1
    #α₅ :: FT = 10^10.0       # kg mol-1 s-1
    #β₅ :: FT = 9.0           # s-1
    #α₆ :: FT = 1.3*10^-3     # kg mol-1 s-1
    #β₆ :: FT = 2.8*10^10.0   # kg mol-1 s-1
    #α₇ :: FT = 20.0          # s-1
    #β₇ :: FT = 10^10.0       # kg mol-1 s-1
#end

@kwdef struct DynamicCarbonateChemistry{FT} <: AbstractContinuousFormBiogeochemistry
    α₁ :: FT = 0.0036        # s-1
    β₁ :: FT = 3.3*10^4      # kg mol-1 s-1
    α₄ :: FT = 8500.0        # kg mol-1 s-1
    β₄ :: FT = 3.7*10^-4     # s-1
    α₅ :: FT = 10^10.0       # kg mol-1 s-1
    β₅ :: FT = 9.0           # s-1
    α₆ :: FT = 1.3*10^-3     # kg mol-1 s-1
    β₆ :: FT = 2.8*10^10.0   # kg mol-1 s-1
    α₇ :: FT = 20.0          # s-1
    β₇ :: FT = 10^10.0       # kg mol-1 s-1
end


required_biogeochemical_tracers(::DynamicCarbonateChemistry) = (:c₁, :c₂, :c₃,:c₄,:c₅,:c₆,:c₇)


@inline function (bgc::DynamicCarbonateChemistry)(::Val{:c₁}, x, y, z, t, c₁, c₂, c₃,c₄,c₅,c₆,c₇)
   co2_evolution=(bgc.β₁*c₄+bgc.β₄)*c₂ - (bgc.α₁+bgc.α₄*c₅)*c₁
  # println(co2_evolution)

   return co2_evolution
end

@inline function (bgc::DynamicCarbonateChemistry)(::Val{:c₂}, x, y, z, t, c₁, c₂, c₃,c₄,c₅,c₆,c₇)
   HCO3_evolution=(bgc.α₁*c₁) - (bgc.β₁*c₄*c₂) + (bgc.α₄*c₁*c₅) - (bgc.β₄*c₂) + (bgc.α₅*c₄*c₃) - (bgc.β₅*c₂)

   return HCO3_evolution
end

@inline function (bgc::DynamicCarbonateChemistry)(::Val{:c₃}, x, y, z, t, c₁, c₂, c₃,c₄,c₅,c₆,c₇)
   CO3_evolution=(bgc.β₅*c₂) - (bgc.α₅*c₄*c₃)
   #println(CO3_evolution)

   return CO3_evolution
end

@inline function (bgc::DynamicCarbonateChemistry)(::Val{:c₄}, x, y, z, t, c₁, c₂, c₃,c₄,c₅,c₆,c₇)
   H_evolution=(bgc.β₅-bgc.β₁*c₄)*c₂ - (bgc.α₅*c₄*c₃) + (bgc.α₁*c₁) - (bgc.α₁*c₄*c₃) + (bgc.α₆) - (bgc.β₆*c₄*c₅) + (bgc.α₇*c₆) - (bgc.β₇*c₇*c₄)
   println(H_evolution)
   return H_evolution
end

@inline function (bgc::DynamicCarbonateChemistry)(::Val{:c₅}, x, y, z, t, c₁, c₂, c₃,c₄,c₅,c₆,c₇)
   OH_evolution=(bgc.β₄*c₂) - (bgc.α₄*c₁*c₅) + bgc.α₆ - bgc.β₆*c₄*c₅
   return OH_evolution
end

@inline function (bgc::DynamicCarbonateChemistry)(::Val{:c₆}, x, y, z, t, c₁, c₂, c₃,c₄,c₅,c₆,c₇)
   BOH3_evolution= - (bgc.α₇*c₆) + (bgc.β₇*c₇*c₄)
   
   return BOH3_evolution
end

@inline function (bgc::DynamicCarbonateChemistry)(::Val{:c₇}, x, y, z, t, c₁, c₂, c₃,c₄,c₅,c₆,c₇)
   B_evolution=(bgc.α₇*c₆) - (bgc.β₇*c₇*c₄)
   return B_evolution
end


biogeochemistry = Biogeochemistry(DynamicCarbonateChemistry())

clock = Clock(; time = 0.0)

@inline temp(t) = 2.4 * cos(t * 2π / year + 50days) + 26
const year = years = 365days


model = BoxModel(; biogeochemistry,
                   prescribed_tracers = (; T = temp),
                   clock)
               
#set!(model, c₁ = 10 ,c₂ = 10,c₃ = 10, c₄ = 10,c₅ = 10,c₆ = 10,c₇ = 10)
set!(model, 
    c₁ = 1.5e-1,    # CO₂ (scaled up from 1.5e-5)
    c₂ = 1.9e-3,    # HCO₃⁻ (scaled up from 1.9e-3)
    c₃ = 2.5e-4,    # CO₃²⁻ (scaled up from 2.5e-4)
   c₄ = 3.16e-5,   # H⁺ (scaled up from 3.16e-8)
    c₅ = 3.16e-5,   # OH⁻ (scaled up from 3.16e-7)
    c₆ = 3.75e-8,   # B(OH)₃ (scaled up from 3.75e-4)
    c₇ = 1.25e-6    # B(OH)₄⁻ (scaled up from 1.25e-4)
)

#set!(model, 
#    c₁ = 7.57*10^-6,    # CO₂ (scaled up from 1.5e-5)
 #   c₂ = 1.67*10^-3,    # HCO₃⁻ (scaled up from 1.9e-3)
 #   c₃ = 3.15*10^-4,    # CO₃²⁻ (scaled up from 2.5e-4)
  #  c₄ = 6.31*10^-9,   # H⁺ (scaled up from 3.16e-8)
   # c₅ = 9.60*10^-6,   # OH⁻ (scaled up from 3.16e-7)
    #c₆ = 2.97e-4,   # B(OH)₃ (scaled up from 3.75e-4)
   # c₇ = 1.19e-4    # B(OH)₄⁻ (scaled up from 1.25e-4)
#)

simulation = Simulation(model; Δt = 0.0000001seconds, stop_time = 0.01seconds)

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "box_np.jld2", schedule = TimeInterval(0.00001seconds), overwrite_existing = true)

using CairoMakie
#println(model.fields)


# ## Run the model (should only take a few seconds)
@info "Running the model..."
run!(simulation)

c1 = FieldTimeSeries("box_np.jld2", "c₁")
c2= FieldTimeSeries("box_np.jld2", "c₂")
c3 = FieldTimeSeries("box_np.jld2", "c₃")
c4 = FieldTimeSeries("box_np.jld2", "c₄")
c5 = FieldTimeSeries("box_np.jld2", "c₅")
c6 = FieldTimeSeries("box_np.jld2", "c₆")
c7 = FieldTimeSeries("box_np.jld2", "c₇")




times = c1.times

# ## And plot

fig = Figure(size = (1200, 1200), fontsize = 20)

#axN= Axis(fig[1, 1], ylabel = "Nutrient \n(mmol N / m³)")
#lines!(axN, times / year, P[1, 1, 1, :], linewidth = 3)

#axP = Axis(fig[1, 2], ylabel = "Phytoplankton \n(mmol N / m³)")
#lines!(axP, times / year, N[1, 1, 1, :], linewidth = 3)
axN = Axis(fig[1, 1], 
    ylabel = "CO2 \n(mmol/m³)",
    xlabel = "Time (days)")
   # limits = (0, 1e-8)
ylims!(axN, 0, 0.1)
lines!(axN, times, c1[1, 1, 1, :], linewidth = 3)

axP = Axis(fig[1, 2], 
    ylabel = "HCO3 \n(mmol/m³)",
    xlabel = "Time (days)")
   # limits = (0, 1e-8)
ylims!(axP, 0, 0.01)
lines!(axP, times, c2[1, 1, 1, :], linewidth = 3)

axZ = Axis(fig[2, 1], 
    ylabel = "CO₃²⁻ \n(mmol/m³)",
    xlabel = "Time (days)")
   # limits = (0, 1e-8)
ylims!(axZ, 0, 0.0005)
lines!(axZ, times, c3[1, 1, 1, :], linewidth = 3)

axD = Axis(fig[2, 2], 
    ylabel = "H⁺ \n(mmol/m³)",
    xlabel = "Time (days)")
   # limits = (0, 1e-8)
ylims!(axD, 0, 1e-8)
lines!(axD, times, c4[1, 1, 1, :], linewidth = 3)

axD = Axis(fig[3, 1], 
    ylabel = "OH⁻ \n(mmol/m³)",
    xlabel = "Time (days)")
   # limits = (0, 1e-8)
ylims!(axD, 0, 0.0001)
lines!(axD, times, c5[1, 1, 1, :], linewidth = 3)

axe = Axis(fig[3, 2], 
    ylabel = "B(OH)₃ \n(mmol/m³)",
    xlabel = "Time (days)")
   # limits = (0, 1e-8)
ylims!(axe, 0, 0.00001)
lines!(axe, times, c6[1, 1, 1, :], linewidth = 3)

axf = Axis(fig[4, 1], 
    ylabel = "B(OH)₄⁻ \n(mmol/m³)",
    xlabel = "Time (days)")
   # limits = (0, 1e-8)
ylims!(axf, 0, 0.00001)
lines!(axf, times, c7[1, 1, 1, :], linewidth = 3)


display(fig)

#axPAR= Axis(fig[2, 1], ylabel = "PAR (einstein / m² / s)", xlabel = "Time (years)")
#lines!(axPAR, times / year, PAR_func.(times), linewidth = 3)

#axT = Axis(fig[2, 2], ylabel = "Temperature (°C)", xlabel = "Time (years)")
#lines!(axT, times / minutes, temp.(times), linewidth = 3)

#display(fig)

sleep(10)