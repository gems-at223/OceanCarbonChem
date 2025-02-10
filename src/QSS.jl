using OceanBioME, Oceananigans
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Units

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity

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
   println(co2_evolution)

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

