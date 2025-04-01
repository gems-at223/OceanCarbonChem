using OceanBioME, Oceananigans
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Units

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity

@kwdef struct DynamicCarbonateChemistry3{FT} <: AbstractContinuousFormBiogeochemistry
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


required_biogeochemical_tracers(::DynamicCarbonateChemistry3) = (:c₁, :c₂, :c₃,:c₆,:c₇)


@inline function (bgc::DynamicCarbonateChemistry3)(::Val{:c₁}, x, y, z, t, c₁, c₂, c₃,c₆,c₇)
    a = bgc.β₁*c₂*bgc.β₆ + bgc.α₅*c₃*bgc.β₆ + bgc.β₇*c₇*bgc.β₆
    
    b = bgc.β₁*c₂*bgc.α₄*c₁ + bgc.α₅*c₃*bgc.α₄*c₁ + bgc.β₆*bgc.β₄*c₂ + 
        bgc.β₆*bgc.α₆ + bgc.β₇*c₇*bgc.α₄*c₁ - bgc.β₅*c₂*bgc.β₆ - 
        bgc.α₁*c₁*bgc.β₆ - bgc.α₆*bgc.β₆ - bgc.α₇*c₆*bgc.β₆
    
   c = -(bgc.β₅*c₂*bgc.α₄*c₁ + bgc.α₁*c₁*bgc.α₄*c₁ + bgc.α₆*bgc.α₄*c₁ + bgc.α₇*c₆*bgc.α₄*c₁) 
   discriminant = b^2 - 4*a*c
   root1 = (-b + sqrt(discriminant)) / (2*a)
   root2 = (-b - sqrt(discriminant)) / (2*a)
   if root1 >= 0
    c₄ = root1
   elseif root2 >= 0
    c₄ = root2
    elseif discriminant==0
        c₄=-b/(2*a)
        if c₄<0
            error("No positive root found for c₄ in discriminant in CO2")
        end
    else
        println(root1)
        println(root2)
        error("No positive root found for c₄ in CO2")
        println(c₄)
   end  
   c₅=(bgc.β₄*c₂+bgc.α₆)/(bgc.α₄*c₁+bgc.β₆*c₄)
   co2_evolution=(bgc.β₁*c₄+bgc.β₄)*c₂ - (bgc.α₁+bgc.α₄*c₅)*c₁
   println(co2_evolution)

   return co2_evolution
end

@inline function (bgc::DynamicCarbonateChemistry3)(::Val{:c₂}, x, y, z, t, c₁, c₂, c₃,c₆,c₇)

    a = bgc.β₁*c₂*bgc.β₆ + bgc.α₅*c₃*bgc.β₆ + bgc.β₇*c₇*bgc.β₆
    
    b = bgc.β₁*c₂*bgc.α₄*c₁ + bgc.α₅*c₃*bgc.α₄*c₁ + bgc.β₆*bgc.β₄*c₂ + 
        + bgc.β₇*c₇*bgc.α₄*c₁ - bgc.β₅*c₂*bgc.β₆ - 
        bgc.α₁*c₁*bgc.β₆  - bgc.α₇*c₆*bgc.β₆
    
   c = -(bgc.β₅*c₂*bgc.α₄*c₁ + bgc.α₁*c₁*bgc.α₄*c₁ + bgc.α₆*bgc.α₄*c₁ + bgc.α₇*c₆*bgc.α₄*c₁) 
   discriminant = b^2 - 4*a*c
   root1 = (-b + sqrt(discriminant)) / (2*a)
   root2 = (-b - sqrt(discriminant)) / (2*a)
   if root1 > 0
    c₄ = root1
   elseif root2 > 0
    c₄ = root2
   elseif discriminant==0
    c₄=-b/(2*a)
    else 
        println(root1)
        println(root2)
        error("No positive root found for c₄ HCO3")
        println(c₄)
   end  
   c₅=(bgc.β₄*c₂+bgc.α₆)/(bgc.α₄*c₁+bgc.β₆*c₄)
   #c₄= (bgc.β₅*c₂ + bgc.α₁*c₁ + bgc.α₆ + bgc.α₇*c₆)/(bgc.β₁*c₂ + bgc.α₅*c₃+bgc.β₆*c₅+bgc.β₇*c₇)
   HCO3_evolution=(bgc.α₁*c₁) - (bgc.β₁*c₄*c₂) + (bgc.α₄*c₁*c₅) - (bgc.β₄*c₂) + (bgc.α₅*c₄*c₃) - (bgc.β₅*c₂)

   return HCO3_evolution
end

@inline function (bgc::DynamicCarbonateChemistry3)(::Val{:c₃}, x, y, z, t, c₁, c₂, c₃,c₆,c₇)
   a = bgc.β₁*c₂*bgc.β₆ + bgc.α₅*c₃*bgc.β₆ + bgc.β₇*c₇*bgc.β₆
    
   b = bgc.β₁*c₂*bgc.α₄*c₁ + bgc.α₅*c₃*bgc.α₄*c₁ + bgc.β₆*bgc.β₄*c₂ + 
        bgc.β₆*bgc.α₆ + bgc.β₇*c₇*bgc.α₄*c₁ - bgc.β₅*c₂*bgc.β₆ - 
        bgc.α₁*c₁*bgc.β₆ - bgc.α₆*bgc.β₆ - bgc.α₇*c₆*bgc.β₆
    
   c = -(bgc.β₅*c₂*bgc.α₄*c₁ + bgc.α₁*c₁*bgc.α₄*c₁ + bgc.α₆*bgc.α₄*c₁ + bgc.α₇*c₆*bgc.α₄*c₁) 
   discriminant = b^2 - 4*a*c
   root1 = (-b + sqrt(discriminant)) / (2*a)
   root2 = (-b - sqrt(discriminant)) / (2*a)
   if root1 > 0
    c₄ = root1
   elseif root2 > 0
    c₄ = root2
   elseif discriminant==0
    c₄=-b/(2*a)
    else
        println(root1)
        println(root2)
        error("No positive root found for c₄ in CO3")
        
   end  
   #c₄= (bgc.β₅*c₂ + bgc.α₁*c₁ + bgc.α₆ + bgc.α₇*c₆)/(bgc.β₁*c₂ + bgc.α₅*c₃+bgc.β₆*c₅+bgc.β₇*c₇)
   CO3_evolution=(bgc.β₅*c₂) - (bgc.α₅*c₄*c₃)
   #println(CO3_evolution)

   return CO3_evolution
end


@inline function (bgc::DynamicCarbonateChemistry3)(::Val{:c₆}, x, y, z, t, c₁, c₂, c₃,c₆,c₇)
    a = bgc.β₁*c₂*bgc.β₆ + bgc.α₅*c₃*bgc.β₆ + bgc.β₇*c₇*bgc.β₆
    
    b = bgc.β₁*c₂*bgc.α₄*c₁ + bgc.α₅*c₃*bgc.α₄*c₁ + bgc.β₆*bgc.β₄*c₂ + 
        bgc.β₆*bgc.α₆ + bgc.β₇*c₇*bgc.α₄*c₁ - bgc.β₅*c₂*bgc.β₆ - 
        bgc.α₁*c₁*bgc.β₆ - bgc.α₆*bgc.β₆ - bgc.α₇*c₆*bgc.β₆
    
    c = -(bgc.β₅*c₂*bgc.α₄*c₁ + bgc.α₁*c₁*bgc.α₄*c₁ + bgc.α₆*bgc.α₄*c₁ + bgc.α₇*c₆*bgc.α₄*c₁) 
    discriminant = b^2 - 4*a*c
    root1 = (-b + sqrt(discriminant)) / (2*a)
    root2 = (-b - sqrt(discriminant)) / (2*a)
    if root1 > 0
     c₄ = root1
    elseif root2 > 0
     c₄ = root2
    elseif discriminant==0
        c₄=-b/(2*a)
    else
         println(root1)
         println(root2)
         error("No positive root found for c₄ in BOH3")
         println(c₄)
    end  
   BOH3_evolution= - (bgc.α₇*c₆) + (bgc.β₇*c₇*c₄)
   
   return BOH3_evolution
end

@inline function (bgc::DynamicCarbonateChemistry3)(::Val{:c₇}, x, y, z, t, c₁, c₂, c₃,c₆,c₇)
    a = bgc.β₁*c₂*bgc.β₆ + bgc.α₅*c₃*bgc.β₆ + bgc.β₇*c₇*bgc.β₆
    
    b = bgc.β₁*c₂*bgc.α₄*c₁ + bgc.α₅*c₃*bgc.α₄*c₁ + bgc.β₆*bgc.β₄*c₂ + 
        bgc.β₆*bgc.α₆ + bgc.β₇*c₇*bgc.α₄*c₁ - bgc.β₅*c₂*bgc.β₆ - 
        bgc.α₁*c₁*bgc.β₆ - bgc.α₆*bgc.β₆ - bgc.α₇*c₆*bgc.β₆
    
    c = -(bgc.β₅*c₂*bgc.α₄*c₁ + bgc.α₁*c₁*bgc.α₄*c₁ + bgc.α₆*bgc.α₄*c₁ + bgc.α₇*c₆*bgc.α₄*c₁) 
    discriminant = b^2 - 4*a*c
    root1 = (-b + sqrt(discriminant)) / (2*a)
    root2 = (-b - sqrt(discriminant)) / (2*a)
    if root1 > 0
     c₄ = root1
    elseif root2 > 0
     c₄ = root2
    elseif discriminant==0
        c₄=-b/(2*a)
     else
         println(root1)
         println(root2)
         error("No positive root found for c₄ in BOH4")
         println(c₄)
    end  
   B_evolution=(bgc.α₇*c₆) - (bgc.β₇*c₇*c₄)
   return B_evolution
end