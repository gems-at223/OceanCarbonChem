using OceanBioME, Oceananigans
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Units
using CairoMakie

const year = 365 * 24 * 60 * 60  # seconds


α1 = 0.037 # s^(-1)
β1 = 2.66e4 # kg/mol*s
α2 = 4.05e3 # kg/mol*s
β2 = 1.76e-4 # s^(-1)
α3 = 5.0e10 # kg/mol*s
β3 = 59.4 # s^(-1)
α4 = 6.0e9 # kg/mol*s
β4 = 3.06e5
α5 = 1.40e-3 # kg/mol*s
β5 = 2.31e-10 # kg/mol*s
α6 = 1.04e7 # kg/mol*s
β6 = 249.0 # s^(-1)
α7 = 6.92e6 # kg/mol*s
β7= 3.26e6 # kg/mol*s

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

@kwdef struct CarbonChemistry{FT} <: AbstractContinuousFormBiogeochemistry
            α1 :: FT = 0.037 # s^(-1)
            β1 :: FT = 2.66e4 # kg/mol*s
            α2 :: FT = 4.05e3 # kg/mol*s
            β2 :: FT = 1.76e-4 # s^(-1)
            α3 :: FT = 5.0e10 # kg/mol*s
            β3 :: FT = 59.4 # s^(-1)
            α4 :: FT = 6.0e9 # kg/mol*s
            β4 :: FT = 3.06e5
            α5 :: FT = 1.40e-3 # kg/mol*s
            β5 :: FT = 2.31e-10 # kg/mol*s
            α6 :: FT = 1.04e7 # kg/mol*s
            β6 :: FT = 249.0 # s^(-1)
            α7 :: FT = 6.92e6 # kg/mol*s
            β7 :: FT = 3.26e6 # kg/mol*s
end

required_biogeochemical_tracers(::CarbonChemistry) = (:c1, :c2, :c3, :c6, :c7)



function c4_solution(c1, c2, c3, c6, c7, α1, α2, α3, α4, α5, α6, β1, β2, β3, β4, β5, β6)
    A = β5 * (β1 * c2 + α3 * c3)
    B = (β1 * c2 + α3 * c3) * (α2 * c1 + α4 * c2 + α6 * c6) + β5 * (β2 * c2 + β4 * c3 + β6 * c7 - α1 * c1 - β3 * c2)
    C = -(α1 * c1 + β3 * c2 + α5) * (α2 * c1 + α4 * c2 + α6 * c6)
    
    discriminant = B^2 - 4 * A * C

    if discriminant == 0
        c4 = -B / (2 * A)
        println(c4)
        return c4
    elseif discriminant > 0
        root1 = (-B + sqrt(discriminant)) / (2 * A)
        root2 = (-B - sqrt(discriminant)) / (2 * A)

        c4 = max(root1, root2, 0)
        println(c4)
        return c4
    else
        error("c4 is imaginary, cannot compute c5")
    end
end

function c5_solution(c1, c2, c3, c4, c6, c7, α2, α4, α5, α6, β2, β4, β5, β6)

    return (β2 * c2 + β4 * c3 + α5 + β6 * c7) / (α2 * c1 + α4 * c2 + β5 * c4 + α6 * c6)
end

@inline function (bgc::CarbonChemistry)(::Val{:c1}, x, y, z, t, c1, c2, c3, c6, c7)
    c4 = c4_solution(c1, c2, c3, c6, c7, bgc.α1, bgc.α2, bgc.α3, bgc.α4, bgc.α5, bgc.α6, bgc.β1, bgc.β2, bgc.β3, bgc.β4, bgc.β5, bgc.β6)
    if c4 === nothing
        error("c4 is nothing, cannot compute dc1")
    end
    c5 = c5_solution(c1, c2, c3, c4, c6, c7, bgc.α2, bgc.α4, bgc.α5, bgc.α6, bgc.β2, bgc.β4, bgc.β5, bgc.β6)
    dc1 = -(bgc.α1 + bgc.α2 * c5) * c1 + (bgc.β1 * c4 + bgc.β2) * c2
    return dc1
end

@inline function (bgc::CarbonChemistry)(::Val{:c2}, x, y, z, t, c1, c2, c3, c6, c7)
    c4 = c4_solution(c1, c2, c3, c6, c7, bgc.α1, bgc.α2, bgc.α3, bgc.α4, bgc.α5, bgc.α6, bgc.β1, bgc.β2, bgc.β3, bgc.β4, bgc.β5, bgc.β6)
    if c4 === nothing
        error("c4 is nothing, cannot compute dc2")
    end
    c5 = c5_solution(c1, c2, c3, c4, c6, c7, bgc.α2, bgc.α4, bgc.α5, bgc.α6, bgc.β2, bgc.β4, bgc.β5, bgc.β6)
    dc2 = (bgc.α1 + bgc.α2 * c5) * c1 - (bgc.β1 * c4 + bgc.β2 + bgc.β3 + bgc.α4 * c5 + bgc.β7 * c7) * c2 + (bgc.α3 * c4 + bgc.β4 + bgc.α7 * c6) * c3
    return dc2
end

@inline function (bgc::CarbonChemistry)(::Val{:c3}, x, y, z, t, c1, c2, c3, c6, c7)
    c4 = c4_solution(c1, c2, c3, c6, c7, bgc.α1, bgc.α2, bgc.α3, bgc.α4, bgc.α5, bgc.α6, bgc.β1, bgc.β2, bgc.β3, bgc.β4, bgc.β5, bgc.β6)
    if c4 === nothing
        error("c4 is nothing, cannot compute dc3")
    end
    c5 = c5_solution(c1, c2, c3, c4, c6, c7, bgc.α2, bgc.α4, bgc.α5, bgc.α6, bgc.β2, bgc.β4, bgc.β5, bgc.β6)
    dc3 = (bgc.β3 + bgc.α4 * c5 + bgc.β7 * c7) * c2 - (bgc.α3 * c4 + bgc.β4 + bgc.α7 * c6) * c3
    return dc3
end

@inline function (bgc::CarbonChemistry)(::Val{:c6}, x, y, z, t, c1, c2, c3, c6, c7)
    c4 = c4_solution(c1, c2, c3, c6, c7, bgc.α1, bgc.α2, bgc.α3, bgc.α4, bgc.α5, bgc.α6, bgc.β1, bgc.β2, bgc.β3, bgc.β4, bgc.β5, bgc.β6)
    if c4 === nothing
        error("c4 is nothing, cannot compute dc6")
    end
    c5 = c5_solution(c1, c2, c3, c4, c6, c7, bgc.α2, bgc.α4, bgc.α5, bgc.α6, bgc.β2, bgc.β4, bgc.β5, bgc.β6)
    dc6 = bgc.β7 * c7 * c2 - bgc.α7 * c6 * c3 - (bgc.α6 * c5 * c6 - bgc.β6 * c7)
    return dc6
end

@inline function (bgc::CarbonChemistry)(::Val{:c7}, x, y, z, t, c1, c2, c3, c6, c7)
    c4 = c4_solution(c1, c2, c3, c6, c7, bgc.α1, bgc.α2, bgc.α3, bgc.α4, bgc.α5, bgc.α6, bgc.β1, bgc.β2, bgc.β3, bgc.β4, bgc.β5, bgc.β6)
    if c4 === nothing
        error("c4 is nothing, cannot compute dc7")
    end
    c5 = c5_solution(c1, c2, c3, c4, c6, c7, bgc.α2, bgc.α4, bgc.α5, bgc.α6, bgc.β2, bgc.β4, bgc.β5, bgc.β6)
    dc7 = - bgc.β7 * c7 * c2 + bgc.α7 * c6 * c3 + (bgc.α6 * c5 * c6 - bgc.β6 * c7)
    return dc7
end


using Oceananigans.Fields: FunctionField

clock = Clock(; time = 0.0)
@inline temp(t) = 2.4 * cos(t * 2π / year + 50days) + 26

biogeochemistry = Biogeochemistry(CarbonChemistry())

model = BoxModel(; biogeochemistry, prescribed_tracers = (; T = temp), clock)

set!(model, c1 = 7.57e-6, c2 = 1.67e-3, c3 = 3.15e-4, c6 = 2.97e-4, c7 = 1.19e-4) # mol/kg

simulation = Simulation(model; Δt = 1e-3seconds, stop_time = 1seconds)

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "box_np.jld2", schedule = TimeInterval(1e-3seconds), overwrite_existing = true)

@info "Running the model..."
@info "Before run: Time = $(model.clock.time)"
run!(simulation)
@info "After run: Time = $(model.clock.time)"

c1 = FieldTimeSeries("box_np.jld2", "c1")
c2 = FieldTimeSeries("box_np.jld2", "c2")
c3 = FieldTimeSeries("box_np.jld2", "c3")
c6 = FieldTimeSeries("box_np.jld2", "c6")
c7 = FieldTimeSeries("box_np.jld2", "c7")

times = c1.times

fig = Figure(size = (1500, 700), fontsize = 15, tight = true)

ax_c1 = Axis(fig[1, 1], ylabel = "CO₂ Conc. [mol/kg]", xlabel = "Time [s]")
lines!(ax_c1, times, c1[:], linewidth = 3)


ax_c2 = Axis(fig[1, 2], ylabel = "HCO₃⁻ Conc. [mol/kg]", xlabel = "Time [s]")
limits!(ax_c2, 0, 1e-3, 0, 2e-3)
lines!(ax_c2, times, c2[:], linewidth = 3)

ax_c3 = Axis(fig[2, 1], ylabel = "CO₃²⁻ Conc. [mol/kg]", xlabel = "Time [s]")
limits!(ax_c3, 0, 1e-3, 0, 1e-3)
lines!(ax_c3, times, c3[:], linewidth = 3)

ax_c4 = Axis(fig[2, 2], ylabel = "H⁺ Conc. [mol/kg]", xlabel = "Time [s]")
limits!(ax_c4, 0, 1e-3, 0, 1e-8)
lines!(ax_c4, times, c4[:], linewidth = 3)

ax_c5 = Axis(fig[3, 1], ylabel = "OH⁻ Conc. [mol/kg]", xlabel = "Time [s]")
limits!(ax_c5, 0, 1e-3, 0, 3e-5)
lines!(ax_c5, times, c5[:], linewidth = 3)

ax_c6 = Axis(fig[3, 2], ylabel = "B(OH)₃ Conc. [mol/kg]", xlabel = "Time [s]")
limits!(ax_c6, 0, 1e-3, 0, 5e-4)
lines!(ax_c6, times, c6[:], linewidth = 3)

ax_c7 = Axis(fig[4, 1], ylabel = "B(OH)₄⁻ Conc. [mol/kg]", xlabel = "Time [s]")
limits!(ax_c7, 0, 1e-3, 0, 2e-4)
lines!(ax_c7, times, c7[:], linewidth = 3)

fig
save("carbon_chemistry_plot_3.png", fig, dpi = 300)