include("../DynamicCarbonateChemistry/evolutionequations.jl")
include("../DynamicCarbonateChemistry/QSS.jl")
include("../DynamicCarbonateChemistry/QSS2.jl")


using OceanBioME, Oceananigans
using Oceananigans.Units
using CairoMakie

function setup_model(c1, c2, c3, c4, c5, c6, c7, biogeochemistry, include_c4=true,include_c5=true)
    clock = Clock(; time = 0.0)
    @inline temp(t) = 2.4 * cos(t * 2π / year + 50days) + 26

    model = BoxModel(; 
        biogeochemistry,
        prescribed_tracers = (; T = temp),
        clock
    )
    set!(model, 
        c₁ = c1,    # CO₂
        c₂ = c2,    # HCO₃⁻
        c₃ = c3,    # CO₃²⁻
        #c₅ = c5,    # OH⁻
        c₆ = c6,    # B(OH)₃
        c₇ = c7     # B(OH)₄⁻
    )
    if include_c4
        set!(model, c₄ = c4)   # H⁺
    end
    if include_c5
        set!(model, c₅ = c5)   # OH⁻
    end
    return model
end

function run_simulation(model, filename,timestep,stop_time,time_int)
    simulation = Simulation(model; 
        Δt = timestep, 
        stop_time = stop_time
    )
    
    simulation.output_writers[:fields] = JLD2OutputWriter(
        model, 
        model.fields;
        filename = filename,
        schedule = TimeInterval(time_int),
        overwrite_existing = true
    )
    
    run!(simulation)
    return simulation
end

function plot_combined_results(c1_1, c2_1, c3_1, c4_1, c5_1, c6_1, c7_1, times_1,
                               c1_2, c2_2, c3_2, c4_2, c5_2, c6_2, c7_2, times_2, include_c4_1=true, include_c4_2=true, include_c5_1=true, include_c5_2=true)
    fig = Figure(size = (1200, 1200), fontsize = 20)
    
    # CO2 plot
    ax1 = Axis(fig[1, 1], ylabel = "CO2 \n(mol/kg)", xlabel = "Time (seconds)")
    ylims!(ax1, 0, 1e-4)
    lines!(ax1, times_1, c1_1[1, 1, 1, :], linewidth = 3, color = :blue, label = "Full model Simulation")
    lines!(ax1, times_2, c1_2[1, 1, 1, :], linewidth = 3, color = :red, linestyle = :dash, label = " Reduced Simulation 2")
    axislegend(ax1)

    # HCO3 plot
    ax2 = Axis(fig[1, 2], ylabel = "HCO3 \n(mol/kg)", xlabel = "Time (seconds)")
    ylims!(ax2, 0, 0.01)
    lines!(ax2, times_1, c2_1[1, 1, 1, :], linewidth = 3, color = :blue, label = "Full model Simulation ")
    lines!(ax2, times_2, c2_2[1, 1, 1, :], linewidth = 3, color = :red, linestyle = :dash, label = "Reduced Simulation 2")
    axislegend(ax2)

    # CO3 plot
    ax3 = Axis(fig[2, 1], ylabel = "CO₃²⁻ \n(mol/kg)", xlabel = "Time (seconds)")
    ylims!(ax3, 0, 0.0005)
    lines!(ax3, times_1, c3_1[1, 1, 1, :], linewidth = 3, color = :blue, label = "Full model Simulation")
    lines!(ax3, times_2, c3_2[1, 1, 1, :], linewidth = 3, color = :red, linestyle = :dash, label = "Reduced Simulation 2")
    axislegend(ax3)

    # H⁺ plot
    if include_c4_1 && include_c4_2
        ax4 = Axis(fig[2, 2], ylabel = "H⁺ \n(mol/kg)", xlabel = "Time (seconds)")
        ylims!(ax4, 0, 1e-7)
        lines!(ax4, times_1, c4_1[1, 1, 1, :], linewidth = 3, color = :blue, label = "Full model Simulation ")
        lines!(ax4, times_2, c4_2[1, 1, 1, :], linewidth = 3, color = :red, linestyle = :dash, label = "Reduced Simulation 2")
        axislegend(ax4)
    elseif include_c4_1
        ax4 = Axis(fig[2, 2], ylabel = "H⁺ \n(mol/kg)", xlabel = "Time (seconds)")
        ylims!(ax4, 0, 1e-7)
        lines!(ax4, times_1, c4_1[1, 1, 1, :], linewidth = 3, color = :blue, label = "Full model Simulation")
        axislegend(ax4)
    elseif include_c4_2
        ax4 = Axis(fig[2, 2], ylabel = "H⁺ \n(mol/kg)", xlabel = "Time (seconds)")
        ylims!(ax4, 0, 1e-7)
        lines!(ax4, times_2, c4_2[1, 1, 1, :], linewidth = 3, color = :red, linestyle = :dash, label = "Reduced Simulation 2")
        axislegend(ax4)
    end

    # OH⁻ plot
    if include_c5_1 && include_c5_2
        ax5 = Axis(fig[3, 1], ylabel = "OH⁻ \n(mol/kg)", xlabel = "Time (seconds)")
        ylims!(ax5, 0, 1e-4)
        lines!(ax5, times_1, c5_1[1, 1, 1, :], linewidth = 3, color = :blue, label = "Full model Simulation")
        lines!(ax5, times_2, c5_2[1, 1, 1, :], linewidth = 3, color = :red, linestyle = :dash, label = "Reduced Simulation 2")
        axislegend(ax5)
    elseif include_c5_1
        ax5 = Axis(fig[3, 1], ylabel = "OH⁻ \n(mol/kg)", xlabel = "Time (seconds)")
        ylims!(ax5, 0, 1e-4)
        lines!(ax5, times_1, c5_1[1, 1, 1, :], linewidth = 3, color = :blue, label = "Full model Simulation")
        axislegend(ax5)
    elseif include_c5_2
        ax5 = Axis(fig[3, 1], ylabel = "OH⁻ \n(mol/kg)", xlabel = "Time (seconds)")
        ylims!(ax5, 0, 1e-4)
        lines!(ax5, times_2, c5_2[1, 1, 1, :], linewidth = 3, color = :red, linestyle = :dash, label = "Reduced Simulation 2")
        axislegend(ax5)
    end

    # B(OH)₃ plot
    ax6 = Axis(fig[3, 2], ylabel = "B(OH)₃ \n(mol/kg)", xlabel = "Time (seconds)")
    ylims!(ax6, 0, 1e-3)
    lines!(ax6, times_1, c6_1[1, 1, 1, :], linewidth = 3, color = :blue, label = "Full model Simulation")
    lines!(ax6, times_2, c6_2[1, 1, 1, :], linewidth = 3, color = :red, linestyle = :dash, label = "Reduced Simulation 2")
    axislegend(ax6)

    # B(OH)₄⁻ plot
    ax7 = Axis(fig[4, 1], ylabel = "B(OH)₄⁻ \n(mol/kg)", xlabel = "Time (seconds)")
    ylims!(ax7, 0, 1e-3)
    lines!(ax7, times_1, c7_1[1, 1, 1, :], linewidth = 3, color = :blue, label = "Full model Simulation")
    lines!(ax7, times_2, c7_2[1, 1, 1, :], linewidth = 3, color = :red, linestyle = :dash, label = "Reduced Simulation 2")
    axislegend(ax7)

    ax8 = Axis(fig[4, 2], ylabel = "Total carbon \n(mol/kg)", xlabel = "Time (seconds)")
    ylims!(ax8, 0, 1e-2)
    lines!(ax8, times_1, c1_1[1, 1, 1, :]+c2_1[1, 1, 1, :]+c3_1[1, 1, 1, :], linewidth = 3, color = :blue, label = "Full model Simulation")
    lines!(ax8, times_2, c1_2[1, 1, 1, :]+c2_2[1, 1, 1, :]+c3_2[1, 1, 1, :], linewidth = 3, color = :red, linestyle = :dash, label = "Reduced Simulation 2")
    axislegend(ax8)

    return fig
end

function main(include_c4_1=false, include_c4_2=false, include_c5_1=false, include_c5_2=true)
    # First simulation
    model1 = setup_model(7.57*10^-6, 1.67*10^-3, 3.15*10^-4, 6.31*10^-9, 9.60*10^-6, 2.97e-4, 1.19e-4, DynamicCarbonateChemistry3(), include_c4_1, include_c5_1)
    simulation1 = run_simulation(model1, "box_np1.jld2", 0.004seconds, 60seconds, 0.5seconds)
    
    # Second simulation with different initial conditions
    model2 = setup_model(7.57*10^-6, 1.67*10^-3, 3.15*10^-4, 6.31*10^-9, 9.60*10^-6, 2.97e-4, 1.19e-4, DynamicCarbonateChemistry2(), include_c4_2, include_c5_2)
    simulation2 = run_simulation(model2, "box_np2.jld2", 0.004seconds, 60seconds, 0.5seconds)
    
    # Load results from both simulations
    c1_1 = FieldTimeSeries("box_np1.jld2", "c₁")
    c2_1 = FieldTimeSeries("box_np1.jld2", "c₂")
    c3_1 = FieldTimeSeries("box_np1.jld2", "c₃")
    c4_1 = include_c4_1 ? FieldTimeSeries("box_np1.jld2", "c₄") : nothing
    c5_1 = include_c5_1 ? FieldTimeSeries("box_np1.jld2", "c₅") : nothing
    c6_1 = FieldTimeSeries("box_np1.jld2", "c₆")
    c7_1 = FieldTimeSeries("box_np1.jld2", "c₇")
    
    c1_2 = FieldTimeSeries("box_np2.jld2", "c₁")
    c2_2 = FieldTimeSeries("box_np2.jld2", "c₂")
    c3_2 = FieldTimeSeries("box_np2.jld2", "c₃")
    c4_2 = include_c4_2 ? FieldTimeSeries("box_np2.jld2", "c₄") : nothing
    c5_2 = include_c5_2 ? FieldTimeSeries("box_np2.jld2", "c₅") : nothing
    c6_2 = FieldTimeSeries("box_np2.jld2", "c₆")
    c7_2 = FieldTimeSeries("box_np2.jld2", "c₇")
    
    # Plot combined results
    fig = plot_combined_results(c1_1, c2_1, c3_1, c4_1, c5_1, c6_1, c7_1, c1_1.times,
                                c1_2, c2_2, c3_2, c4_2, c5_2, c6_2, c7_2, c1_2.times, include_c4_1, include_c4_2, include_c5_1, include_c5_2)
    save("combined_model_results2.png", fig)
    display(fig)
end

main()