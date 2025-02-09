# src/visualization.jl
using CairoMakie

function plot_results(c1, c2, c3, c4, c5, c6, c7, times)
    fig = Figure(size = (1200, 1200), fontsize = 20)
    
    # CO2 plot
    ax1 = Axis(fig[1, 1], 
        ylabel = "CO2 \n(mmol/m³)",
        xlabel = "Time (days)")
    ylims!(ax1, 0, 0.1)
    lines!(ax1, times, c1[1, 1, 1, :], linewidth = 3)
    
    # HCO3 plot
    ax2 = Axis(fig[1, 2], 
        ylabel = "HCO3 \n(mmol/m³)",
        xlabel = "Time (days)")
    ylims!(ax2, 0, 0.01)
    lines!(ax2, times, c2[1, 1, 1, :], linewidth = 3)

    # CO3 plot
    ax3 = Axis(fig[2, 1], 
        ylabel = "CO₃²⁻ \n(mmol/m³)",
        xlabel = "Time (days)")
    ylims!(ax3, 0, 0.0005)
    lines!(ax3, times, c3[1, 1, 1, :], linewidth = 3)


    axD = Axis(fig[2, 2], 
    ylabel = "H⁺ \n(mmol/m³)",
    xlabel = "Time (seconds)")

    ylims!(axD, 0, 1e-8)
    lines!(axD, times, c4[1, 1, 1, :], linewidth = 3)

    axD = Axis(fig[3, 1], 
    ylabel = "OH⁻ \n(mmol/m³)",
    xlabel = "Time (seconds)")
   # limits = (0, 1e-8)
    ylims!(axD, 0, 0.0001)
    lines!(axD, times, c5[1, 1, 1, :], linewidth = 3)

    axe = Axis(fig[3, 2], 
    ylabel = "B(OH)₃ \n(mmol/m³)",
    xlabel = "Time (seconds)")
   # limits = (0, 1e-8)
    ylims!(axe, 0, 0.00001)
    lines!(axe, times, c6[1, 1, 1, :], linewidth = 3)

    axf = Axis(fig[4, 1], 
    ylabel = "B(OH)₄⁻ \n(mmol/m³)",
    xlabel = "Time (seconds)")
   # limits = (0, 1e-8)
    ylims!(axf, 0, 0.00001)
    lines!(axf, times, c7[1, 1, 1, :], linewidth = 3)


    display(fig)
end