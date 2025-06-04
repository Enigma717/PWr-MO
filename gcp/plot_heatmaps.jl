using Plots

function plot_heatmaps()
    cx_probabilites::Vector{Float64} = [0.3, 0.4, 0.5, 0.6, 0.7]
    mx_probabilites::Vector{Float64} = [0.1, 0.15, 0.2, 0.25, 0.3]

    subgroup5::Matrix{Float64} = [
        8.0 8.0 8.0 7.8 8.0;
        8.0 8.0 8.0 8.0 8.0;
        8.0 8.0 8.0 8.0 8.0;
        8.0 8.0 8.0 7.8 8.0;
        8.0 7.8 8.0 8.0 8.0;
    ]

    subgroup6::Matrix{Float64} = [
        8.8 9.0 8.8 9.0 9.4;
        9.0 9.0 8.8 8.8 8.8;
        9.2 9.2 8.8 9.2 9.4;
        9.4 9.2 9.8 9.0 9.8;
        9.4 9.2 9.0 9.2 9.4;
    ]

    subgroup7::Matrix{Float64} = [
        8.8 8.8 8.8 8.6 8.8;
        8.6 8.8 8.0 9.0 9.0;
        9.0 9.0 8.8 9.0 9.0;
        9.2 8.4 8.2 9.0 9.2;
        9.6 9.2 8.6 9.2 9.6;
    ]

    subgroup8::Matrix{Float64} = [
        8.8 8.4 8.8 8.8 8.4;
        9.2 8.8 8.2 8.4 8.6;
        9.2 9.0 8.6 8.2 8.4;
        9.4 8.8 8.4 8.4 8.4;
        9.4 8.8 8.2 8.6 8.6;
    ]

    ga_plots_data::Vector{Matrix{Float64}} = [subgroup5, subgroup6, subgroup7, subgroup8]

    for (index, subgroup) in enumerate(ga_plots_data)
        heatmap_plot = heatmap(
            mx_probabilites,
            cx_probabilites,
            subgroup,
            xlabel = "Prawdopodobieństwo mutacji",
            ylabel = "Prawdopodobieństwo krzyżowania",
            yflip=true,
            clims=(7.5, 10),
            color = cgrad(:matter, rev = true),
            plot_title = "Wielkość grupy turniejowej $(4+index)")

        savefig(heatmap_plot, "./plots/tuning/subgroup$(4+index)_heatmap.pdf")
    end

    temperatures::Vector{Float64} = [1000000.0, 500000.0, 100000.0, 50000.0, 10000.0]
    cooling_rates::Vector{Float64} = [0.9999, 0.9995, 0.999, 0.995, 0.99]

    sa_data = [
        16.9 8.0 8.1 8.0 7.9;
        16.5 7.9 8.0 8.0 8.0;
        16.2 8.2 8.0 8.0 8.0;
        15.8 8.0 8.0 7.7 8.0;
        13.9 8.2 8.0 7.9 8.0;
    ]

    heatmap_plot = heatmap(
        sa_data,
        xlabel = "Temperatura początkowa",
        ylabel = "Współczynnik schładzania",
        yflip=true,
        clims=(7.5, 17.5),
        color = cgrad(:matter, rev = true),
        plot_title = "Symulowane wyżarzanie")

    savefig(heatmap_plot, "./plots/tuning/sa_heatmap.pdf")
end

plot_heatmaps()
