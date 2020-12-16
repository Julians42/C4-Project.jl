using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, JLD2, Statistics, Glob, ColorSchemes, Plots.PlotUtils, HDF5

using Pkg 
ENV["GR"] = ""
Pkg.build("GR")

# ColorScheme
σs = 0.5:0.1:2
normal_x = -5:0.01:5
normal_y = [exp.(-normal_x.^2 / (2σ^2)) / (2π * σ^2) for σ in σs];
loadcolorscheme(:cm_maxamp,ColorSchemes.gist_heat.colors[end-50:-1:1], "maxamp color", "for waveform plot");

# Extract information and corr_data for a single station
function create_corr(file::HDF5File, station::String, stacktype::String, component::String)
    S = CorrData()
    # get corr data
    corr = file[station][component][stacktype][1:end]
    S.corr = reshape(corr,size(corr,1),1)
    S.name = station
    S.dist = read(CHN[names(CHN)[1]]["meta"]["dist"])
    S.loc.lat = read(CHN[names(CHN)[1]]["meta"]["lat"])
    S.loc.lon = read(CHN[names(CHN)[1]]["meta"]["lon"])
    S.loc.el = read(CHN[names(CHN)[1]]["meta"]["el"])
    S.azi = read(CHN[names(CHN)[1]]["meta"]["azi"])
    S.baz = read(CHN[names(CHN)[1]]["meta"]["baz"])
    S.fs = 20. # For the 20 Hz product
    return S
end

# get information for a set of filtered stations, default to all stations
function get_corrs(file::HDF5File, stacktype::String="linear",
    component::String="ZZ", filter::String="")
    filtered = [name for name in names(file) if occursin(filter, name)] 
    corrs = [create_corr(file, station, stacktype, component) for station in filtered]
    return corrs
end 
# Plotting
function vert_plot(corrs_passed::Array{CorrData,1}, sta::String="", comp_iter::String="ZZ", line::String="", stack_type::String="linear")
    # Add filtered correlations to appropriate plots
    T = collect(-300.:1/corrs_passed[1].fs:300.)
    y_max = length(corrs_passed)*5+9
    for j in 1:length(frequency_plots)
        #Select frequency
        fmin = frequency_plots[j][1]
        fmax = frequency_plots[j][2]
        # Define plots 
        ZZ_plot = plot(xlims = (0,100), ylims = (0,y_max), 
                        yticks=[],xlabel = "Time (s)", ylabel = "Stations (Near to Far)", 
                        xtickfontsize=5,ytickfontsize=5,fontsize=5, xguidefontsize = 10, yguidefontsize = 10,
                        legendfontsize = 15)
        # Adjust scaling by comparison by first corr
        Cstack = bandpass(corrs_passed[1], fmin, fmax)
        id_mid = round(Int, size(Cstack.corr, 1))
        scaled = scale*maximum(broadcast(abs, Cstack.corr))
        max_amplitudes = Array{Float64, 1}(undef,0)
        plot_process = []
        # Stack over days and add to plot
        for i in 1:length(corrs_passed)
            Cstack = bandpass(corrs_passed[i], fmin, fmax)
            push!(max_amplitudes, maximum(Cstack.corr/scaled))
            push!(plot_process, Cstack.corr/scaled .+5*i)

        end
        title!("$(sta) S$line $(comp_iter) $fmin to $fmax ($stack_type)", fontsize=5)
        plot!(ZZ_plot, -T, plot_process, color=:cm_maxamp, colorbar_title="Normalized Maximum Amplitude", 
            line_z=max_amplitudes', fmt = :png, linewidth = lw, reuse = false, legend = false)
        plot!(size=(250,400),dpi=1000)
        filepath = joinpath(rootdir,"stack_plots/$stack_type/S$(line)_$(sta)_$(comp_iter)_$(fmin)to$(fmax).png")
        # ensure filepath is valid 
        DIR = dirname(filepath)
        if !isdir(DIR)
            mkpath(DIR)
        end
        png(ZZ_plot,filepath) #"S$(line)_$(sta)_$(comp_iter)_$(fmin)to$(fmax).png"
    end
end

# select example parameters 
name, stacktype, component, filter = "CHN", "linear", "ZZ", "NO.B4"
frequency_plots = [[0.1,0.2],[0.2,0.5],[0.5,1.]]
lw = 0.5 #Decrease line thickness by half from default for cleaner plots
scale = 0.1 # inversely related to plot amplitude - may need to scale accordingly
rootdir = "/Users/julianschmitt/Desktop/SeisPlots" # chose root directory for plots

# Open file ex 
fpath = "/Users/julianschmitt/Downloads/2017/2017_CI.CHN.h5"
file = h5open(fpath,"r")

# filter and extract correlations and plot
fcorrs = get_corrs(file, stacktype, component, filter)
vert_plot(fcorrs, name, component, filter, stacktype) 