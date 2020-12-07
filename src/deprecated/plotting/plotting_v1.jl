# Install Packages and build GR for plotting
using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays, ColorSchemes, Plots.PlotUtils

using Pkg 
ENV["GR"] = ""
Pkg.build("GR")


#coeffs
cc_step, cc_len = 3600, 3600
maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9
half_win, water_level = 30, 0.01

# Select Source Station
source_station = "SVD"
component = "TT" # pick out just node data with "XX" component



# Select plotting parameters
frequency_plots = [[0.1,0.5],[0.5,1.0],[0.1,1.]] #[[0.2,0.3],[0.3,0.4],[0.4,0.5]]  #[[0.5,1.],[1.,2.],[0.1,0.2]]
lw = 0.5 #Decrease line thickness by half from default

#read file names and load corrs
#pth = "/Users/julianschmitt/Downloads/corr_data/2017/$(source_station)"
pth = "/Users/julianschmitt/Downloads/corr_data/2017/CJM"
corr_data = joinpath.(pth,readdir(pth)[2:end])
# Load corrs into array
corrs = Array{CorrData,1}(undef,0)
for elt in corr_data
    corr = load_corr(elt, elt[end-6:end-5])
    push!(corrs,corr)
end
println(corrs)


# Returns corrs which are of component comp
function comp_corrs(corrs::Array{CorrData,1}, comp::String)
    corrs_comp = Array{CorrData,1}(undef,0)
    for (index, corr) in enumerate(corrs)
        if corr.comp == comp
            push!(corrs_comp, corr)
        end
    end
    return corrs_comp
end
#Returns B40XX node data only
function filter_nodes(corrs::Array{CorrData,1})
    corrs_nodes = Array{CorrData,1}(undef,0)
    for (index, corr) in enumerate(corrs)
        if occursin("NO.B4", corr.name) # Select node series
            push!(corrs_nodes, corr)
        end
    end
    return corrs_nodes
end

node_data = filter_nodes(comp_corrs(corrs, component)); # Do filtering


# Plot and test custom color scheme - to include lighter colors lower 50
σs = 0.5:0.1:2
normal_x = -5:0.01:5
normal_y = [exp.(-normal_x.^2 / (2σ^2)) / (2π * σ^2) for σ in σs]
loadcolorscheme(:cm_maxamp,ColorSchemes.gist_heat.colors[end-50:-1:1], "maxamp color", "for waveform plot")
plot(normal_x, normal_y, color=:cm_maxamp, line_z = σs')

# VERTICAL PLOTS: Plotting function for node data 
function plot_1comp(corrs_passed::Array{CorrData,1}, sta::String, comp_iter::String)
    # Add filtered correlations to appropriate plots
    T = collect(-corrs_passed[1].maxlag:1/corrs_passed[1].fs:corrs_passed[1].maxlag)
    y_max = length(corrs_passed)*5+9
    for j in 1:length(frequency_plots)
        #Select frequency
        fmin = frequency_plots[j][1]
        fmax = frequency_plots[j][2]
        # Define plots 
        ZZ_plot = plot(xlims = (0,100), ylims = (0,y_max), 
                        yticks=[],xlabel = "Time (s)", ylabel = "South (Bottom) to North (Top)", 
                        xtickfontsize=5,ytickfontsize=5,fontsize=5, xguidefontsize = 10, yguidefontsize = 10,
                        legendfontsize = 15)
        
        # Adjust scaling by comparison by first corr
        Cstack = bandpass(SeisNoise.stack(corrs_passed[1]), fmin, fmax)
        id_mid = round(Int, size(Cstack.corr, 1))
        #println(id_mid)
        scaled = 0.15*maximum(broadcast(abs, shorten(Cstack,200.).corr))
        #println(scaled)
        max_amplitudes = Array{Float64, 1}(undef,0)
#         corr_processed = Array{CorrData, 1}(undef, 0)
#         for (index, corr) in enumerate(corrs_passed)
#             processed_comp = Cstack = bandpass(SeisNoise.stack(corrs_passed[i]), fmin, fmax).corr/scaled
        plot_process = []
        for i in 1:length(corrs_passed)

            Cstack = bandpass(SeisNoise.stack(corrs_passed[i]), fmin, fmax)
            push!(max_amplitudes, maximum(shorten(Cstack,200.).corr/scaled))
            push!(plot_process, Cstack.corr/scaled .+5*i)
            #println(maximum(Cstack.corr / scaled))
                #plot!(ZZ_plot, T, Cstack.corr / scaled .+ 5*i, fmt = :png, linewidth = lw, reuse = false, legend = false)
        end
        title!("$(sta) G1 $(comp_iter)", fontsize=5)
        plot!(ZZ_plot, -T, plot_process, color=:cm_maxamp, colorbar_title="Normalized Maximum Amplitude", line_z=max_amplitudes', fmt = :png, linewidth = lw, reuse = false, legend = false)
        plot!(size=(250,400),dpi=1000)
        png(ZZ_plot,"/Users/julianschmitt/Desktop/SeisPlots/nodeG1_$(sta)_$(comp_iter)_$(fmin)to$(fmax).png")
    end
end

# Loop for VERTICAL PLOTS
# source_stations = ["SVD", "IPT","TA2"]
# components = ["RR","TT","ZZ"] # pick out just node data with "XX" component
source_stations = ["SVD"]
components = ["ZZ"]
for sta in source_stations
    # Load Corrs
    pth = "/Users/julianschmitt/Downloads/corr_data/2017/$(sta)"
    corr_data = joinpath.(pth,readdir(pth)[2:end])
    # Load corrs into array
    corrs = Array{CorrData,1}(undef,0) # Load can be moved outside of only doing one station
    for elt in corr_data
        corr = load_corr(elt, elt[end-6:end-5])
        push!(corrs,corr)
    end
    for i in 1:length(components)
        comp_iter = components[i]
        #Filter correct data
        node_data2 = filter_nodes(comp_corrs(corrs, comp_iter)); 
        println(length(node_data2))
        # Plot 
        plot_1comp(node_data2, sta, comp_iter)
    end
end








# Messy but working 3 by 3 moveout plot
function plot_corrs(corrs::Array{CorrData,1})
    # Add filtered correlations to appropriate plots
    T = collect(-corrs[1].maxlag:1/corrs[1].fs:corrs[1].maxlag)
    for j in 1:length(frequency_plots)
        #Select frequency
        fmin = frequency_plots[j][1]
        fmax = frequency_plots[j][2]
        
        y_max = 103.
        lim_x = (-100,100)
        # Define plots 
        RR_plot = plot(title = "RR", xlims = lim_x, ylims = (-1,y_max), color = :black)
        RT_plot = plot(title = "RT", xlims = lim_x, ylims = (-1,y_max), color = :black)
        RZ_plot = plot(title = "RZ", xlims = lim_x, ylims = (-1,y_max), color = :black)
        TR_plot = plot(title = "TR", xlims = lim_x, ylims = (-1,y_max), color = :black)
        TT_plot = plot(title = "TT", xlims = lim_x, ylims = (-1,y_max), color = :black)
        TZ_plot = plot(title = "TZ", xlims = lim_x, ylims = (-1,y_max), color = :black)
        ZR_plot = plot(title = "ZR", xlims = lim_x, ylims = (-1,y_max), color = :black)
        ZT_plot = plot(title = "ZT", xlims = lim_x, ylims = (-1,y_max), color = :black)
        ZZ_plot = plot(title = "ZZ", xlims = lim_x, ylims = (-1,y_max), color = :black)


        for i in 1:length(corrs)
            comp = corrs[i].comp
            if comp == "RR"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(RR_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif comp == "RT"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(RT_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif comp == "RZ"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(RZ_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif comp == "TR"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(TR_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif comp == "TT"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(TT_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif comp == "TZ"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(TZ_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif comp == "ZR"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(ZR_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif comp == "ZT"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(ZT_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif comp == "ZZ"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(ZZ_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            end
        end
        big_plot = plot(RR_plot, RT_plot, RZ_plot, TR_plot, TT_plot, TZ_plot, ZR_plot, ZT_plot, ZZ_plot, layout=9, dpi = 1000)#, title="Moveout Plots filtered $(fmin)-$(fmax) Hz")
        png(big_plot,"/Users/julianschmitt/Desktop/SeisPlots/bigmoveouts/SVD$(fmin)_$(fmax).png")
    end
end


# Careful Heatplot applies bandpass to the corrs passed!
function heatmapper(corrs2::Array{CorrData,1}, maxbin::Float64, bin_size::Float64, minfreq::Float64, maxfreq::Float64)
    """
        Function takes an array of correlation data, filters the data between minfreq and maxfreq

        maxbin is the cutoff for the maximum distance between correlations (in km), larger distances are discarded,
            the bin_size intuitively gives the range of correlation distances which are stacked together.

        A good starting place to check visibility of large-scale correlations is: 
            maxbin, bin_size = 600., 10. 
            minfreq, maxfreq = 0.1, 0.2
    """
    bins = collect(0.:bin_size:maxbin-bin_size)
    bins1 = collect(bin_size:bin_size:maxbin)
    counts = zeros(length(bins))
    # freqmin = 1.
    # freqmax = 2.
    C = corrs2[1]
    Cmat = zeros(eltype(C.corr),size(C.corr,1),length(bins))
    for ii = 1:length(corrs)
        C = corrs[ii]
        ind = 0
        try
            ind = findall((C.dist .>= bins) .& (C.dist .<= bins1))[1]
        catch
            continue
        end
        bandpass!(C, minfreq,maxfreq)
        clean_up!(C,minfreq,maxfreq)
        abs_max!(C)
        Cmat[:,ind] .+= C.corr[:]
        counts[ind] += 1
        #println("Reading file $ii")
    end

    maxdist = maxbin # maxdist and maxbin are both 100. 
    maxlag = 300.
    ind = findall((bins .> 0.15) .& (bins .<= maxdist))
    Cmat = Cmat[:,ind]
    counts = counts[ind]
    bins = bins[ind]
    Cmat ./= counts'
    lags = -C.maxlag:1/C.fs:C.maxlag
    lagind = findall(abs.(lags) .<= maxlag)
    Cmat = Cmat[lagind,:]
    abs_max!(Cmat)
    # plot envolopes at 5 km/s and 1.3 km/s
    heatmap(
        lags[lagind],
        bins,
        Cmat',
        c=:balance,
        title="SVD ZZ Heatmap",
        xlabel="Lag [s]",
        ylabel="Inter-station Distance [km]",
        legend = :none,
        dpi=500,
    )
    plot!([0.,80.],[0.1,400],line=(1.5,:dash),color=:black,label="")
    plot!([0.,-80.],[0.1,400],line=(1.5,:dash),color=:black,label="")
    plot!([0.,267.],[0.03,400],line=(1.5,:dash),color=:black,label="")
    plot!([0.,-267.],[0.03,400],line=(1.5,:dash),color=:black,label="")
    annotate!(12,85, Plots.text("5 km/s", 14, :dark, rotation = 78 ))
    annotate!(63,85, Plots.text("1.5 km/s", 14, :dark, rotation = 60 ))
    ylims!((minimum(bins),maximum(bins)))
    #xlims!((-maxlag,maxlag))
    xlims!((-100,100))
    CORROUT = expanduser("plots/moveouts/")
    if isdir(CORROUT) == false
        mkpath(CORROUT)
    end
    png(expanduser("/Users/julianschmitt/Desktop/SeisPlots/heatplot_SVD_$(minfreq)_to_$(maxfreq).png"))
end