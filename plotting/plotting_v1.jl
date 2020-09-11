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
source_station = "CJM"
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
        if occursin("NO.B4", corr.name)
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

# Plotting function for node data vertical plots 
function plot_1comp(corrs_passed::Array{CorrData,1})
    # Add filtered correlations to appropriate plots
    T = collect(-corrs_passed[1].maxlag:1/corrs_passed[1].fs:corrs_passed[1].maxlag)
    for j in 1:length(frequency_plots)
        #Select frequency
        fmin = frequency_plots[j][1]
        fmax = frequency_plots[j][2]
        # Define plots 
        ZZ_plot = plot(xlims = (0,80), ylims = (0,130), 
                        yticks=[],xlabel = "Time (s)", ylabel = "South (Bottom) to North (Top)", 
                        xtickfontsize=5,ytickfontsize=5,fontsize=5, xguidefontsize = 10, yguidefontsize = 10,
                        legendfontsize = 15)
        
        # Adjust scaling by comparison by first corr
        Cstack = bandpass(SeisNoise.stack(corrs_passed[1]), fmin, fmax)
        id_mid = round(Int, size(Cstack.corr, 1))
        println(id_mid)
        scaled = 0.15*maximum(broadcast(abs, shorten(Cstack,200.).corr))
        println(scaled)
        max_amplitudes = Array{Float64, 1}(undef,0)
        # corr_processed = Array{CorrData, 1}(undef, 0)
        # for (index, corr) in enumerate(corrs_passed)
        #     processed_comp = Cstack = bandpass(SeisNoise.stack(corrs_passed[i]), fmin, fmax).corr/scaled
        plot_process = []
        for i in 1:length(corrs_passed)
            comp = corrs_passed[i].comp
            if comp == component
                Cstack = bandpass(SeisNoise.stack(corrs_passed[i]), fmin, fmax)
                push!(max_amplitudes, maximum(shorten(Cstack,200.).corr/scaled))
                push!(plot_process, Cstack.corr/scaled .+5*i)
                println(maximum(Cstack.corr / scaled))
                #plot!(ZZ_plot, T, Cstack.corr / scaled .+ 5*i, fmt = :png, linewidth = lw, reuse = false, legend = false)
            end
        end
        title!("ZZ", fontsize=5)
        plot!(ZZ_plot, -T, plot_process, color=:cm_maxamp, colorbar_title="Normalized maximum amplitude", line_z=max_amplitudes', fmt = :png, linewidth = lw, reuse = false, legend = false)
        plot!(size=(250,400),dpi=1000)
        png(ZZ_plot,"/Users/julianschmitt/Desktop/SeisPlots/nodenew_$(component)_$(fmin)to$(fmax).png")
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
        
        #Determine the y scale axis
        # distances = []
        # for i in 1:length(corrs)
        #     push!(distances, corrs[i].dist)
        # end
        y_max = 300.
        # Define plots 
        RR_plot = plot(title = "RR Plot", xlims = (-300,300), ylims = (-1,y_max), color = :black)
        RT_plot = plot(title = "RT Plot", xlims = (-300,300), ylims = (-1,y_max), color = :black)
        RZ_plot = plot(title = "RZ Plot", xlims = (-300,300), ylims = (-1,y_max), color = :black)
        TR_plot = plot(title = "TR Plot", xlims = (-300,300), ylims = (-1,y_max), color = :black)
        TT_plot = plot(title = "TT Plot", xlims = (-300,300), ylims = (-1,y_max), color = :black)
        TZ_plot = plot(title = "TZ Plot", xlims = (-300,300), ylims = (-1,y_max), color = :black)
        ZR_plot = plot(title = "ZR Plot", xlims = (-300,300), ylims = (-1,y_max), color = :black)
        ZT_plot = plot(title = "ZT Plot", xlims = (-300,300), ylims = (-1,y_max), color = :black)
        ZZ_plot = plot(title = "ZZ Plot", xlims = (-300,300), ylims = (-1,y_max), color = :black)


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
        big_plot = plot(RR_plot, RT_plot, RZ_plot, TR_plot, TT_plot, TZ_plot, ZR_plot, ZT_plot, ZZ_plot, layout=9)#, title="Moveout Plots filtered $(fmin)-$(fmax) Hz")
        png(big_plot,"/Users/julianschmitt/Downloads/Seismo/big_plot$(fmin)_$(fmax).png")
    end
end