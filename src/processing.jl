export fftderivative, loc_max_amp, create_corr, get_corrs, mov_avg, get_node_files,
        plot_NS_dist, plot_tseries, plot_simple_nodes, wrap_plotter

# extracting and processing correlations 
function fftderivative(A::CorrData)
    """ Computes discrete derivative via FFT"""
    C = deepcopy(A)
    detrend!(C)
    taper!(C)
    F = fft(C.corr,1)
    F .*= fftfreq(length(C.corr)).* 1im .* 2π
    C.corr = real.(ifft(F,1))/(C.fs^2)
    return C
end
# gets location and maximum amplitude for corrs in array
function loc_max_amp(corr::CorrData, filt::Array{Float64,1}=[0.1,0.5])
    """ Computes essential correlation metrics:
    latitude, longitude, max amplitude, max amplitude for anticausal side,
     max amplitude for causal side, SNR (signal-to-noise ratio), SNR (anticausal side), 
        SNR (causal), correlation name
    """
    lat, lon = corr.loc.lat, corr.loc.lon
    bpass_corr = bandpass(fftderivative(fftderivative(corr)),filt[1],filt[2])
    half_len = convert(Int64,round(length(bpass_corr.corr)/2))
    max_ampl = maximum(abs.(bpass_corr.corr))
    max_ampl_anti = maximum(abs.(bpass_corr.corr[1:half_len]))
    max_ampl_causal = maximum(abs.(bpass_corr.corr[half_len:end]))
    snr = maximum(abs.(bpass_corr.corr))/std(bpass_corr.corr)
    snr_anti = maximum(abs.(bpass_corr.corr[1:half_len]))/std(bpass_corr.corr[1:half_len])
    snr_causal = maximum(abs.(bpass_corr.corr[half_len:end]))/std(bpass_corr.corr[half_len:end])
    return (lat, lon, max_ampl, max_ampl_anti, max_ampl_causal, snr, snr_anti, snr_causal, corr.name)
end
function create_corr(file::HDF5File, station::String, stacktype::String, component::String)
    """ Returns correlation with metadata from stacked HDF5 file"""
    try
        S = CorrData()
        # get corr data
        corr = file[station][component][stacktype][1:end]
        S.corr = reshape(corr,size(corr,1),1)
        S.name = station
        S.dist = read(file[station]["meta"]["dist"])
        S.loc.lat = read(file[station]["meta"]["lat"])
        S.loc.lon = read(file[station]["meta"]["lon"])
        S.loc.el = read(file[station]["meta"]["el"])
        S.azi = read(file[station]["meta"]["azi"])
        S.baz = read(file[station]["meta"]["baz"])
        S.fs = 20. # For the 20 Hz product
        return S
    catch
        return nothing
    end
end
# get information for a set of filtered stations, default to all stations
function get_corrs(file::HDF5File, stacktype::String="linear",
    component::String="ZZ", filter::String="")
    """Returns array of correlations from specified filters"""
    filtered = [name for name in names(file) if occursin(filter, name)] 
    corrs = [create_corr(file, station, stacktype, component) for station in filtered]
    corrs = [corr for corr in corrs if !isnothing(corr)] # strip nothings - catches meta file etc which errors
    return corrs
end 


# plotting stuff
function mov_avg(df::DataFrame, range::Int64=30) 
    """ Returns dataframe with `moving_avg` column box function filter"""
    df_new = df[:,["date","max_amp"]]
    amps = df_new.max_amp
    avg, len = similar(amps), length(amps)
    for i =1:len
        first = max(1, convert(Int64, round(i-range/2)))
        last = min(len, convert(Int64, round(i+range/2)))
        window = amps[first:last]
        sd, mn = std(window), mean(window)
        avg[i] = mean(filter(x -> (x < mn + 3*sd) && (x > mn-3*sd), window))
    end
    df_new.moving_avg = avg
    return df_new
end
function get_node_files(source::String, deployment_filt::String, int_range=Array{Int64, 1})
    return filter(x -> occursin(source, x) & occursin(deployment_filt,x) & (parse(Int64, split(x, ".")[end-2][end-1:end]) in collect(int_range[1]:int_range[2])), nodes)
end

function plot_NS_dist(corrs::Array{CorrData,1}, freqs::Array{Array{Float64,1},1}, plot_sorttype::String, 
    attr::Array{String,1}=["","","",""])
    lw, scale = 0.5, 3
    sort_method, p = [], []
    if plot_sorttype == "Latitude"
    sort_method = sort([corr.loc.lat for corr in corrs])
    p = sortperm([corr.loc.lat for corr in corrs]) # sort based on latitude 
    elseif plot_sorttype == "Longitude"
    sort_method = sort([corr.loc.lon for corr in corrs])
    p = sortperm([corr.loc.lon for corr in corrs]) # sort based on longitude
    else # sort by distance 
    plot_sorttype == "distance"
    println("Distance is not currently a supported stacktype")
    end

    lmin, lmax = minimum(sort_method), maximum(sort_method)
    new_corrs, n_corrs = corrs[p], length(corrs)
    y_labels = (round.(collect(range(0, 5*n_corrs, length=10)),digits=3),round.(collect(range(lmin, lmax, length=10)),digits=3))
    for freq_pair in freqs
        NS_plot = plot(xlims = (0,100), ylims = (-5,n_corrs*5+5), 
                yticks=y_labels,xlabel = "Time (s)", ylabel = "Station Nodes by $plot_sorttype",  
                xtickfontsize=5,ytickfontsize=5,fontsize=5, xguidefontsize = 10, yguidefontsize = 10,
                legendfontsize = 10, title="$(attr[3]) $(attr[2]) from $(attr[1]): $(freq_pair[1]) to $(freq_pair[2])",
                titlefontsize= 10)
        # bandpass correlations 
        corrs_processed, max_amp = [], []
        for corr in new_corrs
        # bandpass, scale, shift, and add to plot
        processed = bandpass(corr, freq_pair[1], freq_pair[2]).corr ./ scale
        push!(max_amp, maximum(processed[1:convert(Int64, round(length(processed)/2))]))
        if plot_sorttype == "Latitude"
            push!(corrs_processed, processed .+ ((corr.loc.lat -lmin)/(lmax-lmin)*n_corrs*5))
        else
            push!(corrs_processed, processed .+ ((corr.loc.lon -lmin)/(lmax-lmin)*n_corrs*5))
        end
    end

    # plotting
    T = collect(-300.:1/new_corrs[1].fs:300.)
    plot!(NS_plot, -T, corrs_processed, color=:cm_maxamp, colorbar_title="Normalized Maximum Amplitude", 
        line_z=max_amp', fmt = :png, linewidth = 0.5, reuse = false, legend = false)    
    plot!(size=(250,400),dpi=500)
    filepath = joinpath(rootdir,
    "stack_plots/$(n_deriv_names[n_derivatives+1])_$(attr[4])/$(attr[1])/$(attr[3])_from_$(attr[1])_$(attr[2])_$(freq_pair[1])to$(freq_pair[2]).png")
    # ensure filepath is valid 
    DIR = dirname(filepath)
    if !isdir(DIR)
    mkpath(DIR)
    end
    png(NS_plot,filepath)
    end
end

function plot_tseries(perm_station::String, node_stations_::Array{String,1}, source_, dir::String)
    test = DataFrame(CSV.File(perm_station))
    #get moving averages and naming schemes
    times, amplitudes = test.date, test.max_amp
    times_avg, amps_avg = mov_avg(test,90).date, mov_avg(test,90).moving_avg
    perm_name = join(split(split(perm_station,"/")[end],".")[3:4],".")
    node_range = join([split(node_stations_[1],".")[end-2],split(node_stations_[end],".")[end-2][end-1:end]],"-")
    name = "$perm_name and $node_range Peak Amplitude Time Series"
    
    # Plot permanent station moving average
    tseries = plot(times_avg, amps_avg, title = name, label="$perm_name Mov Avg", color="gray", 
        linewidth=2, xtickfontsize=10, xlabel="Time", ylabel="Correlation Velocity Max Amp - filt: 0.1-0.5 Hz",
        legend = :outertopright, size=(900, 400), dpi=300, leftmargin= 7mm, bottommargin=5mm, linealpha=.6)
    
    # plot nodes
    for (ind, file) in enumerate(node_stations_)
        node_name, color = join(split(split(file,"/")[end],".")[3:4],"."), ""
        if occursin("NO.B1", node_name) # color B1 nodes green and intersecting nodes blue
            color = "green"
        else
            color = "blue"
        end
        test2 = DataFrame(CSV.File(file))
        test3 = mov_avg(test2, 20)
        times2, amplitudes2 = test2.date, test2.max_amp
        times3, amplitudes3 = test3.date, test3.moving_avg
        plot!(tseries, times3, amplitudes3, color=color, linealpha=.6, linewidth=2, label="$node_name Mov Avg")
    end
    save_path = joinpath(plot_root, "tseries_combined/$dir/CI.$(source_)_$(perm_name)_to_nodes.png")
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end
    png(save_path)
end
function plot_simple_nodes(stations::Array{String,1}, source::String, vert_name::String, dir::String)
    name = "CI.$source Peak Amplitude Time Series at Intersect - $vert_name and B1"
    tseries = plot(title = name,linewidth=3, xtickfontsize=10, xlabel="Time", 
        ylabel="Correlation Velocity Max Amp - filt: 0.1-0.5 Hz",
        legend = :outertopright, size=(900, 400), dpi=300, leftmargin= 7mm, bottommargin=5mm)
    for (ind, file) in enumerate(stations)
        rec_name = join(split(split(file,"/")[end],".")[3:4],".")
        series = DataFrame(CSV.File(file))
        series_m = mov_avg(series, 5)
        times, amplitudes = series.date, series.max_amp
        times_m, amplitudes_m = series_m.date, series_m.moving_avg
        color=""
        if occursin("B1", rec_name)
            color = "green"
        else
            color = "blue"
        end
        println(rec_name)
        plot!(tseries, 1:length(amplitudes_m), amplitudes_m, color=color, linealpha=.6, 
                linewidth=2, label="$rec_name Mov Avg")
    end
    save_path = joinpath(plot_root,"tseries_combined/$dir/CI.$(source)_B1_$(vert_name)_intersect.png")
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end
    png(save_path)
end
function wrap_plotter(vnodes::Array{String,1}, bnodes::Array{String,1}, reciever::String, vert_name, dir)
    nodes = vcat(vnodes, bnodes)
    for source_ in ["CHN","CJM","DEV","IPT","LPC","SNO","SVD","TA2"]
        try
            nodes_source = filter(x -> occursin(source_, split(x, ".")[2]), readdir("Tseries_nodes"))
            node_intersect = filter(x-> any(occursin.(join(split(x,".")[3:4],"."), nodes)), nodes_source)

            plot_simple_nodes(joinpath.("Tseries_nodes",node_intersect), source_, vert_name, dir)
            rec = "Tseries/CI.$(source_).$reciever.csv"
            println(rec)
            nodes = ["Tseries_nodes/CI.$(source_).$(node).jld2.csv" for node in nodes]
            println(nodes)
            plot_tseries(rec, joinpath.("Tseries_nodes", node_intersect), source_, dir)
        catch e
            println(e)
        end
    end
end