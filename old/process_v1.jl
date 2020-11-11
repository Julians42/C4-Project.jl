T = @elapsed using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays, AWSS3

#coeffs
cc_step, cc_len = 3600, 3600
maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9
half_win, water_level = 30, 0.01

using Pkg 
ENV["GR"] = ""
Pkg.build("GR")

#Add procs to access multiple cores
addprocs()

function indexer(ar::Array{String,1}, stations::Array{String,1})
    """
        Returns Indices of stations and recievers (others) in ar using occursin
    """
    indices = Array{Int64,1}(undef,0)
    for i in 1:length(stations)
        for j in 1:length(ar)
            if occursin(stations[i],ar[j]) == true
                push!(indices, j)
            end
        end
    end
    not_indices = setdiff(1:length(ar), indices)
    return([indices,not_indices])
end

function plot_corrs(corrs::Array{CorrData,1})
    # Add filtered correlations to appropriate plots
    T = collect(-corrs[1].maxlag:1/corrs[1].fs:corrs[1].maxlag)
    for j in 1:length(frequency_plots)
        #Select frequency
        fmin = frequency_plots[j][1]
        fmax = frequency_plots[j][2]
        
        #Determine the y scale axis
        distances = []
        for i in 1:length(corrs)
            push!(distances, corrs[i].dist)
        end
        y_max = maximum(distances) +0.5

        # Define plots 
        EE_plot = plot(title = "EE Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)
        EN_plot = plot(title = "EN Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)
        EZ_plot = plot(title = "EZ Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)
        NN_plot = plot(title = "NN Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)
        NZ_plot = plot(title = "NZ Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)
        ZZ_plot = plot(title = "ZZ Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)

        for i in 1:length(corrs)
            if corrs[i].comp == "EE"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(EE_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif corrs[i].comp == "EN"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(EN_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif corrs[i].comp == "EZ"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(EZ_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif corrs[i].comp == "NN"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(NN_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif corrs[i].comp == "NZ"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(NZ_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif corrs[i].comp == "ZZ"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(ZZ_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            end
        end
        png(EE_plot,"plots/moveouts/EE_plot$(fmin)_$(fmax).png")
        png(EN_plot,"plots/moveouts/EN_plot$(fmin)_$(fmax).png")
        png(EZ_plot,"plots/moveouts/EZ_plot$(fmin)_$(fmax).png")
        png(NN_plot,"plots/moveouts/NN_plot$(fmin)_$(fmax).png")
        png(NZ_plot,"plots/moveouts/NZ_plot$(fmin)_$(fmax).png")
        png(ZZ_plot,"plots/moveouts/ZZ_plot$(fmin)_$(fmax).png")
    end
end
maxbin, bin_size = 300., 10. 
minfreq, maxfreq = 0.1, 0.2
#Heatmap plot function - could be improved, some arguments are inbedded
function heatmapper(corrs::Array{CorrData,1}, maxbin::Float64, bin_size::Float64, minfreq::Float64, maxfreq::Float64)
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
    C = corrs[1]
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

    maxdist = 300.
    ind = findall((bins .> 0.15) .& (bins .<= maxdist))
    Cmat = Cmat[:,ind]
    counts = counts[ind]
    bins = bins[ind]
    Cmat ./= counts'
    lags = -C.maxlag:1/C.fs:C.maxlag
    maxlag = 300.
    lagind = findall(abs.(lags) .<= maxlag)
    Cmat = Cmat[lagind,:]
    abs_max!(Cmat)
    # plot envolopes at 5 km/s and 1.3 km/s
    heatmap(
        lags[lagind],
        bins,
        Cmat',
        c=:balance,
        xlabel="Lag [s]",
        ylabel="Inter-station Distance [km]",
        legend = :none,
        dpi=500,
    )
    plot!([0.,80.],[0.1,400],line=(1.5,:dash),color=:black,label="")
    plot!([0.,-80.],[0.1,400],line=(1.5,:dash),color=:black,label="")
    plot!([0.,267.],[0.03,400],line=(1.5,:dash),color=:black,label="")
    plot!([0.,-267.],[0.03,400],line=(1.5,:dash),color=:black,label="")
    annotate!(2.5,17, Plots.text("5 km/s", 14, :dark, rotation = 78 ))
    annotate!(200,300, Plots.text("1.5 km/s", 14, :dark, rotation = 48 ))
    ylims!((minimum(bins),maximum(bins)))
    xlims!((-maxlag,maxlag))
    CORROUT = expanduser("plots/moveouts/")
    if isdir(CORROUT) == false
        mkpath(CORROUT)
    end
    png(expanduser("plots/moveouts/heatplot_test.png"))
end

function list_corrs(csv_filelist::String, sources::Array{String,1}, path::String)
    """
        Takes a seisbasin csv filename, and a list of sources 
        Returns an array of day-specific corr_data filepaths ready for download 
        via SeisNoise's ec2download fucntion
    """
    corr_info = CSV.read(csv_filelist) # Read in csv file 
    #println(corr_info)
    names = Array{String,1}(undef,0)
    for i in 1:size(corr_info)[1]
        location = corr_info.location[i]
        if ismissing(location) == true #Deal with empty location information
            location = ""
        end
        name = join([corr_info.network[i],corr_info.station[i], location,corr_info.channel[i]],".")
        push!(names, name)
    end
    #println(names)
    #Index for correlation pairs
    source_index = indexer(names,sources)
    len_src, len_rec = length(source_index[1]), length(source_index[2])

    # Create DataFrame to store the names of the correlations with components 
    corr_pairs = DataFrame(corr_name = String[], comp = String[]) #Pairs is length of sources * receivers
    for i in 1:len_src
        for j in 1:len_rec
            pth = strip(join(vcat(split(names[source_index[1][i]],".")[1:3],split(names[source_index[2][j]],".")[1:3]),"."),'.')
            comp = join([names[source_index[1][i]][end],names[source_index[2][j]][end]],"")
            push!(corr_pairs, Dict(:corr_name => "$pth", :comp => "$comp"))
        end
    end
    #println(corr_pairs)
    # Make file paths by joining pth and comp from corr_pairs and file routing 
    file_paths = Array{String, 1}(undef, length(corr_pairs.corr_name)); 
    for i in 1:length(corr_pairs.corr_name)
        pth = "corr_data/$(corr_pairs.corr_name[i])/$(corr_pairs.comp[i])/$(path)_$(corr_pairs.corr_name[i]).jld2"
        file_paths[i] = pth
    end
    #Return list of filepaths
    return file_paths
end

@everywhere begin
    using SCEDC, AWSCore, Dates 
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    startdate = "2017-01-26" # Select Start Date
    enddate = "2017-02-28" # Select End Date
    days = Date(startdate):Day(1):Date(enddate)
    network = "CI"
    channel = "BH?"
    OUTDIR = "~/data"
end


# Get list of CSVs to transfer corr_index/2017/2017_026_correlation_index.csv
list_csv = Array{String,1}(undef, length(days))
for (index, dy) in enumerate(days)
    year = Dates.year(dy)
    path = join([Dates.year(dy),lpad(Dates.dayofyear(dy),3,"0")],"_") # Yeilds "YEAR_JDY"
    list_csv[index] = "corr_index/$year/$(path)_correlation_index.csv"
end

#Download CSVs
@eval @everywhere list_csv = $list_csv # Retain SCEDC download functionality, parallel download not actually needed
ec2download(aws, bucket2, list_csv, "~/") # Retain seisbasin filepathing






source_station = "CJM"
t1 = Dates.now()
filelist = nothing
for (index, dy) in enumerate(days[1:7])
    year = Dates.year(dy)
    path = join([Dates.year(dy),lpad(Dates.dayofyear(dy),3,"0")],"_") # Yeilds "YEAR_JDY"
    next_filelist = list_corrs("/home/ubuntu/corr_index/$year/$(path)_correlation_index.csv",[source_station], path) #IPT, SVD
    if index == 1
        filelist =next_filelist
    else
        filelist = vcat(filelist, next_filelist)
    end
end
println(length(filelist))


function load_crap2(filelist::Array{String,1}) #Load corr filelist into array with error handling
    corrs = Array{CorrData,1}(undef,0)
    for (index, corr) in enumerate(filelist)
        try 
            elt = load_corr(corr, convert(String,split(corr,"/")[3]))
            push!(corrs, elt)
        catch 
            deleteat!(filelist, index)
        end
    end
    return corrs
end
#corrs = load_crap2(filelist)

println(length(filelist))

#join(split(filelist[i],"/")[1:3],"/")

# split elmts in filelist and return unique corr_pair and comps for stacking
truncated = Array{String,1}(undef,0)
for i in 1:length(filelist)
    elt = join(split(filelist[i],"/")[1:3],"/")
    push!(truncated, elt)
end
unq_path = unique(truncated)


@everywhere begin
    function cc_medianmute!(C::CorrData, cc_medianmute_α::Float64 = 10.0)
        C.corr, inds = cc_medianmute(C.corr, cc_medianmute_α)
        C.t = remove_medianmute(C, inds)
        return nothing
    end
    function cc_medianmute(A::AbstractArray, cc_medianmute_α::Float64 = 10.0)
        #1. compute median of maximum amplitude of all corrs
        T, N = size(A)
        cc_maxamp = vec(maximum(abs.(A), dims=1))
        cc_medianmax = median(cc_maxamp)
        inds = findall(x-> x <= cc_medianmute_α*cc_medianmax,cc_maxamp)
        #NOTE: you cannot unbind entire array, so remove_nanandzerocol! is not used here.
        return A[:, inds], inds
    end
    remove_medianmute(C::CorrData, inds) = (return C.t[inds])
    function load_sum_stack(paths::String, len::Float64)
        try
            directory = joinpath.(paths,readdir(paths))
            single_corrs = load_corr.(directory, convert(String,split(paths,"/")[end]))
            single_temp = sum(shorten.(single_corrs,len))
            #cc_medianmute!(single_temp,2.)
            pair_stack = SeisNoise.stack!(single_temp, allstack=true)
            return pair_stack
        catch e
            println(e)
        end
    end
end
corrs = pmap(x -> load_sum_stack(x, 300.), unq_path) #returns array of corr_data stacked over days

truncated2 = Array{String,1}(undef,0) # Get unique station pairs
for i in 1:length(filelist)
    elt = join(split(filelist[i],"/")[1:2],"/")
    push!(truncated2, elt)
end
unq_pair = unique(truncated2)

#sta_pairs = findall(x -> occursin(unq_pair[1], x), unq_path)# unq_path is pair + comp, while unq_pair is just station pairs

big_array = Array{Array{CorrData,1}}(undef,0) # nested array of corr stacks per station pair
for (index, pair) in enumerate(unq_pair)
    try
        sta_pairs = findall(x -> occursin(pair, x), unq_path)# unq_path is pair + comp, while unq_pair is just station pairs
        if length(sta_pairs)==9
            println(corrs[sta_pairs[1]].name)
            push!(big_array, corrs[sta_pairs])
        end
    catch e
        println(e)
    end
end

rotated = Array{Array{CorrData,1}}(undef,0) #rotate all corrs
for sta_group in big_array
    try
        rtt = SeisNoise.rotate(sta_group, sta_group[1].azi, sta_group[1].baz)
        push!(rotated, rtt)
    catch e
        println(e)
    end
end

processed_corrs = Array{CorrData,1}(undef,0) #all corrs, 1 array 
for elt in rotated
    for corr in elt
        push!(processed_corrs, corr)
    end
end

#Save processed corrs to disk (for transfer to bucket/local)
function save_corr3(C::CorrData, CORROUT::String)
    # check if CORRDIR exists
    CORROUT = expanduser(CORROUT)
    if isdir(CORROUT) == false
        mkpath(CORROUT)
    end

    #name is YEAR_JDY_PAIRNAME
    yr,j_day = Dates.year(Date(C.id)), lpad(Dates.dayofyear(Date(C.id)),3,"0")
    p_name = strip(join(deleteat!(split(C.name,"."),[4,8]),"."),'.')
    name = join([yr,j_day,p_name,C.comp],"_")

    # create JLD2 file and save correlation
    filename = joinpath(CORROUT,"$(name).jld2")
    file = jldopen(filename, "a+")
    if !(C.comp in keys(file))
        group = JLD2.Group(file, C.comp)
        group[C.id] = C
    else
        file[C.comp][C.id] = C
    end
    close(file)
end
for i in 1:length(processed_corrs) # Save all correlations 
    save_corr3(processed_corrs[i], "processed/$source_station")
end
t2 = Dates.now()
t_diff = t2-t1
println("Corr Processing of $(length(processed_corrs)) completed in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(t_diff)))).")


scp -i /Users/julianschmitt/Downloads/my_aws_key.pem -r ubuntu@ec2-54-187-31-190.us-west-2.compute.amazonaws.com:/home/ubuntu/processed/CJM /Users/julianschmitt/Downloads/corr_data/
# transfer
aws s3 cp s3://seisbasin/corr_data /home/ubuntu/corr_data --recursive  --exclude "*"  --include "CI.CJM*NO.B4*.jld2"
