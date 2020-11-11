T = @elapsed using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays, AWSS3, ColorSchemes, Plots.PlotUtils

#coeffs
cc_step, cc_len = 3600, 3600
maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9
half_win, water_level = 30, 0.01

# Select plotting parameters
frequency_plots = [[0.1,0.5],[0.5,1.0],[0.1,1.]] #[[0.2,0.3],[0.3,0.4],[0.4,0.5]]  #[[0.5,1.],[1.,2.],[0.1,0.2]]
lw = 0.5 #Decrease line thickness by half from default

σs = 0.5:0.1:2
normal_x = -5:0.01:5
normal_y = [exp.(-normal_x.^2 / (2σ^2)) / (2π * σ^2) for σ in σs];
loadcolorscheme(:cm_maxamp,ColorSchemes.gist_heat.colors[end-50:-1:1], "maxamp color", "for waveform plot");

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

function index_download(files::Array{String,1}, wanted_stations::Array{String,1})
    """
        Returns Indices of stations and recievers (others) in ar using occursin
    """
    source_stations = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV"
                        ,"VINE", "ROPE", "ARNO", "LUCI", "ROUF", "KUZD", "ALLI"]
    all_indices = Array{Int64,1}(undef,0)
    wanted_indices = Array{Int64,1}(undef,0)
    for i in 1:length(files)
        if any(occursin.(source_stations, files[i])) == true
            push!(all_indices, i)
        end
        if any(occursin.(wanted_stations, files[i])) == true
            push!(wanted_indices, i)
        end
    end
    reciever_indices = setdiff(1:length(files), all_indices)
    return([wanted_indices,reciever_indices])
end

function list_corrs(csv_filelist::String, sources::Array{String,1}, path::String) # sources ex: ["IPT"]
    """
        Takes a seisbasin csv filename, and a list of sources 
        Returns an array of day-specific corr_data filepaths ready for download 
        via SeisNoise's ec2download fucntion
    """
    corr_info = CSV.read(csv_filelist) # Read in csv file 
    names = Array{String,1}(undef,0)
    for i in 1:size(corr_info)[1]
        location = corr_info.location[i]
        if ismissing(location) == true #Deal with empty location information
            location = ""
        end
        name = join([corr_info.network[i],corr_info.station[i], location,corr_info.channel[i]],".")
        push!(names, name)
    end
    #Index for correlation pairs - sources MUST be the same as sources used during correlations, selected sources should be a subset array 
    source_index = index_download(names, sources)

    # Create DataFrame to store the names of the correlations with components 
    corr_pairs = DataFrame(corr_name = String[], comp = String[]) #Pairs is length of sources * receivers
    for (i, source) in enumerate(source_index[1])
        for (j, reciever) in enumerate(source_index[2])
            pth = strip(join(vcat(split(names[source],".")[1:3],split(names[reciever],".")[1:3]),"."),'.')
            comp = join([names[source][end],names[reciever][end]],"")
            push!(corr_pairs, Dict(:corr_name => "$pth", :comp => "$comp"))
        end
    end
    # Make seisbasin file path
    file_paths = ["corr_data/$(corr_pairs.corr_name[i])/$(corr_pairs.comp[i])/$(path)_$(corr_pairs.corr_name[i]).jld2" for i in 1:length(corr_pairs.corr_name)]
    #Return list of filepaths
    return file_paths
end
@everywhere begin
    using SeisNoise, SeisIO, Statistics
    function load_filelist(file::String) # load correlations
        try 
            elt = load_corr(file, convert(String,split(file,"/")[3]))
            bool = true
            return [elt, bool]
        catch e
            print(e)
            bool = false
            return [nothing, bool]
        end
    end
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
            single_corrs = [corr for corr in single_corrs if !isnothing(corr)]
            single_temp = sum(shorten.(single_corrs,len))
            #cc_medianmute!(single_temp,2.)
            pair_stack = SeisNoise.stack!(single_temp, allstack=true)
            return pair_stack
        catch e
            println(e)
        end
    end
    function load_sum_stack_robust(paths::String, len::Float64)
        try
            directory = joinpath.(paths,readdir(paths))
            single_corrs = load_corr.(directory, convert(String,split(paths,"/")[end]))
            single_corrs = [corr for corr in single_corrs if !isnothing(corr)]
            single_temp = sum(shorten.(single_corrs,len))
            #cc_medianmute!(single_temp,2.)
            pair_stack = SeisNoise.robuststack!(single_temp)
            return pair_stack
        catch e
            println(e)
        end
    end
    function load_sum_stack_pws(paths::String, len::Float64)
        try
            directory = joinpath.(paths,readdir(paths))
            single_corrs = load_corr.(directory, convert(String,split(paths,"/")[end]))
            single_corrs = [corr for corr in single_corrs if !isnothing(corr)]
            single_temp = sum(shorten.(single_corrs,len))
            #cc_medianmute!(single_temp,2.)
            pair_stack = SeisNoise.pws!(single_temp)
            return pair_stack
        catch e
            println(e)
        end
    end
end


#### add in some error checking bro
function s3_file_map(aws::AWSConfig,bucket::String,filein::String,fileout::String)
    try
        s3_get_file(aws, bucket, filein, fileout)
        println("Downloading file: $filein       \r")
    catch 
        print("File $filein does not exist!")
    end
	return nothing
end

@everywhere begin
    using SCEDC, AWSCore, Dates 
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    startdate = "2017-06-30" # Select Start Date
    enddate = "2017-12-31" # Select End Date
    days = Date(startdate):Day(1):Date(enddate)
end

#@everywhere begin # helper functions for safe correlation download
 
"""
ec2download(aws,bucket,filelist,OUTDIR)
Download files using pmap from S3 to EC2.
# Arguments
- `aws::AWSConfig`: AWSConfig configuration dictionary
- `bucket::String`: S3 bucket to download from.
- `filelist::Array{String}`: Filepaths to download in `bucket`.
- `OUTDIR::String`: The output directory on EC2 instance.
# Keywords
- `v::Int=0`: Verbosity level. Set v = 1 for download progress.
- `XML::Bool=false`: Download StationXML files for request. Downloads StationXML to
    `joinpath(OUTDIR,"XML")`.
"""
function safe_download(aws::AWSConfig,seisbucket::String,file_list::Array{String},OUTDIR::String;
                        v::Int=0,XML::Bool=false)

    # check being run on AWS
    tstart = now()
    !localhost_is_ec2() && @warn("Running locally. Run on EC2 for maximum performance.")

    println("Starting Download...      $(now())")
    println("Using $(nworkers()) cores...")


    # create output files
    OUTDIR = expanduser(OUTDIR)
    outfiles = [joinpath(OUTDIR,f) for f in file_list]
    filedir = unique([dirname(f) for f in outfiles])
    for ii = 1:length(filedir)
        if !isdir(filedir[ii])
            mkpath(filedir[ii])
        end
    end

    # get XMLmesse
    if XML
        XMLDIR = joinpath(OUTDIR,"XML")
        getXML(aws,seisbucket,file_list,XMLDIR,v=v)
    end

    # send outfiles everywhere
    @eval @everywhere outfiles=$outfiles
    # do transfer to ec2
    startsize = diskusage(OUTDIR)
    if v > 0
        pmap(
            s3_file_map,
            fill(aws,length(outfiles)),
            fill(seisbucket,length(outfiles)),
            file_list,
            outfiles
        )
    else
        pmap(
            s3_get_file,
            fill(aws,length(outfiles)),
            fill(seisbucket,length(outfiles)),
            file_list,
            outfiles,
        )
    end

    println("Download Complete!        $(now())          ")
    tend = now()
    # check data in directory
    endsize = diskusage(OUTDIR)
    downloadsize = endsize - startsize
    downloadseconds = (tend - tstart).value / 1000
    println("Download took $(Dates.canonicalize(Dates.CompoundPeriod(tend - tstart)))")
    println("Download size $(formatbytes(downloadsize))")
    println("Download rate $(formatbytes(Int(round(downloadsize / downloadseconds))))/s")
    return nothing
end

function s3_file_map(aws::AWSConfig,seisbucket::String,filein::String,fileout::String)
    try
        s3_get_file(aws, seisbucket, filein, fileout)
        println("Downloading file: $filein       \r")
        return nothing
    catch e
        println("Unable to download file $filein")
    end
end

function s3_get_seed(
    aws::AWSConfig,seisbucket::String,
    filein::String,
    demean::Bool,
    detrend::Bool,
    msr::Bool,
    prune::Bool,
    rr::Bool,
    taper::Bool,
    ungap::Bool,
    unscale::Bool,
    resample::Bool,
    fs::Real,)
    f = s3_get(aws, seisbucket, filein)
    S = parseseed(f)

    # remove empty channels
    if prune == true
        prune!(S)
    end

    # Get list of channels with sane instrument codes
    CC = get_seis_channels(S)

    if msr == true
        @warn("Getting response not implemented yet.")
    end

    # unscale
    if unscale == true
        unscale!(S, chans=CC)
    end

    # Demean
    if demean == true
    demean!(S, chans=CC)
    end

    # Taper
    if taper == true
    taper!(S, chans=CC)
    end

    # Ungap
    if ungap == true
    ungap!(S, chans=CC)
    end

    # resample data
    if resample == true && fs != 0
        resample!(S, chans=CC, fs=fs)
    end

    # Remove response
    # need to implement attaching response
    if rr == true
    @warn("Removing response not implemented yet.")
    end

    return S
end

function formatbytes(bytes::Real, digits::Int=1)
    units = ["B", "KB", "MB", "GB", "TB","PB"]
    bytes = max(bytes,0)
    pow = Int(floor((bytes > 0 ? log(bytes) : 0) / log(1024)))
    pow = min(pow,length(units))
    powind = pow < length(units) ? pow + 1 : pow
    return string(round(bytes / 1024 ^ pow,digits=digits)) * units[powind]
end

function diskusage(dir)
    s = read(`du -k $dir`, String)
    kb = parse(Int, split(s)[1])
    return 1024 * kb
end

"""
parseseed(f)
Convert uint8 data to SeisData.
"""
function parseseed(f::AbstractArray)
    S = SeisData()
    SeisIO.SEED.parsemseed!(S,IOBuffer(f),SeisIO.KW.nx_new,SeisIO.KW.nx_add,false,0)
    return S
end
#end

# plotting functions
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
function filter_nodes(corrs::Array{CorrData,1}, corr_group)
    corrs_nodes = Array{CorrData,1}(undef,0)
    for (index, corr) in enumerate(corrs)
        if occursin("NO.$corr_group", corr.name)
            push!(corrs_nodes, corr)
        end
    end
    return corrs_nodes
end
function vert_plot(corrs_passed::Array{CorrData,1}, sta::String, comp_iter::String, line::String, stack_type::String)
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
        scaled = 0.15*maximum(broadcast(abs, shorten(Cstack,200.).corr))
        max_amplitudes = Array{Float64, 1}(undef,0)
        plot_process = []
        # Stack over days and add to plot
        for i in 1:length(corrs_passed)
            Cstack = bandpass(SeisNoise.stack(corrs_passed[i]), fmin, fmax)
            push!(max_amplitudes, maximum(shorten(Cstack,200.).corr/scaled))
            push!(plot_process, Cstack.corr/scaled .+5*i)
        end
        title!("$(sta) S$line $(comp_iter) using $stack_type", fontsize=5)
        plot!(ZZ_plot, -T, plot_process, color=:cm_maxamp, colorbar_title="Normalized Maximum Amplitude", 
            line_z=max_amplitudes', fmt = :png, linewidth = lw, reuse = false, legend = false)
        plot!(size=(250,400),dpi=1000)
        filepath = "~/stack_plots/$stack_type/S$(line)_$(sta)_$(comp_iter)_$(fmin)to$(fmax).png"
        # ensure filepath is valid 
        DIR = dirname(filepath)
        if !isdir(DIR)
            mkpath(DIR)
        end
        png(ZZ_plot,"stack_plots/$stack_type/S$(line)_$(sta)_$(comp_iter)_$(fmin)to$(fmax).png")
    end
end

############################# Index Download #####################################

# Get list of CSVs to transfer corr_index/2017/2017_026_correlation_index.csv
month_index_unq = unique([(Dates.year(d), Dates.monthname(d)) for d in days])
month_fnames = ["month_index/$(ind[1])/$(ind[2]).csv" for ind in month_index_unq]

#Download CSVs
@eval @everywhere month_fnames = $month_fnames # Retain SCEDC download functionality, parallel download not actually needed
safe_download(aws, bucket2, month_fnames, "~/") # Retain seisbasin filepathing

stations = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV"
            ,"VINE", "ROPE", "ARNO", "LUCI", "ROUF", "KUZD", "ALLI", "CHN", "USB", "Q0048"]


########################### Get unique source stations ###############

repeat_pairs = Array{String, 1}(undef, 0)
all_files = Array{String, 1}(undef, 0)
for name in month_fnames
    month_file = CSV.read(name)
    m_files = month_file["Files"]
    a_pth = month_file["paths"]
    for (ind, pair) in enumerate(m_files) # could use something like python's ravel, but not spending more time here
        push!(repeat_pairs, pair)
        push!(all_files, a_pth[ind])
    end
end

unq_sta_pairs = unique(repeat_pairs)


####################### Loop through station pairs ############################

for unq_pair in unq_station_pairs
    # Get correct jld2 files
    corr_files = all_files[occursin.(unq_pair, all_files)]
    # Transfer files
    safe_download(aws, bucket2, corr_files, "~/")
    # Read in correlations -- USE GLOB to get over years/months




############################# Loop through stations ############################
for station in stations
    source_station = [station]
    t1 = Dates.now()
    filelist = nothing
    for (index, dy) in enumerate(days)
        year = Dates.year(dy)
        path = join([Dates.year(dy),lpad(Dates.dayofyear(dy),3,"0")],"_") # Yeilds "YEAR_JDY"
        next_filelist = list_corrs("/home/ubuntu/corr_index/$year/$(path)_correlation_index.csv",source_station, path) 
        if index == 1
            filelist = next_filelist
        else
            filelist = vcat(filelist, next_filelist)
        end
    end
    if length(filelist) == 0 # if no correlations from that source station continue 
        continue
    end
    println("$(length(filelist)) correlations to be processed for station $(source_station[1])")

    safe_download(aws, bucket2, filelist, "~/") # data download

    ################################# Load Data and Index ###################################

    # correlations = joinpath.("/home/ubuntu/corrs", filelist)
    #  # load station files
    # corrs = pmap(x -> load_filelist(x), filelist)
    # bools = [corr[2] for corr in corrs]
    # corr_data = [corr[1] for corr in corrs]

    # split elmts in filelist and return unique corr_pair and comps for stacking
    truncated = Array{String,1}(undef,0) # unique stations and components 
    truncated2 = Array{String,1}(undef,0) # Get unique station pairs
    for i in 1:length(filelist)
        elt = join(split(filelist[i],"/")[1:3],"/")
        elt2 = join(split(filelist[i],"/")[1:2],"/")
        push!(truncated, elt)
        push!(truncated2, elt2)
    end
    unq_path = unique(truncated)
    unq_pair = unique(truncated2)

    ######################################### Stack Correlations ###################################
    println("Stacking correlations .....")
    corrs_stacked_mean = pmap(x -> load_sum_stack(x, 300.), unq_path) #returns array of corr_data stacked over days
    corrs_stacked_robust = pmap(x -> load_sum_stack_robust(x, 300.), unq_path) # for robust stacking 
    #corrs_stacked_pws = pmap(x -> load_sum_stack_pws(x, 300.), unq_path) # for phase-weighted stacking 
    println("Correlations stacked!")
    function corr_rotate(corrs_stacked::Array{CorrData,1}, sta_pair, sta_pair_comp)
        big_array = Array{Array{CorrData,1}}(undef,0) # nested array of corr stacks per station pair
        for (index, pair) in enumerate(sta_pair)
            try
                sta_pairs = findall(x -> occursin(pair, x), sta_pair_comp)# unq_path is pair + comp, while unq_pair is just station pairs
                if length(sta_pairs)==9
                    println(corrs_stacked[sta_pairs[1]].name)
                    push!(big_array, corrs_stacked[sta_pairs])
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

        processed_corrs = Array{CorrData,1}(undef,0) #flattens array - could be a one liner 
        for elt in rotated
            for corr in elt
                push!(processed_corrs, corr)
            end
        end
        return processed_corrs
    end
    println("Rotating correlations....")
    processed_mean = corr_rotate(corrs_stacked_mean, unq_pair, unq_path)
    processed_robust = corr_rotate(corrs_stacked_robust, unq_pair, unq_path)
    #processed_pws = corr_rotate(corrs_stacked_pws, unq_pair, unq_path)
    println("Correlations rotated!")
    ##################################### Save and Transfer Stacked Correlations ##############################

    #Save processed corrs to disk (for transfer to bucket/local)
    function save_processed(C::CorrData, CORROUT::String)
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
    types = ["mean", "robust", "pws"]
    for (ind, stacks) in enumerate([processed_mean, processed_robust])#, processed_pws])
        stack_type = types[ind]
        map(x -> save_processed(stacks, "processed/$stack_type/$source_station"))
    end
    #Transfer stacked correlations
    stacked_files = join.("processed/$stack_type/$source_station", readdir("processed/$stack_type/$source_station"))
    @elapsed pmap(x -> s3_put(aws, bucket2, "processed/$stack_type/$source_station", read(x)), stacked_files)

    ################################### Make and Transfer Node Plots ########################################
    println("Plotting nodal data now!")
    types = ["mean", "robust", "pws"]
    nodes = ["B4", "G1", "G2"]
    components = ["ZZ", "TT","RR"]
    for (ind, stacks) in enumerate([processed_mean, processed_robust])#, processed_pws])
        stack_type = types[ind]
        for (j,k) in Iterators.product(1:length(nodes), 1:length(components))
            node, component = nodes[j], components[k]
            filtered_nodes = filter_nodes(comp_corrs(stacks, component), node)
            vert_plot(filtered_nodes, source_station, component, node, stack_type)
        end
    end

    # Transfer node plots
    stacked_plots = glob("*/*$(source_station)*.png", "stack_plots") # Transfer plots from this station only
    @elapsed pmap(x -> s3_put(aws, bucket2, "stack_plots", read(x)), stacked_plots)
end