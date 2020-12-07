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

@everywhere begin
    using SCEDC, AWSCore, Dates, DataFrames, AWSS3
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    startdate = "2017-01-01" # Select Start Date
    enddate = "2017-06-30" # Select End Date
    days = Date(startdate):Day(1):Date(enddate)
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
df = CSV.File("files/all_locations_socal.csv") |> DataFrame! 
@eval @everywhere df = $df
job_name = "linear_2017_h5"
@eval @everywhere job_name = $job_name
############################# Get files for unique station pair ####################

to_download = unique(Iterators.flatten([DataFrame(CSV.File(file)).paths for file in month_fnames]))
safe_download(aws, bucket2, to_download, "~/")

function csv_merge(large_index = Array{String, 1})
    """Returns dataframe with unique station pairs and csvs containing that station pair"""
    all_pair_paths = DataFrame(pair =String[], files = Array[])
    unq_pairs = unique(Iterators.flatten([DataFrame(CSV.File(file)).Files for file in large_index]))
    for pair in unq_pairs
        pair_paths = Array{String, 1}(undef, 0)
        for csv in large_index
            try
                df = DataFrame(CSV.File(csv))
                path = df[(findall(x -> x==pair, df.Files)),:].paths[1]
                push!(pair_paths, path)
            catch 
                # That CSV doesn't have that pair - not really a problem!
            end
        end
        #print(pair_paths)
        push!(all_pair_paths, [pair, pair_paths])
    end
    return all_pair_paths
end

pair_paths_df = csv_merge(month_fnames)

@everywhere begin 
    using DataFrames, JLD2
    # get lat, lon, and elevation (LLE) for station 
    function LLE(station, df)
        row = df[(findall(x -> x==station, df.station)),:]
        lat, lon, el = row.latitude[1], row.longitude[1], row.elevation[1]
        return lat, lon, el
    end
    # Add location, distance and azi to a corr
    function add_corr_locations(corr::CorrData, df::DataFrame)
        # Get names
        source_station, reciever_station= split(corr.name, ".")[2], split(corr.name, ".")[end-2]
        # Get lat lons of each station
        lat1, lon1, el1 = LLE(source_station, df)
        lat2, lon2, el2 = LLE(reciever_station, df)
        # Convert to geoloc object 
        source_loc = GeoLoc(lat=lat1, lon=lon1, el=float(el1))
        reciever_loc = GeoLoc(lat=lat2, lon=lon2, el=float(el2))
        # compute necessary params and add to corr. 
        dist, azi, baz = get_dist(source_loc, reciever_loc), get_azi(source_loc, reciever_loc), get_baz(source_loc, reciever_loc)
        corr.loc = source_loc
        corr.dist, corr.azi, corr.baz = dist, azi, baz 
    end
    # gets array of corrs for particular file
    function corr_load(corr_large, key)
        jld = jldopen(corr_large, "r")
        files = keys(jld[key])
        corrs_comp = [jld["$key/$file"] for file in files]
        close(jld)
        return corrs_comp
    end
    function save_hdf5(C::CorrData, name)
        S=SeisData(1)
        fill!(S.id,C.id)
        fill!(S.name,C.name)
        fill!(S.loc,C.loc)
        fill!(S.fs,C.fs)
        fill!(S.misc,Dict([ ("corr_type",C.corr_type), ("cc_len",C.cc_len),("cc_step",C.cc_len),
        ("whitened",C.whitened), ("time_norm",C.time_norm), ("notes",C.notes),("dist",C.dist),
        ("azi",C.azi),("baz",C.baz),("maxlag",C.maxlag)]))
        # here i should convert the *t* as a date into a unix time.
        S.t[1] = [1 convert(Int,C.t[1] * 1e6); size(C.corr,1) 0]
        fill!(S.x,C.corr[:])
        write_hdf5(name,S)
    end
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
        filename = joinpath(CORROUT,"$(name).h5")
        save_hdf5(C, filename)
    end

    components =["EE", "EN", "EZ", "NE", "NN", "NZ", "ZE", "ZN", "ZZ"]
    # load, sum, stack and rotate correlations
    function postprocess_corrs(ar_corr_large, df)
        comp_9 = Array{CorrData,1}(undef, 0)
        # iterate and stack across components
        for key in components
            # Get all correlations for particular component
            corr_singles = Iterators.flatten([corr_load(corr_large, key) for corr_large in ar_corr_large])
            # stack and append to array
            corr_summed = SeisNoise.stack(sum(shorten.(corr_singles, 300.)), allstack=true)
            push!(comp_9, corr_summed)
        end
        # Location patch - only add to 6 stacked corrs
        map(x -> add_corr_locations(x, df), comp_9)
        SeisNoise.rotate!(comp_9, comp_9[1].azi, comp_9[1].baz) # rotate
        # write to disk
        map(C-> save_processed(C, "processed/$job_name"), comp_9)
        return comp_9 # to check work - probably don't need to return unless plotting
    end
end

# process raw correlations 
post_corrs = pmap(x -> postprocess_corrs(x, df), pair_paths_df.files)

#Transfer stacked correlations
stacked_files = joinpath.("processed/$job_name", readdir("processed/$job_name"))
@elapsed pmap(x -> s3_put(aws, bucket2, x, read(x)), stacked_files)
println("Done")





# Plotting functionality not yet available - coming soon 
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