T = @elapsed using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays, AWSS3

#coeffs
cc_step, cc_len = 3600, 3600
maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9
half_win, water_level = 30, 0.01

source_stations = ["IPT"] # Determine which stations to download

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

@everywhere begin
    using SCEDC, AWSCore, Dates 
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    startdate = "2017-01-17" # Select Start Date
    enddate = "2017-04-17" # Select End Date
    days = Date(startdate):Day(1):Date(enddate)
    network = "CI"
    channel = "BH?"
    OUTDIR = "~/data"
end


# Get list of CSVs to transfer corr_index/2017/2017_026_correlation_index.csv
list_csv = Array{String,1}(undef, length(days));
for (index, dy) in enumerate(days)
    year = Dates.year(dy)
    path = join([Dates.year(dy),lpad(Dates.dayofyear(dy),3,"0")],"_") # Yeilds "YEAR_JDY"
    list_csv[index] = "corr_index/$year/$(path)_correlation_index.csv"
end

#Download CSVs
@eval @everywhere list_csv = $list_csv # Retain SCEDC download functionality, parallel download not actually needed
ec2download(aws, bucket2, list_csv, "~/") # Retain seisbasin filepathing

# Write Filelist 
source_station = ["IPT"]
t1 = Dates.now()
filelist = nothing
for (index, dy) in enumerate(days)
    year = Dates.year(dy)
    path = join([Dates.year(dy),lpad(Dates.dayofyear(dy),3,"0")],"_") # Yeilds "YEAR_JDY"
    next_filelist = list_corrs("/home/ubuntu/corr_index/$year/$(path)_correlation_index.csv",source_station, path) 
    if index == 1
        filelist =next_filelist
    else
        filelist = vcat(filelist, next_filelist)
    end
end
println(length(filelist))

@everywhere begin
    function get_corr(file::String)
        # Transfers correlation with error catching
        try
            s3_get_file(aws, bucket2, file, file)
        catch 
            println("Could not transfer file $file")
        end
    end
end
function make_dir(file::String) 
    # This one doesn't like parallelization for some reason - "cannot serialize a running Task"
    OUTDIR = dirname(file)
    if !isdir(OUTDIR)
        mkpath(OUTDIR)
    end
end

# make directories for corr objects and transfer files
map(x -> make_dir(x), filelist)
T = @elapsed map(x -> get_corr(x), filelist) # This works under map but not pmap!
println(T)
   

correlations = joinpath.("/home/ubuntu/corrs", filelist)


function load_crap2(filelist::Array{String,1}) #Load corr filelist into array with error handling
    corrs = Array{CorrData,1}(undef,0)
    for (index, corr) in enumerate(filelist)
        try 
            elt = load_corr(corr, convert(String,split(corr,"/")[end -1]))
            push!(corrs, elt)
        catch 
            deleteat!(filelist, index)
        end
    end
    return corrs
end
corrs = load_crap2(correlations)

println(length(filelist))

@everywhere begin
    function load_filelist(file::String)
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
end

corrs = pmap(x -> load_filelist(x), filelist)
bools = [corr[2] for corr in corrs]
corr_data = [corr[1] for corr in corrs]
#join(split(filelist[i],"/")[1:3],"/")

# split elmts in filelist and return unique corr_pair and comps for stacking
truncated = Array{String,1}(undef,0)
for i in 1:length(filelist)
    elt = join(split(filelist[i],"/")[1:3],"/")
    push!(truncated, elt)
end
unq_path = unique(truncated)



corrs_stacked = pmap(x -> load_sum_stack(x, 300.), unq_path) #returns array of corr_data stacked over days

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


# example terminal transfer in ec2 to s3:    aws s3 cp /home/ubuntu/processed/[\"IPT\"]/ s3://seisbasin/processed/3month2017/IPT/ --recursive



scp -i /Users/julianschmitt/Downloads/my_aws_key.pem -r ubuntu@ec2-54-187-31-190.us-west-2.compute.amazonaws.com:/home/ubuntu/processed/CJM /Users/julianschmitt/Downloads/corr_data/
# transfer
aws s3 cp s3://seisbasin/corr_data /home/ubuntu/corr_data --recursive  --exclude "*"  --include "CI.IPT*.jld2"