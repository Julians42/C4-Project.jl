T = @elapsed using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays, AWSS3

using Pkg 
ENV["GR"] = ""
Pkg.build("GR")

#Add procs to access multiple cores
addprocs()

#coeffs - send to all cores
@everywhere begin 
    cc_step, cc_len = 3600, 3600
    maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
    freqmin, freqmax = 0.05, 9.9
    half_win, water_level = 30, 0.01
end

@everywhere using SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWSCore, StructArrays, AWSS3, Statistics, JLD2
@everywhere begin 
    function is_window(C::SeisData, cc_len::Int64)
        """ Returns true if data has large enough ungapped window for correlation """
        windows = u2d.(SeisIO.t_win(C.t[1], C.fs[1]) .* 1e-6)
        startend = windows[:,2] .- windows[:,1] .> Second(cc_len)
        bool = any(startend)
        return bool
    end
    function load_file(file::String)
        """ Load raw data file """
        data = read_data("mseed", file)
        gaps = size(data.t[1])[1] # N-2 gaps (eg gaps = 12 tests for 10 gaps)
        pts = size(data.x[1])[1]
        fs_temp = data.fs[1]
        if gaps < 12 && is_window(data, cc_len) == true # If few gaps and sufficient (2/3) data present, data is clean
            return [data, file]
        end
    end
    function preprocess(S::SeisData, fs::Float64, freqmin::Float64, freqmax::Float64, cc_step::Int64, cc_len::Int64, half_win::Int64, water_level::Float64)
        """
            Pre-process raw seismic data object.
            - Removes mean from `S`.
            - Detrends each channel in `S`.
            - Tapers `S`
            - We recommend including `bandpass!` and `coherence` to improve signal
            - Downsamples data to sampling rate `fs`
            - Phase-shifts data to begin at 00:00:00.0
        """
        try
            process_raw!(S, fs)
            R = RawData(S,cc_len,cc_step)
            SeisNoise.detrend!(R)
            SeisNoise.taper!(R)
            bandpass!(R,freqmin,freqmax,zerophase=true)
            FFT = compute_fft(R) # Compute Fourier Transform
            coherence!(FFT,half_win, water_level)
            bool = true
            return [FFT, bool] 
        catch  e # Handle error catching if earlier filters haven't caught it yet
            println(e)
            bool = false
            return [nothing, bool]
        end
    end

    function cc_medianmute(A::AbstractArray, cc_medianmute_α::Float64 = 10.0)
        """
            Remove noisy correlation windows before stacking
            - Remove if average noise is greater than 10x the average
        """
        T, N = size(A)
        cc_maxamp = vec(maximum(abs.(A), dims=1))
        cc_medianmax = median(cc_maxamp)
        inds = findall(x-> x <= cc_medianmute_α*cc_medianmax,cc_maxamp)
        return A[:, inds], inds
    end
    remove_medianmute(C::CorrData, inds) = (return C.t[inds])
    function cc_medianmute!(C::CorrData, cc_medianmute_α::Float64 = 10.0)
        C.corr, inds = cc_medianmute(C.corr, cc_medianmute_α)
        C.t = remove_medianmute(C, inds)
        return nothing
    end
    function name_corr(C::CorrData)
        """ Returns corr name string: CH1.STA1.LOC1.CH2.STA2.LOC2 """
        return strip(join(deleteat!(split(C.name,"."),[4,8]),"."),'.')
    end
    function correlate_pair(ffts::Array{FFTData,1}, maxlag::Float64)
        """ 
            Correlation function for pair of fft data
            - noise filter and stacking
            - saves directly to disk: ensure directory is correct if not on AWS
        """
        C = correlate(ffts[1],ffts[2],maxlag)
        cc_medianmute!(C, 10.) # remove correlation windows with high noise
        stack!(C)
        pair, comp = name_corr(C), C.comp
        save_named_corr(C,"CORR/$pair/$comp")
    end
    function save_named_corr(C::CorrData, CORROUT::String)
        """ Implements custom naming scheme for project """
        CORROUT = expanduser(CORROUT) # ensure file directory exists
        if isdir(CORROUT) == false
            mkpath(CORROUT)
        end

        yr,j_day = Dates.year(Date(C.id)), lpad(Dates.dayofyear(Date(C.id)),3,"0") # YEAR, JULIAN_DAY
        p_name = name_corr(C) 
        name = join([yr,j_day,p_name],"_") #YEAR_JDY_CORRNAME

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
end


#Add location from dataframe to array 
function add_locations(ar::Array{SeisData,1},df::DataFrame)
    for i in 1:length(ar)
        try
            station_name = split(ar[i].name[1],".")[2]
            for j in 1:size(df)[1]
                if station_name == df[j,2]
                    ar[i].loc[1].lat = df[j,3]
                    ar[i].loc[1].lon = df[j,4]
                    ar[i].loc[1].el = df[j,5]
                    break
                end
            end
        catch
            println("Cannot add location to station $i.")
        end
    end
end

#Returns indices of source stations at index 1 and non-source stations at index 2
function index_sources(ar::Array{String,1}, stations::Array{String,1})
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

function station_pairs(indices::Array{Array{Int64,1}}) # returns array of station pairs to correlate
    pairs = Array{Array{Int64,1}}(undef, 0) #array of station pairs to correlate
    for i in 1:length(indices[1])
        for j in 1:length(indices[2])
            push!(pairs, [indices[1][i],indices[2][j]])
        end
    end
    return pairs
end

function load_lin(filelist::Array{String,1})
    ar_files = Array{SeisData,1}(undef,0)
    for i in 1:length(filelist)
        data = read_data("mseed", filelist[i])
        gaps = size(data.t[1])[1] # N-2 gaps (eg gaps = 12 tests for 10 gaps)
        pts = size(data.x[1])[1]
        fs_temp = data.fs[1]
        if gaps < 12 && pts > fs_temp*24*3600 -1 # If few gaps and sufficient (2/3) data present, data is clean
            push!(ar_files, [data, filelist[i]])
        end
    end
    return(ar_files)
end


function make_csv(ar::Array{SeisData, 1}, day_index::Int64, path::String)
    sta_index = DataFrame(ms_filename = String[],network = String[],station = String[], location = String[], channel = String[],lat = String[], lon=String[], el = String[], sample_rate=String[])
    for j in 1:length(ar)
        f_name = "$(ar[j].name).ms"
        sta_name = split(ar[j].name[1],'.')
        push!(sta_index, Dict(:ms_filename => "$f_name", :network => "$(sta_name[1])", 
                            :station => "$(sta_name[2])", :location=> "$(sta_name[3])",
                            :channel => "$(sta_name[4])", :lat => "$(ar[j].loc[1].lat)", 
                            :lon => "$(ar[j].loc[1].lon)", :el => "$(ar[j].loc[1].el)", :sample_rate => "$fs"))
    end
    CORROUT = "files/index/$(Dates.year(days[day_index]))" # ensure file directory exists
    if isdir(CORROUT) == false
        mkpath(CORROUT)
    end
    CSV.write("files/index/$(Dates.year(days[day_index]))/$(path)_correlation_index.csv", sta_index)
end

# Read in station locations and list source stations
df = CSV.File("files/all_locations_socal.csv") |> DataFrame! 
#sources = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV"]
sources = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV"
        ,"VINE", "ROPE", "ARNO", "LUCI", "ROUF", "KUZD", "ALLI"]

@everywhere begin
    using SCEDC, AWSCore, Dates 
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    startdate = "2017-01-17" # Select Start Date
    enddate = "2017-03-14" # Select End Date
    days = Date(startdate):Day(1):Date(enddate)
    network = "CI"
    channel = "BH?"
    OUTDIR = "~/data"
end
num_corrs = 0
for i in 1:length(days)
    yr = Dates.year(days[i])
    path = join([Dates.year(days[i]),lpad(Dates.dayofyear(days[i]),3,"0")],"_") # Yeilds "YEAR_JDY"
    # load with ec2stream ( preprocess within)
    filelist_scedc = s3query(aws, days[i], enddate = days[i], network=network, channel=channel)
    @eval @everywhere filelist_scedc=$filelist_scedc
    # Dowload Data and read file paths - replace with ec2stream when available
    data_avail = false
    try
        ec2download(aws, bucket, filelist_scedc, OUTDIR)
        data_avail = true
    catch
        println("Error Downloading SCEDC data for $path. Potentially no data available.")
    end

    #S3 query call for data file information
    dict = collect(s3_list_objects(aws, "seisbasin", "continuous_waveforms/$(yr)/$(path)/", max_items=1000))
    filelist_basin = Array{String,1}(undef,length(dict))
    #Index to filepath given by the "Key" element of the dictionary
    for i in 1:length(dict)
        filelist_basin[i] = dict[i]["Key"]
    end
    @eval @everywhere filelist_basin=$filelist_basin

    try
        ec2download(aws, bucket2, filelist_basin, OUTDIR)
        data_avail = true
    catch
        println("Error Downloading seisbasin data for $path. Potentially no data available.")
    end
    if data_avail == false # If there isn't data for that day from niether SCEDC or Seisbasin proceed to next day
        continue 
    end

    # read in data
    fpaths = readdir("/home/ubuntu/data/continuous_waveforms/$(Dates.year(days[i]))/$(path)")
    files = joinpath.("/home/ubuntu/data/continuous_waveforms/$(Dates.year(days[i]))/$(path)",fpaths)
    
    T_load = @elapsed ar_file = pmap(x->load_file(x), files) # load data 
    ar = [elt[1] for elt in ar_file if isnothing(elt)==false]
    add_locations(ar, df) # add source locations
    clean_files = [elt[2] for elt in ar_file if isnothing(elt)==false]
    println("Data for $(days[i]) loaded in $T_load seconds. $(length(ar)) channels to be correlated, with $(length(ar_file)-length(ar)) channels discarded.")

    T_preprocess = @elapsed fft_raw = pmap(x -> preprocess(x, fs, freqmin, freqmax, cc_step, cc_len, half_win, water_level), ar)

    ffts = [fft[1] for fft in fft_raw] # extract ffts from raw array
    bools = [fft[2] for fft in fft_raw] # extract processing success/failure bools from raw array
    println("Data for $(days[i]) preprocessed in $T_preprocess seconds.")
    ffts, ar, clean_files = ffts[bools], ar[bools], clean_files[bools]
    ar = ar[bools]
    clean_files = clean_files[bools]
    if any(bools == 0) == true # report if some channels failed in preprocessing
        num_bad_preprocess = length([bool for bool in bools if bool ==1])
        println("$num_bad_preprocess channels failed in preprocessing. $(length(ar)) channels to be correlated.")
    end
    
    indices = index_sources(clean_files, sources) # returns indices of sources
    sta_pairs = station_pairs(indices)

    # correlate
    T2 = @elapsed pmap(x -> correlate_pair(x, maxlag), map(y -> ffts[y], sta_pairs))

    # push to bucket 
    println("Data for $(days[i]) correlated in $(T2) seconds!")

    # transfer that day of data 
    make_csv(ar, i, path) # make csv to store station info of stations correlated 
    corr_paths = glob("*/*/$(path)*.jld2","CORR")
    # Transfer CSV and Data
    full_acl = "bucket-owner-full-control"
    s3_put(aws, "seisbasin", "corr_index/$(Dates.year(days[i]))/$(path)_correlation_index.csv",read("files/index/$(Dates.year(days[i]))/$(path)_correlation_index.csv"))
    Transfer = @elapsed pmap(x ->s3_put(aws, "seisbasin", "corr_data/$(join(deleteat!(split(x, "/"),[1]),"/"))", read(x)), corr_paths)
    println("$(length(corr_paths)) correlations transfered to $(bucket2) in $Transfer seconds!")
    num_corrs += length(corr_paths)
    # Perform cleanup of instance
    rm("CORR", recursive=true) # Remove transfered correlation data 
    rm("data/continuous_waveforms", recursive=true) # Remove raw data
end