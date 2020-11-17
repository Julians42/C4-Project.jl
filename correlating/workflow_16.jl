T = @elapsed using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays, AWSS3

using Pkg 
ENV["GR"] = ""
Pkg.build("GR")

#Add procs to access multiple cores
addprocs()

# Read in station locations and list source stations
df = CSV.File("files/all_locations_socal.csv") |> DataFrame! 
sources = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV"
        ,"VINE", "ROPE", "ARNO", "LUCI", "ROUF", "KUZD", "ALLI", "CHN", "USB", "Q0048"]

@everywhere using SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWSCore, StructArrays, AWSS3, Statistics, JLD2, Glob 
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
    function load_corrs(file_list::Array{String,1})
        corrs_in_pair = Array{CorrData,1}(undef, length(file_list))
        for (index, name) in enumerate(file_list)
            comp = string(split(name, "/")[end-1]) # get component
            corrs_in_pair[index] = load_corr(name, comp)
        end
        return corrs_in_pair
    end
    function write_jld2(corr_folder::String, file_dir::String)
        pair_corr_names = glob("$corr_folder/*/*.jld2","CORR")
        pair_corrs = load_corrs(pair_corr_names)
        if !isdir("$file_dir") # ensure filepathing
            mkpath("$file_dir")
        end
        fo = jldopen("$(file_dir)/$(corr_folder).jld2", "w") # name of the station pairs, most computationally expensive
        for (index, value) in enumerate(pair_corrs)
            # continue if no data in corr 
            isempty(value.corr) && continue
            #pair = join(deleteat!(split(corrs[index].name,"."),[3,4,7,8]),".")
            comp = value.comp
            starttime = string(Date(u2d(value.t[1]))) # eg 2017-01-01
            groupname = joinpath(comp, starttime)
            #!haskey(fo, pair) && JLD2.Group(fo, pair) # if it doesn't have pair, make pair layer
            !haskey(fo, comp) && JLD2.Group(fo, comp) # component layer
            !haskey(fo, groupname) && (fo[groupname] = value) # save corr into layer 
        end
        close(fo)
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
        for j in 1:length(indices[2]) # add distance checking - eg don't add anything over 300 km
            push!(pairs, [indices[1][i],indices[2][j]])
        end
    end
    return pairs
end

#coeffs - send to all cores
@everywhere begin 
    cc_step, cc_len = 3600, 3600
    maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
    freqmin, freqmax = 0.05, 9.9
    half_win, water_level = 30, 0.01
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    network = "CI"
    channel = "BH?"
    OUTDIR = "~/data"
end

# Select for months
#2017 minus jan 
#dates = [["2017-02-01","2017-02-28"],["2017-03-01","2017-03-31"],["2017-04-01","2017-04-30"],["2017-05-01","2017-05-31"],["2017-06-01","2017-06-30"],["2017-07-01","2017-07-31"],["2017-08-01","2017-08-31"],["2017-09-01","2017-09-30"],["2017-10-01","2017-10-31"],["2017-11-01","2017-11-30"],["2017-12-01","2017-12-31"]]
#2018 minus feb
#dates = [["2019-01-01","2019-01-31"],["2019-02-01","2019-02-28"],["2019-03-01","2019-03-31"],["2019-04-01","2019-04-30"],["2019-05-01","2019-05-31"],["2019-06-01","2019-06-30"]]#,["2018-07-01","2018-07-31"],["2018-08-01","2018-08-31"],["2018-09-01","2018-09-30"],["2018-10-01","2018-10-31"],["2018-11-01","2018-11-30"],["2018-12-01","2018-12-31"]]
#dates = [["2019-07-01","2019-07-31"],["2019-08-01","2019-08-31"],["2019-09-01","2019-09-30"],["2019-10-01","2019-10-31"],["2019-11-01","2019-11-30"],["2019-12-01","2019-12-31"]]
dates = [["2019-05-01","2019-05-31"], ["2019-11-01","2019-11-30"],["2019-12-01","2019-12-31"]]
for mth in dates
    startdate, enddate = mth[1], mth[2]
    days = Date(startdate):Day(1):Date(enddate)
    @eval @everywhere startdate, enddate = $startdate, $enddate
    @eval @everywhere days = $days
    num_corrs = 0
    T_start = Dates.now()
    for i in 1:length(days)
        try
            yr = Dates.year(days[i])
            path = join([Dates.year(days[i]),lpad(Dates.dayofyear(days[i]),3,"0")],"_") # Yeilds "YEAR_JDY"

            ############################ Data Download ###################################
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

            ####################### Read and Preprocess ####################################
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

            ###################### Index Pairs and Correlate ###############################
            
            indices = index_sources(clean_files, sources) # returns indices of sources
            sta_pairs = station_pairs(indices)

            # correlate
            T2 = @elapsed pmap(x -> correlate_pair(x, maxlag), map(y -> ffts[y], sta_pairs))
 
            println("Data for $(days[i]) correlated in $(T2) seconds!")

            # Perform cleanup of instance
            rm("data/continuous_waveforms", recursive=true) # Remove raw data to prevent memory crash 
        catch e
            println("Difficulty processing $(days[i]). Unexpected error. Continuing to next day.")
        end
    end

    # combine single day data to month files by station pair
    @eval @everywhere station_pair_names, month_, yr = readdir("/home/ubuntu/CORR"), Dates.monthname(Date(startdate)),Dates.year(Date(startdate)) 
    # month = Dates.monthname(Date(startdate)) # get the name of the month
    # yr = Dates.year(Date(startdate))

    if !isdir("home/ubuntu/corr_large/$month_")
        mkpath("home/ubuntu/corr_large/$month_")
    end
    if !isdir("month_index/$yr")
        mkpath("month_index/$yr")
    end

    jld_time = @elapsed pmap(x -> write_jld2(x, "corr_large/$yr/$month_"), station_pair_names) # combine corrs by station pair and write to single file 

    # as station pair names are filenames, we save filenames in a csv to read back during post-process 
    df = DataFrame(Files = station_pair_names, paths =[string("corr_large/$yr/$month_/",elt,".jld2") for elt in station_pair_names])
    CSV.write("month_index/$yr/$(month_).csv",df)

    ################### Transfer to S3 ##########################

    s3_put(aws, "seisbasin", "month_index/$yr/$(month_).csv", read("month_index/$yr/$(month_).csv"))

    month_files = joinpath.("corr_large/$yr/$month_", readdir("corr_large/$yr/$month_"))
    Transfer = @elapsed pmap(x ->s3_put(aws, "seisbasin", x, read(x)), month_files)
    println("$(length(month_files)) correlation files transfered to $(bucket2) in $Transfer seconds!") 


    println("$num_corrs total correlations processed")
    T_end = Dates.now()
    println(T_end-T_start)
    ############# Clean Up #################
    rm("CORR", recursive=true) # Remove single correlation data 
    rm("corr_large/$yr/$month_", recursive=true) #remove large correlation data
end









############### Design goal is below ##################
# currently struggling to implement because of differences between mapping within functions
# and outside of them 

################## Run Jobs ####################
Jan_2017 = @elapsed SeisCore("2017-01-01","2017-01-02")
# Jan_2017 = @elapsed SeisCore("2017-01-01","2017-01-31") # January 2017
# Feb_2017 = @elapsed SeisCore("2017-02-01","2017-02-28") # February 2017
# Mar_2017 = @elapsed SeisCore("2017-03-31","2017-03-31") # March 2017
################ Month Correlate Function ##################### --- DOESN'T WORK
function SeisCore(startdate::String, enddate::String)
    @everywhere startdate = $startdate
    @everywhere enddate = $enddate
    days = Date(startdate):Day(1):Date(enddate)
    @everywhere days = $days
    num_corrs = 0
    for i in 1:length(days) # loop through days 
        try # to catch any uncaught errors - will loose day if these come up (found ~2-3/year)
            yr = Dates.year(days[i])
            path = join([Dates.year(days[i]),lpad(Dates.dayofyear(days[i]),3,"0")],"_") # Yeilds "YEAR_JDY"

            ############################ Data Download ###################################
            filelist_scedc = s3query(aws, days[i], enddate = days[i], network=network, channel=channel)
            @everywhere filelist_scedc=$filelist_scedc
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
            @everywhere filelist_basin=$filelist_basin

            try
                ec2download(aws, bucket2, filelist_basin, OUTDIR)
                data_avail = true
            catch
                println("Error Downloading seisbasin data for $path. Potentially no data available.")
            end
            if data_avail == false # If there isn't data for that day from niether SCEDC or Seisbasin proceed to next day
                continue 
            end
            println("Data Download OK") ##########
            ####################### Read and Preprocess ####################################
            fpaths = readdir("/home/ubuntu/data/continuous_waveforms/$(Dates.year(days[i]))/$(path)")
            files = joinpath.("/home/ubuntu/data/continuous_waveforms/$(Dates.year(days[i]))/$(path)",fpaths)
            print
            T_load = @elapsed ar_file = pmap(x->load_file(x), files) # load data 
            ar = [elt[1] for elt in ar_file if isnothing(elt)==false]
            add_locations(ar, df) # add source locations
            clean_files = [elt[2] for elt in ar_file if isnothing(elt)==false]
            println("Data for $(days[i]) loaded in $T_load seconds. $(length(ar)) channels to be correlated, with $(length(ar_file)-length(ar)) channels discarded.")
            println("Load OK")
            T_preprocess = @elapsed fft_raw = pmap(x -> preprocess(x, fs, freqmin, freqmax, cc_step, cc_len, half_win, water_level), ar)
            println("Preprocess Ran")
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

            ###################### Index Pairs and Correlate ###############################
            
            indices = index_sources(clean_files, sources) # returns indices of sources
            sta_pairs = station_pairs(indices)

            # correlate
            T2 = @elapsed pmap(x -> correlate_pair(x, maxlag), map(y -> ffts[y], sta_pairs))

            # push to bucket 
            println("Data for $(days[i]) correlated in $(T2) seconds!")

            rm("data/continuous_waveforms", recursive=true) # Remove raw data to prevent memory crash 
        catch e
            print(e)
            print("Failure to process correlations for $(days[i]) - unknown error. Continuing to next day")
            continue
        end
    end

    # combine single day data to month files by station pair
    @everywhere station_pair_names, month_, yr = readdir("/home/ubuntu/CORR"), Dates.monthname(Date(startdate)),Dates.year(Date(startdate)) 


    if !isdir("home/ubuntu/corr_large/$yr/$month_")
        mkpath("home/ubuntu/corr_large/$yr/$month_")
    end
    if !isdir("month_index/$yr")
        mkpath("month_index/$yr")
    end

    jld_time = @elapsed pmap(x -> write_jld2(x, "corr_large/$yr/$month_"), station_pair_names) # combine corrs by station pair and write to single file 

    # as station pair names are filenames, we save filenames in a csv to read back during post-process 
    df = DataFrame(Files = station_pair_names, paths =[string("corr_large/$yr/$month_/",elt,".jld2") for elt in station_pair_names])
    CSV.write("month_index/$yr/$(month_).csv",df)

    ################### Transfer to S3 ##########################

    s3_put(aws, "seisbasin", "month_index/$yr/$(month_).csv", read("month_index/$yr/$(month_).csv"))

    month_files = joinpath.("corr_large/$yr/$month_", readdir("corr_large/$month"))
    Transfer = @elapsed pmap(x ->s3_put(aws, "seisbasin", x, read(x)), month_files)
    println("$(length(month_files)) correlation files transfered to $(bucket2) in $Transfer seconds!") 

    println("$num_corrs total correlations processed")

    ############### Corr Cleanup ####################
    rm("CORR", recursive=true) # Remove single correlation data 
    rm("corr_large/$yr/$month_", recursive=true) #remove large correlation data
end

