T = @elapsed using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, 
                    JLD2, Statistics, PyCall, Glob, StructArrays, AWSS3, SeisIO.SEED

#Add procs to access multiple cores
addprocs()
@everywhere using SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWSCore, StructArrays, AWSS3, Statistics, JLD2, Glob, SeisIO.SEED

######################## Meta Data ###############################
#coeffs - send to all cores
@everywhere begin 
    cc_step, cc_len = 3600, 3600
    maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
    freqs = [20., 100.]
    freqmin, freqmax = 0.05, 9.9
    half_win, water_level = 30, 0.01
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    network = "CI"
    channel1 = "BH?"
    channel2 = "HH?"
    OUTDIR = "~/data"
end

# select start/ enddate (Default calculates for entire month: eg start_date on 003 rounds to 001)
# start_date, end_date = "2018-07-01", "2018-12-31" # depreciated in favor of ARGS
mem = true # choose whether to favor smaller I/O (true) or less inter-core communication (false)
chunk_length = 30

# Read in station locations and list source stations
all_stations = DataFrame(CSV.File("files/full_socal.csv"))
##############################################################################

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
        if gaps < 25 && is_window(data, cc_len) == true # If few gaps and sufficient (2/3) data present, data is clean
            return [data, file]
        end
    end
    function LLE_geo(station, df)
        """ Find station matching location and return geoloc object"""
        try
            row = df[(findfirst(x -> x==station, df.station)),:]
            lat, lon, el = row.latitude[1], row.longitude[1], row.elevation[1]
            geo = GeoLoc(lat = float(lat), lon = float(lon), el = float(el))
            return geo
        catch 
            return nothing
        end
    end    
    function add_location(seis::SeisData,df::DataFrame)
        """ Adds locations to array of seisdata from a dataframe """
        name = split(seis.id[1],".")[2]
        geo = LLE_geo(name, df)
        if !isnothing(geo)
            seis.loc[1] = geo
        else
            #println("Station $name doesn't have a location in the dataframe")
        end
    end
    function preprocess(S::SeisData, fs::Float64, freqmin::Float64, freqmax::Float64, cc_step::Int64, 
                        cc_len::Int64, half_win::Int64, water_level::Float64)
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
            if fs <= 20.
                process_raw!(S, fs)
            else
                process_raw!(S, S.fs[1])
            end
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
    function all_preprocess(file, samp_rates::Array{Float64,1}, freqmin::Float64, freqmax::Float64, cc_step::Int64, 
                            cc_len::Int64, half_win::Int64, water_level::Float64)
        """
            Load raw seisdata file and process, saving to fft per frequency
        """
        data = SeisData()
        if occursin("continuous", file[1])
            read_data!(data, "mseed", file[1])
        else
            read_data!(data, "seisio", file[1])
        end
        gaps = size(data.t[1])[1] # N-2 gaps (eg gaps = 12 tests for 10 gaps)
        pts = size(data.x[1])[1]
        fs_temp = data.fs[1]
        if gaps < 25 && is_window(data, cc_len) == true # If few gaps and sufficient (2/3) data present, data is clean
            try
                add_location(data, all_stations)
                # continue with preprocessing 
                for samp_rate in samp_rates
                    println(data.fs[1], samp_rate)
                    S = SeisData()
                    try
                        S = process_raw(data, samp_rate)
                    catch # trying to sample non 100Hz data at 100 Hz - skip resampling 
                        S = process_raw(data, data.fs[1])
                    end
                    R = RawData(S,cc_len,cc_step)
                    SeisNoise.detrend!(R)
                    SeisNoise.taper!(R)
                    bandpass!(R,freqmin,freqmax,zerophase=true)
                    FFT = compute_fft(R) # Compute Fourier Transform
                    coherence!(FFT,half_win, water_level)
                    try # save fft 
                        FFT.gain = file[2]
                        root_fft = "ffts/$path/$(Int(samp_rate))"
                        save_fft(FFT, root_fft)
                    catch e
                        println(e)
                    end
                end
                return data.id[1]
            catch e
                println(e)
            end
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
    function save_named_corr(C::CorrData, CORROUT::String)
        """ Implements custom naming scheme for project """
        CORROUT = expanduser(CORROUT) # ensure file directory exists
        if !isdir(CORROUT)
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
    function get_pairs(len_fft::Int64, len_split::Int64)
        """ Returns indicies for parallelized batch correlation job """
        # define partition
        partition = collect(Iterators.partition(collect(1:len_fft),len_split))
        # make indices
        jobs = Array{Array{Int64,1}}(undef,0)
        for p in 1:length(partition)
            for j in p:length(partition)
                push!(jobs, [p,j])
            end
        end
        fft_groups = [[partition[job[1]],partition[job[2]]] for job in jobs]
        return fft_groups
    end
    function correlate_load(fft_name_ar::Array{Array{String,1},1})
        # Load ffts for sources and receivers
        sources = map(x-> load_fft(x, string(x[end-7:end-5])), fft_name_ar[1])
        receivers = map(x-> load_fft(x, string(x[end-7:end-5])), fft_name_ar[2])

        # Correlate
        #all2all(sources, receivers)
        for s_fft in sources
            for r_fft in receivers
                if (s_fft.gain <= r_fft.gain) && (get_dist(s_fft.loc,r_fft.loc) <= 300.)
                    try
                        C = correlate(s_fft,r_fft, maxlag)
                        cc_medianmute!(C, 10.) # remove correlation windows with high noise
                        stack!(C)
                        pair, comp = name_corr(C), C.comp
                        if s_fft.loc == r_fft.loc # save autocorrs in different rootdir 
                            save_named_corr(C,"AUTOCORR/$pair/$comp")
                        else
                            save_named_corr(C,"CORR/$pair/$comp")
                        end
                    catch e
                        println(name_corr(correlate(s_fft, r_fft, maxlag)))
                    end
                end
            end
        end
    end
    function test_correlate_load(fft_name_ar::Array{Array{String,1},1}, maxdist)
        # Load ffts for sources and receivers
        sources = map(x-> load_fft(x, string(x[end-7:end-5])), fft_name_ar[1])
        receivers = map(x-> load_fft(x, string(x[end-7:end-5])), fft_name_ar[2])

        # Correlate
        indices = []
        for s_fft in sources
            for r_fft in receivers
                if (s_fft.gain <= r_fft.gain) && (get_dist(s_fft.loc,r_fft.loc) <= maxdist)
                    push!(indices, get_dist(s_fft.loc,r_fft.loc))
                elseif (s_fft.gain <= r_fft.gain)
                    push!(indices, nothing)
                end
            end
        end
        return indices
    end
end

function divide_months(t_start::String, t_end::String)
    """ Return all month windows """
    startdate, enddate = Date(t_start), Date(t_end)
    start_ym, end_ym = Dates.yearmonth(startdate), Dates.yearmonth(enddate)
    last_first, last_last = Dates.daysinmonth(startdate), Dates.daysinmonth(enddate)
    m_start1, m_start2 = Date("$(start_ym[1])-$(start_ym[2])-01"), Date("$(end_ym[1])-$(end_ym[2])-01") # start of first/last month
    m_end1, m_end2 = Date("$(start_ym[1])-$(start_ym[2])-$last_first"), Date("$(end_ym[1])-$(end_ym[2])-$last_last") # end of first/last month
    start_range = m_start1:Month(1):m_start2 # month start values
    end_range = m_end1:Month(1):m_end2 # month end values
    dates = [[start_range[ind], end_range[ind]] for ind in 1:length(start_range)] # start and end of each month
    return dates
end


#Returns indices of source stations at index 1 and non-source stations at index 2
function correlate_indices(to_correlate::Array{String, 1}, sources::Array{String,1})
    pairs = Array{Array{Int64,1}}(undef, 0) #array of station pairs to correlate
    source_indices = Array{Int64, 1}(undef, 0)
    # Find stations in successfully preprocessed ffts which are sources
    for (ind, station) in enumerate(to_correlate)
        loc = findfirst(occursin.(station[1:end-1], sources))
        if !isnothing(loc)
            push!(source_indices, ind)
        end
    end
    # get correlation pairs - all sources to all recievers. 
    for source_loc in source_indices
        for (ind, rec_loc) in enumerate(to_correlate)
            push!(pairs, [source_loc, ind])
        end
    end
    return pairs 
end

function get_dict_name(file::String)
    station = convert(String, split(split(file,"/")[end],"_")[1])
    component = split(file,"/")[end][10:12]
    return string(station, "_", component)
end

function foldersize(dir=".")
    """ returns total size of folder in GB """
    size = 0
    for (root, dirs, files) in walkdir(dir)
        size += sum(map(filesize, joinpath.(root, files)))
    end
    return size*10e-10
end

startdate, enddate = "2014-01-01", "2014-01-02"
job_id = Dates.now() # we choose a timestamp for unique job id - only important for job summary output

chunk_length = 30
mem = false # choose whether to favor smaller I/O (true) or less inter-core communication (false). Most efficient depends on system
single = false # mem must be set to false
days = Date(startdate):Day(1):Date(enddate)

i=1
yr = Dates.year(days[i])
path = join([Dates.year(days[i]),lpad(Dates.dayofyear(days[i]),3,"0")],"_") # Yeilds "YEAR_JDY"
@eval @everywhere path, all_stations = $path, $all_stations
# Download SCEDC and Seisbasin (nodes, NCEDC, IRIS) data 

# get SCEDC data
if !isdir("scedc_path/$yr/")
    mkpath("scedc_path/$yr/")
end
s3_get_file(aws, "seisbasin", "scedc_path/$yr/$path.csv", "scedc_path/$yr/$path.csv")
scedc_files = DataFrame(CSV.File("scedc_path/$yr/$path.csv")).Path

ec2download(aws, bucket, scedc_files, OUTDIR)


# get seisbasin data
filelist_basin = Array{String,1}(undef,0)
for data_source in ["continuous","iris","ncedc"]
    dict = collect(s3_list_objects(aws, "seisbasin", "$(data_source)_waveforms/$(yr)/$(path)/", max_items=1000))
    #Index to filepath given by the "Key" element of the dictionary
    for ind in 1:length(dict)
        push!(filelist_basin, dict[ind]["Key"])
    end
end
filter!(x -> x[end] != '/', filelist_basin)
@eval @everywhere filelist_basin=$filelist_basin

try
    ec2download(aws, bucket2, filelist_basin, OUTDIR)
    data_avail = true
catch
    println("Error Downloading seisbasin data for $path. Potentially no data available.")
end

################# All Preprocess Script ##################
files_seis = Glob.glob("data/*_waveforms/$(Dates.year(days[i]))/$(path)/*")
zip_gain = collect(zip(files_seis, collect(1:length(files_seis))))
Total_preprocess = pmap(file_gain -> all_preprocess(file_gain, freqs, freqmin, freqmax, cc_step, cc_len, half_win, water_level), zip_gain)

function list_name(s::String)
    name = convert(String, split(s,"/")[end])
    if occursin("continuous", s)
        net, cha, sta, comp = name[1:2], strip(convert(String, name[3:6]),'_'), 
                            strip(convert(String, name[11:12]),'_'), name[8:10]
        return join([net, cha, sta, comp],".")
    else
        return name
    end
end
seis_names = map(x -> list_name(x), files_seis)
fft_list_20 = ["ffts/$path/20/$name.jld2" for name in seis_names if isfile("ffts/$path/20/$name.jld2")]
fft_list_100 = ["ffts/$path/100/$name.jld2" for name in seis_names if isfile("ffts/$path/100/$name.jld2")]




# run Autocorrs 
#get station data
fft_100_stations = unique([join(split(elt, ".")[1:2],".") for elt in fft_list_100])

@everywhere begin
    function autocorrelate(station_name::String)
        station = map(x-> load_fft(x, string(x[end-7:end-5])), glob("$(station_name)*"))
        for (ind, src) in enumerate(station)
            for rec_ind in ind:length(station)
                C = correlate(src,station[rec_ind], maxlag)
                # Upsample C = upsample!(C)

                cc_medianmute!(C, 10.) # remove correlation windows with high noise
                stack!(C)
                pair, comp = name_corr(C), C.comp
                save_named_corr(C,"AUTOCORR/$pair/$comp")
            end
        end
    end
end
pmap(x-> autocorrelate(x), fft_100_stations)

# run all correlations
pairs = get_pairs(length(fft_list_20), chunk_length); # select chunk length 
fft_names = [[fft_list_20[pair[1]], fft_list_20[pair[2]]] for pair in pairs];


@everywhere begin
    function correlate_load(fft_name_ar::Array{Array{String,1},1})
        # Load ffts for sources and receivers
        sources = map(x-> load_fft(x, string(x[end-7:end-5])), fft_name_ar[1])
        T_load = @elapsed receivers = map(x-> load_fft(x, string(x[end-7:end-5])), fft_name_ar[2])
        println("Recievers loaded in $T_load seconds.")

        # Correlate
        #all2all(sources, receivers)
        for s_fft in sources
            for r_fft in receivers
                T_dist = @elapsed dist = get_dist(s_fft.loc,r_fft.loc) 
                println("Distance in $T_dist with distance $dist")
                if (s_fft.gain <= r_fft.gain) && (dist <= 300.) 
                    try
                        T_correlate = @elapsed C = correlate(s_fft,r_fft, maxlag)
                        T_mm = @elapsed cc_medianmute!(C, 10.) # remove correlation windows with high noise
                        T_stack = @elapsed stack!(C)
                        T_name = @elapsed pair, comp = name_corr(C), C.comp
                        T_save = @elapsed save_named_corr(C,"CORR/$pair/$comp")
                        println("Correlate: $T_correlate. MM: $T_mm. stack: $T_stack. name: $T_name. save:$T_save")
                    catch e
                        println(e)
                    end
                end
            end
        end
    end
end
T_load = @elapsed pmap(x-> correlate_load(x), fft_names)
#rm(root_fft, recursive=true) # cleanup 
println("Pairs correlated for $path from memory in $T_load seconds.")


# Postprocessing and stacking
# first process autocorrelations into h5 files
month = Dates.month(mth[1]) # get month
autocorrs = glob("AUTOCORR/*") # list all autocorrelation sources
autocorr_names = [join(split(auto, ".")[1:2],"."),autocorrs]

autocorrs[1] # gives first autocorrelation folder

name_source = "$(month)_$(autocorrs[i])"
CORROUT = expanduser("autocorrelations/$yr/$(autocorrs[i])/")
if !isdir(CORROUT)
    mkpath(CORROUT)
end
filename = joinpath(CORROUT,"$(name_source).h5") # get output filename

# Get list of files to save
autocorr_list = glob("AUTOCORR/$(autocorrs[i])/*/*.jld2")
C = corr_load(autocorr_list[1], convert(String, split(autocorr_list[1],"/")[1]))
h5open(filename, "cw") do file 
    if !haskey(read(file), "meta")
        write(file, "meta/corr_type", C.corr_type)
        write(file, "meta/cc_len", C.cc_len)
        write(file, "meta/cc_step", C.cc_step)
        write(file, "meta/whitened", C.whitened)
        write(file, "meta/time_norm", C.time_norm)
        write(file, "meta/notes", C.notes)
        write(file, "meta/maxlag", C.maxlag)
        write(file, "meta/starttime", Dates.format(u2d(C.t[1]), "yyyy-mm-dd HH:MM:SS"))
        write(file, "meta/samp_freq", C.fs)
        write(file, "meta/lat", C.loc[1].lat)
        write(file, "meta/lon", C.loc[1].lon)
        write(file, "meta/el", C.loc[1].el)
        write(file, "meta/dep", C.loc[1].dep)
        write(file, "meta/az", C.loc[1].az)
        write(file, "meta/inc", C.loc[1].inc)
    end
    for (ind, name) in enumerate(autocorr_list)
        # load file 
        comp = string(split(name, "/")[end-1]) # get component
        Cl = load_corr(name, comp)
        # write 
        write(file, "$comp/$(Cl.id)", Cl.corr[:])
    end
end

# need to write correlation stacking code ....



# transfer to S3
