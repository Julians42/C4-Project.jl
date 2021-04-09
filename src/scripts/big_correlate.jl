include("../SeisCore.jl")

# Include required packages
using .SeisCore, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWS, Distributed, JLD2, Glob, AWSS3

addprocs()
# add Constants and locations dataframe
@everywhere begin 
    using .SeisCore, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWS, Distributed, JLD2, Glob, AWSS3
    aws = AWS.AWSConfig(region="us-west-2")
    rootdir = ""
    network = "CI"
    channel1 = "BH?"
    channel2 = "HH?"
    OUTDIR = "~/data"
    cc_step, cc_len = 3600, 3600
    maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
    freqmin, freqmax = 0.05, 9.9
    half_win, water_level = 30, 0.01
    locations = DataFrame(CSV.File("root/locations.csv"))
end

# get date range from environment variable
arg = ENV["AWS_BATCH_JOB_ARRAY_INDEX"]
startdate = Date(2004)+Month(arg)
enddate = startdate+Month(1)-Day(1)
days = startdate:Day(1):enddate
println("Processing download for: ", startdate, " to ",enddate)

# go big or go home
function correlate_big(dd::Date, all_stations::DataFrame=locations)
    """ Wrapper function for daily correlations"""
    yr = Dates.year(dd)
    path = join([Dates.year(dd),lpad(Dates.dayofyear(dd),3,"0")],"_") # Yeilds "YEAR_JDY"
    @eval @everywhere path = $path
    # get filepaths for source stations - would be faster if they're available
    scedc_files = nothing
    try # most filelists are stored on seisbasin
        s3_get_file(aws, "seisbasin", "scedc_path/$yr/$path.csv", "scedc_path/$yr/$path.csv")
        scedc_files = DataFrame(CSV.File("scedc_path/$yr/$path.csv")).Path
    catch e # in case file not found/present we use SCEDC lookup
        println(e)
        scedc_files = get_scedc_files(dd, aws)
    end
    filter!(x -> any(occursin.(sources, x)), scedc_files)

    # filepaths for iris, ncedc, and node data
    iris_path  = S3Path("s3://seisbasin/iris_waveforms/$yr/$path/", config=aws)
    ncedc_path = S3Path("s3://seisbasin/ncedc_waveforms/$yr/$path/", config=aws)
    node_path  = S3Path("s3://seisbasin/continuous_waveforms/$yr/$path/", config=aws)
    iris_query, ncedc_query, node_query = convert.(String, readdir(iris_path)), convert.(String, readdir(ncedc_path)), 
                                        convert.(String, readdir(node_path))
    filelist_basin = vcat(joinpath.("iris_waveforms/$yr/$path/", iris_query), joinpath.("ncedc_waveforms/$yr/$path/", ncedc_query),
                                        joinpath.("continuous_waveforms/$yr/$path/", node_query))
    println("There are $(length(filelist_basin)) node files available for $path.")

    # download scedc and seisbasin data
    ec2download(aws, "scedc-pds", scedc_files, OUTDIR)
    ec2download(aws, "seisbasin", filelist_basin, OUTDIR)
    println("Download complete!")

    # preprocess data, integrating node data
    allf = glob("data/*/$yr/$path/*")
    accelerometers = filter(x -> isfile(x), joinpath.("data/continuous_waveforms/$yr/$path/", node_query))
    broadbands = setdiff(allf, accelerometers) # get the rest of the files 

    T_b = @elapsed pmap(f -> preprocess2(f), broadbands[1:50])
    if length(accelerometers) != 0 # save time on pmap allocation overhead
        T_a = @elapsed pmap(f -> preprocess2(f, true), accelerometers) # integrate accelerometers
        println("Preprocessing Completed in $(T_a+T_b) seconds.")
    else
        println("Preprocessing Completed in $T_b seconds.")
    end

    # Run autocorrelations
    fft_list_100 = glob("ffts/$path/100/*.jld2")
    fft_100_stations = unique([join(split(elt, ".")[1:2],".") for elt in fft_list_100]) # get stations to iterate over

    function correlate_pair(src::FFTData, rec::FFTData, maxlag::Float64, pref::String="CORR")
        C = correlate(src, rec, maxlag)
        cc_medianmute!(C, 10.) # remove correlation windows with high noise
        stack!(C)
        save_named_corr(C,"$pref/$(name_corr(C))/$(C.comp)")
    end

    function autocorrelate(station::String, list::Array{String,1}=fft_list_100)
        """ Filter channels by station, load and correlate """
        channels = filter(x -> occursin(station, x), list)
        ffts = map(x -> load_fft(x, string(x[end-7:end-5])), channels)
        pairs = vec([collect(x) for x in Iterators.product(1:length(ffts),1:length(ffts))])
        map(pair -> correlate_pair(pair[1], pair[2], maxlag, "AUTOCORR"), pairs)
    end
    pmap(x -> autocorrelate(x, fft_list_100), fft_100_stations)

    # get indices for block correlation
    fft_paths = sort(glob("ffts/$path/20/*")) # alphabetic sort for all stations

    # Chunk raw ffts into blocks. 30 is best size for IO, but distribute across cores if more cores available
    chunks = collect(Iterators.partition(fft_paths, convert(Int64, minimum(30, ceil(length(fft_paths)/nprocs()))))) 

    # run diagonal chunks
    function diag_chunks(chunk::Array{String,1}, prefix::String = "CORR", filt_dist::Bool = true)
        ffts  = map(x -> load_fft(x, string(x[end-7:end-5])), chunk)
        pairs = vec([collect(x) for x in Iterators.product(1:length(ffts),1:length(ffts))])
        # maybe remove this filter so we get all components 
        #filter!(x -> (x[1] < x[2]) || (x[1]>x[2]), pairs) # filter upper diagonal
        if filt_dist; filter!(x -> get_dist(ffts[x[1]], ffts[x[2]]) <= 300., pairs); end # filter distances
        map(pair -> correlate_pair(ffts[pair[1]], ffts[pair[2]], maxlag, prefix), pairs)
    end
    pmap(chunk -> diag_chunks(convert(Array,chunk)), chunks)

    function offdiag_chunks(chunkpair::Array{Array{String, 1}}, prefix::String = "CORR", filt_dist::Bool = true)
        sources   = map(x -> load_fft(x, string(x[end-7:end-5])), chunkpair[1]) # load
        receivers = map(x -> load_fft(x, string(x[end-7:end-5])), chunkpair[2]) 
        pairs = vec([collect(x) for x in Iterators.product(1:length(sources),1:length(receivers))]) # all pairs src -> rec
        if filt_dist; filter!(x -> get_dist(sources[x[1]].loc, receivers[x[2]].loc) <= 300., pairs); end # filt distances
        map(pair -> correlate_pair(sources[pair[1]], receivers[pair[2]], maxlag, prefix), pairs)
    end
    off_chunks = vec([collect(x) for x in Iterators.product(1:length(chunks),1:length(chunks))])
    filter!(chunk_pair -> chunk_pair[1] < chunk_pair[2], off_chunks) # get upper non-diagonal chunks

    pmap(chunk -> offdiag_chunks(chunk), off_chunks) # run mega correlations

    rm("data/continuous_waveforms", recursive=true) # cleanup raw data
end

# guessing this will process ~ 60 million correlations monthly - hopefully enough space on instance...
pmap(x -> correlate_big(x), days)

# get autocorr names
autocorr_names = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("AUTOCORR/*")]

function stack_auto(name::String, startdate::Date=startdate)
    month = Dates.monthname(startdate) # get month for filename
    yr = Dates.year(startdate)
    # stack autocorrelations 

    CORROUT = expanduser("autocorrelations/$name/")
    if !isdir(CORROUT)
        mkpath(CORROUT)
    end
    filename = joinpath(CORROUT,"$(yr)_$(month)_$name.h5") # get output filename

    # Get list of files to save
    autocorr_list = glob("AUTOCORR/$name*/*/*")
    C = load_corr(autocorr_list[1], convert(String, split(autocorr_list[1],"/")[end-1]))
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
            write(file, "meta/lat", C.loc.lat)
            write(file, "meta/lon", C.loc.lon)
            write(file, "meta/el", C.loc.el)
            write(file, "meta/dep", C.loc.dep)
            write(file, "meta/az", C.loc.az)
            write(file, "meta/inc", C.loc.inc)
        end
        for (ind, name) in enumerate(autocorr_list)
            # load file 
            comp = string(split(name, "/")[end-1]) # get component
            Cl = load_corr(name, comp)
            # write 
            write(file, "$comp/$(Cl.id)", Cl.corr[:])
        end
    end
end
pmap(x -> stack_auto(x, startdate), autocorr_names)

corr_names = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("CORR_20HZ/*")]

function stack_corr(name::String, startdate::Date=startdate, prefix::String = "CORR_20HZ")
    month = Dates.monthname(startdate) # get month for filename
    yr = Dates.year(startdate)
    # stack autocorrelations 

    CORROUT = expanduser("correlations/$name/")
    if !isdir(CORROUT)
        mkpath(CORROUT)
    end
    filename = joinpath(CORROUT,"$(yr)_$(month)_$name.h5") # get output filename
    # r = r"(?<=$(name)).*"
    components = ["EE","EN","EZ", "NE", "NN","NZ", "ZE", "ZN", "ZZ"]

    receivers = glob("$prefix/$name/*/*")
    receivers = Set([join(split(x, ".")[end-2:end-1], ".") for x in receivers]) # just get reciever names
    println(receivers)
    corr_list = glob("$prefix/$name*/*/*")
    C = load_corr(corr_list[1], convert(String, split(corr_list[1],"/")[end-1]))
    source_loc = GeoLoc(lat = C.loc.lat, lon = C.loc.lon, el = C.loc.el)
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
            write(file, "meta/lat", C.loc.lat)
            write(file, "meta/lon", C.loc.lon)
            write(file, "meta/el", C.loc.el)
            write(file, "meta/dep", C.loc.dep)
            write(file, "meta/az", C.loc.az)
            write(file, "meta/inc", C.loc.inc)
        end
        for rec in receivers
            try
                # sample_r = glob("CORR/$tf*$rec*/ZZ/*")[1]
                # Cr = load_corr(sample_r, "ZZ")
                rec_station = split(rec, ".")[2]
                #rec_loc = LLE_geo(rec_station, all_stations)
                rec_loc = load_fft(glob("ffts/*/*/$rec*BHZ*")[1], "BHZ").loc
                if !haskey(read(file), "$rec/meta")
                    write(file, "$rec/meta/lon", rec_loc.lon)
                    write(file, "$rec/meta/lat", rec_loc.lat)
                    write(file, "$rec/meta/el", rec_loc.el)
                    write(file, "$rec/meta/dist", get_dist(source_loc, rec_loc))
                    write(file, "$rec/meta/azi", get_azi(source_loc, rec_loc))
                    write(file, "$rec/meta/baz", get_baz(source_loc, rec_loc))
                end
                for comp in components
                    # load correlations for this receiver by component 
                    files = glob("$prefix/$name/$comp/*$rec*")
                    corrs = [load_corr(f, comp) for f in files]
    
                    # implement various stacktypes
                    corr_mean = SeisNoise.stack(sum(corrs), allstack=true, stacktype=mean)
                    corr_pws = SeisNoise.stack(sum(corrs), allstack=true, stacktype=pws)
                    corr_robust = SeisNoise.stack(sum(corrs), allstack=true, stacktype=robuststack)
                    # save into file 
                    write(file, "$rec/$comp/linear", corr_mean.corr[:])
                    write(file, "$rec/$comp/pws", corr_pws.corr[:])
                    write(file, "$rec/$comp/robust", corr_robust.corr[:])
                end
            catch e
                println(e)
            end
        end
    end
end

pmap(x -> stack_corr(x, startdate), corr_names)


# transfer to s3 
autos = glob("autocorrelations/*/*")
bigcorrs = glob("correlations/*/*")

pmap(x -> s3_put(aws, "seisbasin", x, read(x)), autos)
pmap(x -> s3_put(aws, "seisbasin", x, read(x)), bigcorrs)


### IGNORE ####

function preprocess3(file::String,  accelerometer::Bool=false, rootdir::String="", samp_rates::Array{Float64, 1}=[1., 20., 100.],
    freqmin::Float64=freqmin, freqmax::Float64=freqmax, cc_step::Int64=cc_step, 
    cc_len::Int64=cc_len, half_win::Int64=half_win, water_level::Float64=water_level, 
    all_stations::DataFrame=all_stations, path::String=path)
    """
    Load raw seisdata file and process, saving to fft
    """
    try
        data = SeisData()
        if occursin(".ms", file)
            read_data!(data, "mseed", file)
        else # data is SeisIO native from NCEDC/IRIS
            read_data!(data, "seisio", file)
        end
        println(data)
        gaps = size(data.t[1])[1] # N-2 gaps (eg gaps = 12 tests for 10 gaps)
        pts = size(data.x[1])[1]
        fs_temp = data.fs[1]
        windows = u2d.(SeisIO.t_win(data.t[1], data.fs[1]) .* 1e-6)
        startend = windows[:,2] .- windows[:,1] .> Second(cc_len)
        if gaps < 25 && any(startend) # If few gaps and sufficient (2/3) data present, data is clean
            try
                add_location(data, all_stations)
                for samp_rate in samp_rates
                    S = SeisData()
                    try
                        S = process_raw(data, samp_rate)
                    catch # trying to sample non 100Hz data at 100 Hz - skip resampling 
                        S = process_raw(data, data.fs[1])
                    end
                    println(S)
                    R = RawData(S,cc_len,cc_step)
                    SeisNoise.detrend!(R)
                    SeisNoise.taper!(R)
                    bandpass!(R,freqmin,freqmax,zerophase=true)
                    FFT = nothing
                    if accelerometer # accelerometer - so we integrate
                        FFT = rfft_raw(R,1)
                    else # then its a seismometer
                        FFT = compute_fft(R)
                    end
                    println(FFT)
                    coherence!(FFT,half_win, water_level)
                    try # save fft 
                        root_fft = "ffts/$path/$(Int(samp_rate))/"
                        save_fft(FFT, joinpath(rootdir, root_fft))
                    catch e
                        println(e)
                    end
                    #println("Successfully processed $(data.id[1])")
                end
                return data.id[1]
            catch e
                println(e)
            end
        end
    catch e 
        println(e)
    end
end