include("../SeisCore.jl")

# Include required packages
using .SeisCore, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWS, Distributed, JLD2, Glob, AWSS3

addprocs()
# add Constants and locations dataframe
@everywhere begin 
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
    catch # in case file not found/present we use SCEDC lookup
        scedc_files = get_scedc_files(dd, aws)
    end
    filter!(x -> any(occursin.(sources, x)), scedc_files)

    # filepaths for iris, ncedc, and node data
    iris_path  = S3Path("s3://seisbasin/iris_waveforms/$yr/$path/", config=aws)
    ncedc_path = S3Path("s3://seisbasin/ncedc_waveforms/$yr/$path/", config=aws)
    node_path  = S3Path("s3://seisbasin/continuous_waveforms/$yr/$path/", config=aws)
    iris_query, ncedc_query, node_query = convert.(String, readdir(iris_path)), convert.(String, readdir(ncedc_path)), 
                                        convert.(String, readdir(node_path))
    filelist_basin = vcat(joinpath.("iris_waveforms/$yr/$path/", iris_path), joinpath.("ncedc_waveforms/$yr/$path/", ncedc_query),
                                        joinpath.("continuous_waveforms/$yr/$path/", node_query))
    println("There are $(length(filelist_basin)) node files available for $path.")

    # download scedc and seisbasin data
    ec2download(aws, scedc, scedc_files, OUTDIR)
    ec2download(aws, basin, filelist_basin, OUTDIR)
    println("Download complete!")

    # preprocess data, integrating node data
    allf = glob("data/*/$yr/$path/*")
    accelerometers = filter(x -> isfile(x), joinpath.("data/continuous_waveforms/$yr/$path/", node_query))
    broadbands = setdiff(allf, accelerometers) # get the rest of the files 

    T_b = @elapsed pmap(f -> preprocess(f), broadbands)
    T_a = @elapsed pmap(f -> preprocess(f, true), accelerometers) # integrate accelerometers
    println("Preprocessing Completed")

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
    function diag_chunks(chunk::Array{String,1})
        ffts = map(x -> load_fft(x, string(x[end-7:end-5])), chunk)
        pairs = vec([collect(x) for x in Iterators.product(1:length(ffts),1:length(ffts))])
        filter!(x -> (x[1] < x[2]), pairs) # filter upper diagonal
        filter!(x -> get_dist(ffts[x[1]], ffts[x[2]]) <= 300., pairs) # filter distances
        map(pair -> correlate_pair(pair[1], pair[2], maxlag), pairs)
    end
    pmap(chunk -> diag_chunks(chunk), chunks)

    function offdiag_chunks(chunkpair::Array{Array{String, 1}})
        sources   = map(x -> load_fft(x, string(x[end-7:end-5])), chunkpair[1]) # load
        receivers = map(x -> load_fft(x, string(x[end-7:end-5])), chunkpair[2]) 
        pairs = vec([collect(x) for x in Iterators.product(1:length(sources),1:length(receivers))]) # all pairs src -> rec
        filter!(x -> get_dist(sources[x[1]].loc, receivers[x[2]].loc) <= 300., pairs) # filt distances
        map(pair -> correlate_pair(pair[1], pair[2], maxlag), pairs)
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

function stack_auto(name::String)
    month = Dates.month(mth[1]) # get month for filename
    # stack autocorrelations 

    CORROUT = expanduser("autocorrelations/$yr/$name/")
    if !isdir(CORROUT)
        mkpath(CORROUT)
    end
    filename = joinpath(CORROUT,"$(month)_$name.h5") # get output filename

    # Get list of files to save
    autocorr_list = glob("AUTOCORR/$name*/*/*")
    C = corr_load(autocorr_list[1], convert(String, split(autocorr_list[1],"/")[end-1]))
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
end


