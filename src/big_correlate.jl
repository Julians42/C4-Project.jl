# Functions for large C4 correlation job 
export correlate_pair, autocorrelate, diag_chunks, offdiag_chunks, preprocess2, get_blocks,
        stack_auto, stack_corr, correlate_big, stack_all


# correlation helper functions
function correlate_pair(src::FFTData, rec::FFTData, pref::String="CORR", params::Dict=params)
    try
        C = correlate(src, rec, params["maxlag"])
        cc_medianmute!(C, 10.) # remove correlation windows with high noise
        stack!(C)
        name = join(split(name_corr(C), ".")[1:2],".")
        save_named_corr(C,"$pref/$name/$(C.comp)")
    catch e
        println(e)
    end
end
function autocorrelate(station::String, list::Array{String,1}=fft_list_100, params::Dict=params)
    """ Filter channels by station, load and correlate """
    channels = filter(x -> occursin(station, x), list)
    ffts = map(x -> load_fft(x, string(x[end-7:end-5])), channels)
    pairs = vec([collect(x) for x in Iterators.product(1:length(ffts),1:length(ffts))])
    map(pair -> correlate_pair(ffts[pair[1]], ffts[pair[2]], "AUTOCORR", params), pairs)
end
function diag_chunks(chunk::Array{String,1}, prefix::String = "CORR", filt_dist::Bool = true, params::Dict=params)
    ffts  = map(x -> load_fft(x, string(x[end-7:end-5])), chunk)
    pairs = vec([collect(x) for x in Iterators.product(1:length(ffts),1:length(ffts))])
    # maybe remove this filter so we get all components 
    #filter!(x -> (x[1] < x[2]) || (x[1]>x[2]), pairs) # filter upper diagonal
    if filt_dist; filter!(x -> get_dist(ffts[x[1]], ffts[x[2]]) <= 300., pairs); end # filter distances
    map(pair -> correlate_pair(ffts[pair[1]], ffts[pair[2]], prefix, params), pairs)
end
function offdiag_chunks(chunkpair::Array{Array{String, 1}}, prefix::String = "CORR", filt_dist::Bool = true, 
                        params::Dict=params) 
    sources   = map(x -> load_fft(x, string(x[end-7:end-5])), chunkpair[1]) # load
    receivers = map(x -> load_fft(x, string(x[end-7:end-5])), chunkpair[2]) 
    pairs = vec([collect(x) for x in Iterators.product(1:length(sources),1:length(receivers))]) # all pairs src -> rec
    if filt_dist; filter!(x -> get_dist(sources[x[1]].loc, receivers[x[2]].loc) <= 300., pairs); end # filt distances
    map(pair -> correlate_pair(sources[pair[1]], receivers[pair[2]], prefix, params), pairs)
end
function get_blocks(paths::Array{String,1})
    # Chunk raw ffts into blocks. 30 is best size for IO, but distribute across cores if more cores available
    chunks = collect(Iterators.partition(paths, convert(Int64, minimum([30, ceil(length(paths)/nprocs())])))) 
    off_chunks = vec([collect(x) for x in Iterators.product(1:length(chunks),1:length(chunks))])
    filter!(chunk_pair -> chunk_pair[1] < chunk_pair[2], off_chunks)
    # map chunk indices to filenames
    off_chunk_names = map(chunk -> [collect(chunks[chunk[1]]), collect(chunks[chunk[2]])], off_chunks)
    return chunks, off_chunk_names
end

# preprocessing
function preprocess2(file::String, accelerometer::Bool=false, params::Dict=params)
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
        gaps = size(data.t[1])[1] # N-2 gaps (eg gaps = 12 tests for 10 gaps)
        pts = size(data.x[1])[1]
        fs_temp = data.fs[1]
        windows = u2d.(SeisIO.t_win(data.t[1], data.fs[1]) .* 1e-6)
        startend = windows[:,2] .- windows[:,1] .> Second(params["cc_len"])
        if gaps < 25 && any(startend) # If few gaps and sufficient (2/3) data present, data is clean
            try
                add_location(data, params["all_stations"])
                for samp_rate in params["samp_rates"]
                    S = SeisData()
                    try
                        S = process_raw(data, samp_rate)
                    catch # trying to sample non 100Hz data at 100 Hz - skip resampling 
                        S = process_raw(data, data.fs[1])
                    end
                    R = RawData(S,params["cc_len"],params["cc_step"])
                    SeisNoise.detrend!(R)
                    bandpass!(R,params["freqmin"],params["freqmax"],zerophase=true)
                    SeisNoise.taper!(R)
                    FFT = nothing
                    if accelerometer # accelerometer - so we integrate
                        FFT = rfft_raw(R,1)
                    else # then its a seismometer
                        FFT = compute_fft(R)
                    end
                    coherence!(FFT,params["half_win"], params["water_level"])
                    try # save fft 
                        root_fft = "ffts/$path/$(Int(samp_rate))/"
                        save_fft(FFT, joinpath(params["rootdir"], root_fft))
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

# stacking codes
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
function stack_all()
    autocorr_names = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("AUTOCORR/*")]
    Sauto = @elapsed pmap(x -> stack_auto(x, startdate), autocorr_names)

    corr_names = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("CORR_20HZ/*")]
    S20 = @elapsed pmap(x -> stack_corr(x, startdate), corr_names)

    corr_names_lf = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("CORR_1HZ/*")]
    S1 = @elapsed pmap(x -> stack_corr(x, startdate), corr_names_lf)
    println("Stacking Completed for autocorrs and interstation cross correlations in $(Sauto+S20+S1) seconds")
end

# wrapper function for correlations
function correlate_big(dd::Date, startdate::Date = startdate, params::Dict = params)
    """ Wrapper Function: Computes autocorrelations for a specific day"""
    aws = params["aws"]
    yr = Dates.year(dd)
    path = join([Dates.year(dd),lpad(Dates.dayofyear(dd),3,"0")],"_") # Yeilds "YEAR_JDY"
    @eval @everywhere path = $path

    ###### WRAP THIS MESS #######
    # get filepaths for source stations - would be faster if they're available
    scedc_files = nothing
    if !isdir("scedc_path/$yr/"); mkpath("scedc_path/$yr/"); end
    try # most filelists are stored on seisbasin
        s3_get_file(aws, "seisbasin", "scedc_path/$yr/$path.csv", "scedc_path/$yr/$path.csv")
        scedc_files = DataFrame(CSV.File("scedc_path/$yr/$path.csv")).Path
    catch e # in case file not found/present we use SCEDC lookup
        println(e)
        scedc_files = get_scedc_files(dd, aws)
    end
    #########################

    # filepaths for iris, ncedc, and node data
    iris_path  = S3Path("s3://seisbasin/iris_waveforms/$yr/$path/", config=aws)
    ncedc_path = S3Path("s3://seisbasin/ncedc_waveforms/$yr/$path/", config=aws)
    node_path  = S3Path("s3://seisbasin/continuous_waveforms/$yr/$path/", config=aws)
    iris_query, ncedc_query, node_query = convert.(String, readdir(iris_path)), convert.(String, readdir(ncedc_path)), 
                                        convert.(String, readdir(node_path))
    filelist_basin = vcat(joinpath.("iris_waveforms/$yr/$path/", iris_query), joinpath.("ncedc_waveforms/$yr/$path/", ncedc_query),
                                        joinpath.("continuous_waveforms/$yr/$path/", node_query))
    println("There are $(length(filelist_basin)) node files and $(length(scedc_files)) SCEDC filesavailable for $path.")

    # download scedc and seisbasin data
    ec2download(aws, "scedc-pds", scedc_files, "~/data")
    ec2download(aws, "seisbasin", filelist_basin, "~/data")
    println("Download complete!")


    # preprocess - broadbands and seismometers are processed separately
    allf = glob("data/*/$yr/$path/*")
    accelerometers = filter(x -> isfile(x), joinpath.("data/continuous_waveforms/$yr/$path/", node_query))
    broadbands = setdiff(allf, accelerometers) # get the rest of the files 

    T_b = @elapsed pmap(f -> preprocess2(f, false, params), broadbands)
    if length(accelerometers) != 0 # save time on pmap allocation overhead
        T_a = @elapsed pmap(f -> preprocess2(f, true, params), accelerometers) # integrate accelerometers
        println("Preprocessing Completed in $(T_a+T_b) seconds.")
    else
        println("Preprocessing Completed in $T_b seconds.")
    end

    # autocorrelations    
    fft_list_100 = glob("ffts/$path/100/*.jld2")
    fft_100_stations = unique([join(split(elt, ".")[1:2],".") for elt in fft_list_100]) # get stations to iterate over
    pmap(x -> autocorrelate(x, fft_list_100, params), fft_100_stations)


    fft_paths_20 = sort(glob("ffts/$path/20/*")) # alphabetic sort for all stations so that diagonal files are made correctly
    fft_paths_1 = sort(glob("ffts/$path/1/*")) # alphabetic sort for all stations so that diagonal files are made correctly

    # run 20 HZ correlations
    chunks_20HZ, off_chunk_names_20HZ = get_blocks(fft_paths_20);

    # correlate diagonal chunks
    T20D = @elapsed pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_20HZ", true, params), chunks_20HZ)

    # correlate off-diagonal chunks
    T20O = @elapsed pmap(chunk -> offdiag_chunks(chunk, "CORR_20HZ", true, params), off_chunk_names_20HZ) # run mega correlations



    # run low freq correlations 
    chunks_1HZ, off_chunk_names_1HZ = get_blocks(fft_paths_1);

    # correlate diagonal chunks
    T1D = @elapsed pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_1HZ", false, params), chunks_1HZ)

    # correlate off-diagonal chunks
    T1O = @elapsed pmap(chunk -> offdiag_chunks(chunk, "CORR_1HZ", false, params), off_chunk_names_1HZ) # run mega correlations

    # Correlation summary
    println("All $(length(glob("CORR_*/*/*/$path*"))) Inter-station Correlations computed in $(T20D+T20O + T1D + T1O) seconds")
    rm("data/continuous_waveforms", recursive=true)
end

