export preprocess, correlate_block, correlate_day, stack_h5

####################### Preprocessing routine ##########################
function preprocess(file::String,  accelerometer::Bool=false, params::Dict=params)
    """
        Load raw seisdata file and process, saving to fft per frequency
    """
    # unpack needed params
    rootdir, samp_rate = params["rootdir"], params["fs"]
    freqmin, freqmax, cc_step, cc_len = params["freqmin"], params["freqmax"], params["cc_step"], params["cc_len"]
    half_win, water_level = params["half_win"], water_level["water_level"]

    try
        data = SeisData()
        try # assume mseed
            read_data!(data, "mseed", file)
        catch
            read_data!(data, "seisio", file)
        end
        gaps = size(data.t[1])[1] # N-2 gaps (eg gaps = 12 tests for 10 gaps)
        pts = size(data.x[1])[1]
        fs_temp = data.fs[1]
        windows = u2d.(SeisIO.t_win(data.t[1], data.fs[1]) .* 1e-6)
        startend = windows[:,2] .- windows[:,1] .> Second(cc_len)
        if gaps < 25 && any(startend) # If few gaps and sufficient (2/3) data present, data is clean
            try
                add_location(data, all_stations)
                S = SeisData()
                try
                    S = process_raw(data, samp_rate)
                catch # trying to sample non 100Hz data at 100 Hz - skip resampling 
                    S = process_raw(data, data.fs[1])
                end
                R = RawData(S,cc_len,cc_step)
                SeisNoise.detrend!(R)
                bandpass!(R,freqmin,freqmax,zerophase=true)
                SeisNoise.taper!(R)
                FFT = nothing
                if accelerometer # accelerometer - so we integrate
                    FFT = rfft_raw(R,1)
                else # then its a seismometer
                    FFT = compute_fft(R)
                end
                coherence!(FFT,half_win, water_level) # try no coherence ??
                try # save fft 
                    root_fft = joinpath(rootdir, "ffts/$path/")
                    save_fft(FFT, root_fft)
                catch e
                    println(e)
                end
                return data.id[1]
                println("Successfully processed $(data.id[1])")
            catch e
                println(e)
            end
        end
    catch e 
        println(e)
    end
end

############## Correlate two blocks of day waveforms ###################
function correlate_block(src::Array{String,1}, rec::Array{String,1}, maxlag::Float64)
    """ 
        Correlation function for pair of fft data
        - noise filter and stacking
        - saves directly to disk: ensure directory is correct if not on AWS
    """
    # load sources and recievers
    sources = map(x-> load_fft(x, string(x[end-7:end-5])), src)
    receivers = map(x-> load_fft(x, string(x[end-7:end-5])), rec)
    for source in sources
        for receiver in receivers
            try
                C = correlate(source,receiver,maxlag)
                cc_medianmute!(C, 10.) # remove correlation windows with high noise # try run with medianmute - robust stack on daily
                stack!(C)
                pair, comp = name_corr(C), C.comp
                save_named_corr(C,"CORR/$pair/$comp")
            catch e
                # add counter - iterate num bad; print
                println(e)
            end
        end
    end
end

#################### Downloads and Correlates day of data ###########################
function correlate_day(dd::Date, params::Dict=params)
    """ Wrapper function for daily correlations"""
    path = join([Dates.year(dd),lpad(Dates.dayofyear(dd),3,"0")],"_") # Yeilds "YEAR_JDY"
    @eval @everywhere path = $path

    # unpack needed params
    sources, maxlag, OUTDIR, yr, aws = params["sources"], params["maxlag"], params["OUTDIR"], params["yr"], params["aws"]

    # get filepaths for SCEDC source stations 
    scedc_files = get_scedc_files(dd, params)
    filter!(x -> any(occursin.(sources, x)), scedc_files)

    # filepaths for nodes
    #filelist_b = [f["Key"] for f in s3_list_objects(aws, "seisbasin", "continuous_waveforms/$(yr)/$(path)/")]
    filelist_b = S3Path("s3://seisbasin/continuous_waveforms/$(yr)/$(path)/", config=aws) # use new AWS functions
    filelist_basin = joinpath.("seisbasin/continuous_waveforms/$(yr)/$(path)/", 
                                convert.(String, readdir(filelist_b))) # add directory 
    println(typeof(filelist_basin))
    println("There are $(length(filelist_basin)) node files available for $path")
    # download scedc and seisbasin data
    try
        ec2download(aws, scedc, scedc_files, OUTDIR)
        ec2download(aws, basin, filelist_basin, OUTDIR)
        println("Download complete!")
    catch e
        println("Failed to process for $path. Continuing to next day.")
        return 1
    end
    # preprocess data
    allf = glob("data/continuous_waveforms/$yr/$path/*")
    broadbands = filter(x -> any(occursin.(sources, x)) && !occursin("Q0066",x), allf)
    accelerometers = [f for f in allf if !any(occursin.(f, broadbands))]

    T_b = @elapsed pmap(f -> preprocess(f, false, params), broadbands)
    T_a = @elapsed pmap(f -> preprocess(f, true, params), accelerometers)
    println("Preprocessing Completed in $(T_b + T_a) seconds.")

    # get indices for block correlation
    fft_paths = glob("ffts/$path/*")
    sources = filter(f -> any(occursin.(sources, f)), fft_paths)
    recievers = filter(f -> !any(occursin.(sources, f)), fft_paths)
    reciever_blocks = collect(Iterators.partition(recievers, convert(Int64, ceil(length(recievers)/nprocs())))) # chunk into 25 recievers to streamline IO and speed computation 
    println("Now Correlating $(length(reciever_blocks)) correlation blocks!")
    Tcorrelate = @elapsed pmap(rec_files -> correlate_block(sources, collect(rec_files), maxlag), reciever_blocks)
    rm("data/continuous_waveforms", recursive=true) # cleanup raw data
end

######################## Stacking Routine ############################
function stack_h5(tf::String, postfix::String)
    """Stacks all files for source-reciever pair, saves to h5 file by source"""

    # scrape correlation filenames - get unique sources and receiver combinations
    from_s = glob("CORR/$tf*/*/*")
    found_receivers = unique([join(split(f,".")[4:5],".") for f in from_s])
    components = ["EE","EN","EZ", "NE", "NN","NZ", "ZE", "ZN", "ZZ"]
    fname = join([tf, postfix, ".h5"])
    filename = joinpath("nodestack/$yr", fname)
    if !isdir(dirname(filename))
        mkpath(dirname(filename))
    end
    # create file and open
    h5open(filename, "cw") do file
        source_station = split(tf,".")[2]
        source_loc = LLE_geo(source_station, all_stations)
        samplef = glob("CORR/$tf*/ZZ/*")[1]
        C = load_corr(samplef, "ZZ")
        # add metadata information about correlation processing and source location 
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
            if !isnothing(source_loc)
                write(file, "meta/lon", source_loc.lon)
                write(file, "meta/lat", source_loc.lat)
                write(file, "meta/el", source_loc.el)
            end
        end
        # iterate through node recievers 
        for rec in found_receivers
            try
                # sample_r = glob("CORR/$tf*$rec*/ZZ/*")[1]
                # Cr = load_corr(sample_r, "ZZ")
                rec_station = split(rec, ".")[2]
                rec_loc = LLE_geo(rec_station, all_stations)
                if !haskey(read(file), "$rec/meta")
                    write(file, "$rec/meta/lon", rec_loc.lon)
                    write(file, "$rec/meta/lat", rec_loc.lat)
                    write(file, "$rec/meta/el", rec_loc.el)
                    write(file, "$rec/meta/dist", get_dist(source_loc, rec_loc))
                    write(file, "$rec/meta/azi", get_azi(source_loc, rec_loc))
                    write(file, "$rec/meta/baz", get_baz(source_loc, rec_loc))
                end
                for comp in components
                    try
                        # load correlations for this receiver by component 
                        files = glob("CORR/$tf*$rec*/$comp/*.jld2")
                        if length(files) == 0; continue; end
                        corrs = [load_corr(f, comp) for f in files]
                        corrs = [c for c in corrs if typeof(c) == CorrData]
        
                        # implement linear and phase-weighted stack methods. 
                        corr_sum = sum_corrs(corrs) # combine correlations into single CorrData object 
                        if isnothing(corr_sum); continue; end # if corr_sum fails to return, skip component

                        dims = size(corr_sum.corr)
                        cc_medianmute2!(corr_sum, 3., false) # skip metadata copy which errors on sum_corrs
                        println("$(length(corrs)) day correlations. Dims: $dims. Threw out $(dims[2]-size(corr_sum.corr)[2]) day windows.")

                        # then stack
                        corr_mean = SeisNoise.stack(corr_sum, allstack=true, stacktype=mean)
                        corr_pws = SeisNoise.stack(corr_sum, allstack=true, stacktype=pws)

                        # save to file 
                        write(file, "$rec/$comp/linear", corr_mean.corr[:])
                        write(file, "$rec/$comp/pws", corr_pws.corr[:])
                    catch e 
                        println("Error in stacking $rec $comp: ", e)
                    end
                end
            catch e
                println(e)
            end
        end
    end
end
