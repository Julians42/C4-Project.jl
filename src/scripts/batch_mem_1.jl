using Distributed, Parallelism


addprocs()
# add Constants and locations dataframe
@everywhere begin 
    using Dates, CSV, DataFrames, Statistics, AbstractFFTs, Serialization, Glob
    using  SeisIO, SeisNoise, SCEDC, AWS, AWSS3, HDF5, JLD2
    aws = AWS.AWSConfig(region="us-west-2")
    rootdir = "/scratch" # for docker ecs image we have added 750 GB to docker scratch container: "/scratch"
    network = "CI"
    channel1 = "BH?"
    channel2 = "HH?"
    OUTDIR = "~/data"

    # Figure out dates
    arg = ENV["AWS_BATCH_JOB_ARRAY_INDEX"]
    startdate = Date(2004)+Month(arg)
    enddate = startdate+Month(1)-Day(1)
    days = startdate:Day(1):enddate
    month = Dates.month(startdate)
    
    num_procs = nprocs()
    cc_step, cc_len = 3600, 3600
    maxlag, fs = 1200., 20. # maximum lag time in correlation, sampling frequency
    freqmin, freqmax = 0.05, 9.9
    half_win, water_level = 30, 0.01
    samp_rates = [1., 20., 100.] # for processing
    #all_stations = DataFrame(CSV.File("/home/ubuntu/SeisCore.jl/docs/updated_sources.csv"))
    all_stations = DataFrame(CSV.File("/root/files/CAstations.csv"))
    params = Dict("aws" => aws, "cc_step" => cc_step, "cc_len" => cc_len, "maxlag" => maxlag,
            "fs" => fs, "half_win" => half_win, "water_level" => water_level, "month" => month,
            "all_stations" => all_stations, "samp_rates" => samp_rates, "rootdir" => rootdir,
            "OUTDIR" => OUTDIR, "num_procs"=> num_procs, "freqmin" => freqmin, "freqmax" => freqmax)
end

@everywhere begin
    # correlation helper functions
    function correlate_pair(src::FFTData, rec::FFTData, pref::String="CORR", params::Dict=params)
        try
            C = correlate(src, rec, params["maxlag"])
            cc_medianmute!(C, 10.) # remove correlation windows with high noise
            stack!(C)
            name = join(split(name_corr(C), ".")[1:2],".")
            save_named_corr(C,"$(params["rootdir"])/$pref/$name/$(C.comp)")
        catch e
            println(e)
        end
    end
    function autocorrelate(station::String, list::Array{String,1}=fft_list_100, params::Dict=params)
        """ Filter channels by station, load and correlate """
        channels = filter(x -> occursin(station, x), list)
        ffts = deserialize.(channels)
        pairs = vec([collect(x) for x in Iterators.product(1:length(ffts),1:length(ffts))])
        map(pair -> correlate_pair(ffts[pair[1]], ffts[pair[2]], "AUTOCORR", params), pairs)
    end
    function diag_chunks(chunk::Array{String,1}, prefix::String = "CORR", filt_dist::Bool = true, params::Dict=params)
        ffts  = deserialize.(chunk)
        pairs = vec([collect(x) for x in Iterators.product(1:length(ffts),1:length(ffts))])
        # maybe remove this filter so we get all components 
        #filter!(x -> (x[1] < x[2]) || (x[1]>x[2]), pairs) # filter upper diagonal
        if filt_dist; filter!(x -> get_dist(ffts[x[1]], ffts[x[2]]) <= 300., pairs); end # filter distances
        map(pair -> correlate_pair(ffts[pair[1]], ffts[pair[2]], prefix, params), pairs)
    end
    function offdiag_chunks(chunkpair::Array{Array{String, 1}}, prefix::String = "CORR", filt_dist::Bool = true, 
                            params::Dict=params) 
        sources, receivers = deserialize.(chunkpair[1]), deserialize.(chunkpair[2]) # load
        pairs = vec([collect(x) for x in Iterators.product(1:length(sources),1:length(receivers))]) # all pairs src -> rec
        if filt_dist; filter!(x -> get_dist(sources[x[1]].loc, receivers[x[2]].loc) <= 300., pairs); end # filt distances
        map(pair -> correlate_pair(sources[pair[1]], receivers[pair[2]], prefix, params), pairs)
    end
    function get_blocks(paths::Array{String,1}, params::Dict=params)
        # Chunk raw ffts into blocks. 30 is best size for IO, but distribute across cores if more cores available
        chunks = collect(Iterators.partition(paths, 30))#convert(Int64, minimum([30, ceil(length(paths)/params["num_procs"])]))))
        off_chunks = vec([collect(x) for x in Iterators.product(1:length(chunks),1:length(chunks))])
        filter!(chunk_pair -> chunk_pair[1] < chunk_pair[2], off_chunks)
        # map chunk indices to filenames
        off_chunk_names = map(chunk -> [collect(chunks[chunk[1]]), collect(chunks[chunk[2]])], off_chunks)
        return chunks, off_chunk_names
    end
    # preprocessing
    function preprocess(file::String, params::Dict=params)
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
                    names = f = Array{Union{Nothing,String}, 1}(nothing,3)
                    for (ind, samp_rate) in enumerate(params["samp_rates"])
                        S = SeisData()
                        try
                            S = process_raw(data, samp_rate)
                        catch e# trying to sample non 100Hz data at 100 Hz - skip resampling 
                            println(e)
                            S = process_raw(data, data.fs[1])
                        end
                        R = RawData(S,params["cc_len"],params["cc_step"])
                        SeisNoise.detrend!(R)
                        bandpass!(R,params["freqmin"],params["freqmax"],zerophase=true)
                        # SMOOTHING LINE GOES HERE
                        SeisNoise.taper!(R)
                        FFT = compute_fft(R)
                        coherence!(FFT,params["half_win"], params["water_level"])
                        try # save fft 
                            #root_fft = "ffts/$path/$(Int(samp_rate))/"
                            #save_fft(FFT, joinpath(params["rootdir"], root_fft))
                            serial_fft = "$(Int(samp_rate))_$(FFT.name)"
                            serialize(serial_fft, FFT)
                            names[ind] = serial_fft # update names if serialization is successful
                        catch e
                            println(e)
                        end
                        #println("Successfully processed $(data.id[1])")
                    end
                    return names
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

        CORROUT = expanduser(joinpath(params["rootdir"],"autocorrelations/$name/"))
        if !isdir(CORROUT)
            mkpath(CORROUT)
        end
        filename = joinpath(CORROUT,"$(yr)_$(month)_$name.h5") # get output filename

        # Get list of files to save
        autocorr_list = glob("AUTOCORR/$name*/*/*.jld2", params["rootdir"])
        components = ["EE", "EN", "EZ", "NE", "NN", "NZ", "ZE", "ZN", "ZZ"]

        C = load_corr(autocorr_list[1], convert(String, split(autocorr_list[1],"/")[end-1])) # sample autocorr for meta
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
            # stack per month 
            for comp in components
                #comp_files = filter(f -> convert(String, split(f, "/")[end-1]) == comp, autocorr_list)
                comp_files = glob("AUTOCORR/$name*/$comp/*.jld2", params["rootdir"]) # uni-component files
                autocorrs = [load_corr(f, comp) for f in comp_files]
                println("There are $(length(autocorrs)) separate correlations. Dims are $(size(sum(autocorrs).corr))")
                autocorr_mean = SeisNoise.stack(sum(autocorrs), allstack=true, stacktype=mean)
                autocorr_pws = SeisNoise.stack(sum(autocorrs), allstack=true, stacktype=pws)
                autocorr_robust = SeisNoise.stack(sum(autocorrs), allstack=true, stacktype=robuststack)
                # save into file 
                write(file, "$comp/linear", autocorr_mean.corr[:])
                write(file, "$comp/pws", autocorr_pws.corr[:])
                write(file, "$comp/robust", autocorr_robust.corr[:])
            end
        end
    end
    function stack_corr(name::String, startdate::Date=startdate, prefix::String = "CORR_20HZ")
        month = Dates.monthname(startdate) # get month for filename
        yr = Dates.year(startdate)
        # stack autocorrelations 
    
        CORROUT = expanduser(joinpath(params["rootdir"], "correlations/$name/"))
        if !isdir(CORROUT)
            mkpath(CORROUT)
        end
        filename = joinpath(CORROUT,"$(yr)_$(month)_$name.h5") # get output filename
        components = ["EE","EN","EZ", "NE", "NN","NZ", "ZE", "ZN", "ZZ"]
    
        receivers = glob("$prefix/$name/*/*", params["rootdir"])
        receivers = Set([join(split(x, ".")[end-2:end-1], ".") for x in receivers]) # just get reciever names
        corr_list = glob("$prefix/$name*/*/*", params["rootdir"])
        println("Processing $(length(corr_list)) 20HZ corelations")
        C = load_corr(corr_list[1], convert(String, split(corr_list[1],"/")[end-1]))
        source_loc = C.loc#GeoLoc(lat = C.loc.lat, lon = C.loc.lon, el = C.loc.el)
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
                    rec_network, rec_station = split(rec, ".")[1], split(rec, ".")[2]
                    #rec_loc = LLE_geo(rec_station, all_stations)
                    println(rec_station)
                    rec_loc = LLE_geo(rec_network, rec_station, all_stations)
                    #rec_loc = read_data("mseed", glob("data/*_waveforms/$yr/$yr*/$rec*HZ*", params["rootdir"])[1]).loc # get location from fft - perhaps a faster way to do this. 
                    println("starting the stack!")
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
                        files = glob("$prefix/$name/$comp/*$name..$rec.jld2", params["rootdir"])
                        corrs = [load_corr(f, comp) for f in files]
                        println("There are $(length(corrs)) separate correlations. Dims are $(size(sum(corrs).corr))")
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
                    if e isa(InterruptException)
                        rethrow(e)
                    else
                        println(e)
                    end
                end
            end
        end
    end
    function stack_all(params::Dict=params)
        autocorr_names = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("AUTOCORR/*", params["rootdir"])]
        #@eval @everywhere autocorr_names = $autocorr_names
        Sauto = @elapsed robust_pmap(x -> stack_auto(x, startdate), autocorr_names)

        corr_names = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("CORR_20HZ/*", params["rootdir"])]
        #@eval @everywhere corr_names = $corr_names
        S20 = @elapsed robust_pmap(x -> stack_corr(x, startdate), corr_names)

        corr_names_lf = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("CORR_1HZ/*", params["rootdir"])]
        #@eval @everywhere corr_names_lf = $corr_names_lf
        S1 = @elapsed robust_pmap(x -> stack_corr(x, startdate, "CORR_1HZ"), corr_names_lf)
        println("Stacking Completed for autocorrs and interstation cross correlations in $(Sauto+S20+S1) seconds")
    end
    function LLE_geo(network, station, df)
        """ Find station matching location and return geoloc object"""
        try
            row = df[findfirst(x -> x.Station==station && x.Network==network, eachrow(df)),:]
            lat, lon, el = row.Latitude[1], row.Longitude[1], row.Elevation[1]
            geo = GeoLoc(lat = float(lat), lon = float(lon), el = float(el))
            return geo
        catch 
            try # try lowercase csv format (older csvs)
                row = df[findfirst(x -> x.station==station && x.network==network, eachrow(df)),:]
                lat, lon, el = row.latitude[1], row.longitude[1], row.elevation[1]
                geo = GeoLoc(lat = float(lat), lon = float(lon), el = float(el))
                return geo
            catch
                return nothing
            end
        end
    end
    function add_location(s::SeisData,df::DataFrame)
        """ Adds locations to SeisData object from a dataframe """
        try
            network = split(s.id[1], ".")[1]
            station = split(s.id[1],".")[2]
            println(network, "  ", station)
            geo = LLE_geo(network, station, df)
            if !isnothing(geo)
                s.loc[1] = geo
            else 
                println("Station $network.$station can't be found in the dataframe")
            end
        catch e
            println(e)
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
function correlate_big(dd::Date, startdate::Date = startdate, params::Dict = params)
    """ Wrapper Function: Computes autocorrelations for a specific day"""
    println("Starting correlations for $dd.")
    aws, rootdir = params["aws"], params["rootdir"]
    yr = Dates.year(dd)
    path = join([Dates.year(dd),lpad(Dates.dayofyear(dd),3,"0")],"_") # Yeilds "YEAR_JDY"
    #@eval @everywhere path = $path

    ###### WRAP THIS MESS #######
    # get filepaths for source stations - would be faster if they're available
    scedc_files = nothing

    if !isdir(joinpath(rootdir, "scedc_path/$yr/")); mkpath(joinpath(rootdir, "scedc_path/$yr/")); end
    try # most filelists are stored on seisbasin
        s3_get_file(aws, "seisbasin", "scedc_path/$yr/$path.csv", joinpath(rootdir, "scedc_path/$yr/$path.csv"))
        scedc_files = DataFrame(CSV.File(joinpath(rootdir, "scedc_path/$yr/$path.csv"))).Path
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
    println("There are $(length(filelist_basin)) node files and $(length(scedc_files)) SCEDC files available for $path.")

    # download scedc and seisbasin data
    ec2download(aws, "scedc-pds", scedc_files, joinpath(rootdir,"data"))
    ec2download(aws, "seisbasin", filelist_basin, joinpath(rootdir,"data"))
    println("Download complete!")


    # preprocess - broadbands and seismometers are processed separately
    raw_waveforms = glob("data/*/$yr/$path/*", rootdir)# params["rootdir"])
    println("$(length(raw_waveforms)) waveforms available for processing.")
    #println(readdir("$(params["rootdir"])/data"))

    #@eval @everywhere raw_waveforms = $raw_waveforms
    T_b = @elapsed serial_names = robust_pmap(f -> preprocess(f, params), raw_waveforms)
    println("Preprocessing completed in $T_b seconds. $(length(serial_names)) ffts computed.")

    serial_names = collect(Iterators.flatten(filter(x -> !isnothing(x), serial_names)))
    serial_names = convert(Array{String, 1}, serial_names)

    # autocorrelations    
    fft_list_100 = filter(name -> occursin("100_", name), serial_names)
    fft_list_100 = convert(Array{String, 1}, fft_list_100)
    println("FFT list 100 is $(length(fft_list_100)) stations long")
    fft_100_stations = unique([split(join(split(elt, ".")[1:2],".") ,"_")[2] for elt in fft_list_100]) # get stations to iterate over
    #@eval @everywhere fft_100_stations = $fft_100_stations
    pmap(x -> autocorrelate(x, fft_list_100, params), convert.(String, fft_100_stations))


    fft_paths_20 = sort(filter(name -> occursin("20_", name), serial_names)) # alphabetic sort for all stations so that diagonal files are made correctly
    fft_paths_1 = sort(filter(name -> occursin("1_", name), serial_names)) # alphabetic sort for all stations so that diagonal files are made correctly
    println("There are $(length(fft_paths_20)) fft datas for the 20 HZ correlations")
    # run 20 HZ correlations
    chunks_20HZ, off_chunk_names_20HZ = get_blocks(fft_paths_20, params);

    # correlate diagonal chunks
    #@eval @everywhere chunks_20HZ = $chunks_20HZ
    T20D = @elapsed robust_pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_20HZ", true, params), chunks_20HZ)

    # correlate off-diagonal chunks
    #@eval @everywhere off_chunk_names_20HZ = $off_chunk_names_20HZ
    T20O = @elapsed robust_pmap(chunk -> offdiag_chunks(chunk, "CORR_20HZ", true, params), off_chunk_names_20HZ) # run mega correlations



    # run low freq correlations 
    chunks_1HZ, off_chunk_names_1HZ = get_blocks(fft_paths_1, params);

    # correlate diagonal chunks
    #@eval @everywhere chunks_1HZ = $chunks_1HZ
    T1D = @elapsed robust_pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_1HZ", false, params), chunks_1HZ)

    # correlate off-diagonal chunks
    #@eval @everywhere off_chunk_names_1HZ = $off_chunk_names_1HZ
    T1O = @elapsed robust_pmap(chunk -> offdiag_chunks(chunk, "CORR_1HZ", false, params), off_chunk_names_1HZ) # run mega correlations

    # Correlation summary
    println("All $(length(glob("CORR_*/*/*/$path*", rootdir))) Inter-station Correlations computed in $(T20D + T20O + T1D + T1O) seconds")
    try; rm("$(params["rootdir"])/data/continuous_waveforms/", recursive=true); catch e; println(e); end
    GC.gc()
end


println("Begin processing correlations for: ", startdate, " to ",enddate)

# Big wrapper function for processing and correlating everything
Tcbig = @elapsed map(dd -> correlate_big(dd, startdate, params), days[1:3])
println("All correlations completed in $Tcbig seconds.")

# Stack each of the 100, 20, and 1 Hz data
stack_all()


# transfer to s3
@everywhere begin
    autos = glob("autocorrelations/*/*", params["rootdir"])
    bigcorrs = glob("correlations/*/*", params["rootdir"])
end 

pmap(x -> s3_put(aws, "seisbasin", join(["autocorrelations", convert(String, split(x, "autocorrelations")[2])]), 
                read(x), acl="bucket-owner-full-control"), autos)
pmap(x -> s3_put(aws, "seisbasin", join(["correlations", convert(String, split(x, "correlations")[2])]),
                read(x), acl="bucket-owner-full-control"), bigcorrs)

println("Finished Correlations and upload for: ", startdate, " to ",enddate)