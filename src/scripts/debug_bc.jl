using SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWS, Distributed, JLD2, Glob, AWSS3, HDF5, Statistics

# add AWS, remove and add SCEDC - restart repl to re-precompile updated packages
# using Pkg; Pkg.add(PackageSpec(url="https://github.com/tclements/SCEDC.jl", rev="master"))
 # update AWSS3, AWS
addprocs()
# add Constants and locations dataframe
@everywhere begin 
    using SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWS, JLD2, Glob, AWSS3, HDF5, Statistics
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
    samp_rates = [1., 20., 100.]
    num_procs = nprocs()
    all_stations = DataFrame(CSV.File("files/updated_sources.csv"))
    params = Dict("aws" => aws, "cc_step" => cc_step, "cc_len" => cc_len, "maxlag" => maxlag,
    "fs" => fs, "half_win" => half_win, "water_level" => water_level, 
    "all_stations" => all_stations, "samp_rates" => samp_rates, "rootdir" => rootdir,
    "OUTDIR" => OUTDIR, "num_procs"=> num_procs, "freqmin" => freqmin, "freqmax" => freqmax)
end

@everywhere begin
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
    function get_blocks(paths::Array{String,1}, params::Dict=params)
        # Chunk raw ffts into blocks. 30 is best size for IO, but distribute across cores if more cores available
        chunks = collect(Iterators.partition(paths, 30))#convert(Int64, minimum([30, ceil(length(paths)/params["num_procs"])]))))
        off_chunks = vec([collect(x) for x in Iterators.product(1:length(chunks),1:length(chunks))])
        filter!(chunk_pair -> chunk_pair[1] < chunk_pair[2], off_chunks)
        # map chunk indices to filenames
        off_chunk_names = map(chunk -> [collect(chunks[chunk[1]]), collect(chunks[chunk[2]])], off_chunks)
        return chunks, off_chunk_names
    end
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
    function preprocess2(file::String, accelerometer::Bool=false, params::Dict=params, path::String=path)
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
                        catch e# trying to sample non 100Hz data at 100 Hz - skip resampling 
                            println(e)
                            S = process_raw(data, data.fs[1])
                        end
                        R = RawData(S,params["cc_len"],params["cc_step"])
                        SeisNoise.detrend!(R)
                        bandpass!(R,params["freqmin"],params["freqmax"],zerophase=true)
                        # SMOOTHING LINE GOES HERE
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
                    return data.id
                catch e
                    println(e)
                end
            end
        catch e 
            println(e)
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
    function add_location(s::SeisData,df::DataFrame)
        """ Adds locations to SeisData object from a dataframe """
        try
            name = split(s.id[1],".")[2]
            geo = LLE_geo(name, df)
            if !isnothing(geo)
                s.loc[1] = geo
            else 
                println("Station $name can't be found in the dataframe")
            end
        catch e
            println(e)
        end
    end
    function get_scedc_files(dd::Date, aws::AWSConfig)
        ar_filelist = map(x -> s3query(aws, dd, enddate = dd, network="CI", channel=x),["BH?", "HH?"])
        filelist_scedc_BH = ar_filelist[1]
        filelist_scedc_HH = ar_filelist[2]
        # create dictionary and overwrite HH keys with available BH data
    
        BH_keys = [get_dict_name(file) for file in filelist_scedc_BH]
        HH_keys = [get_dict_name(file) for file in filelist_scedc_HH]
    
        # Convert to dictionary 
        HH_dict = Dict([(name, file) for (name, file) in zip(HH_keys, filelist_scedc_HH)]) 
        BH_dict = Dict([(name, file) for (name, file) in zip(BH_keys, filelist_scedc_BH)]) 
        filelist_dict = merge(HH_dict, BH_dict) # BH dict overwrite HH_dict. This is essentually the union
        filelist_scedc = collect(values(filelist_dict)) # return values as array for download
        return filelist_scedc
    end
    function get_dict_name(file::String)
        """ Helper function for get_scedc_files."""
        station = convert(String, split(split(file,"/")[end],"_")[1])
        component = split(file,"/")[end][10:12]
        return string(station, "_", component)
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
            # stack per month 
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
end
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
    println("There are $(length(filelist_basin)) node files and $(length(scedc_files)) SCEDC files available for $path.")

    # download scedc and seisbasin data
    ec2download(aws, "scedc-pds", scedc_files, "~/data")
    ec2download(aws, "seisbasin", filelist_basin, "~/data")
    println("Download complete!")


    # preprocess - broadbands and seismometers are processed separately
    allf = glob("data/*/$yr/$path/*")
    accelerometers = filter(x -> isfile(x), joinpath.("data/continuous_waveforms/$yr/$path/", node_query))
    broadbands = setdiff(allf, accelerometers) # get the rest of the files 

    T_b = @elapsed pmap(f -> preprocess2(f, false, params, path), broadbands)
    if length(accelerometers) != 0 # save time on pmap allocation overhead
        T_a = @elapsed pmap(f -> preprocess2(f, true, params, path), accelerometers) # integrate accelerometers
        println("Preprocessing Completed in $(T_a+T_b) seconds.")
    else
        println("Preprocessing Completed in $T_b seconds.")
    end

    # autocorrelations    
    fft_list_100 = glob("ffts/$path/100/*.jld2")
    println("FFT list 100 is $(length(fft_list_100)) stations long")
    fft_100_stations = unique([join(split(elt, ".")[1:2],".") for elt in fft_list_100]) # get stations to iterate over
    pmap(x -> autocorrelate(x, fft_list_100, params), fft_100_stations)


    fft_paths_20 = sort(glob("ffts/$path/20/*")) # alphabetic sort for all stations so that diagonal files are made correctly
    fft_paths_1 = sort(glob("ffts/$path/1/*")) # alphabetic sort for all stations so that diagonal files are made correctly
    println("There are $(length(fft_paths_20)) fft datas for the 20 HZ correlations")
    # run 20 HZ correlations
    chunks_20HZ, off_chunk_names_20HZ = get_blocks(fft_paths_20, params);

    # correlate diagonal chunks
    T20D = @elapsed pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_20HZ", true, params), chunks_20HZ)

    # correlate off-diagonal chunks
    T20O = @elapsed pmap(chunk -> offdiag_chunks(chunk, "CORR_20HZ", true, params), off_chunk_names_20HZ) # run mega correlations



    # run low freq correlations 
    chunks_1HZ, off_chunk_names_1HZ = get_blocks(fft_paths_1, params);

    # correlate diagonal chunks
    T1D = @elapsed pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_1HZ", false, params), chunks_1HZ)

    # correlate off-diagonal chunks
    T1O = @elapsed pmap(chunk -> offdiag_chunks(chunk, "CORR_1HZ", false, params), off_chunk_names_1HZ) # run mega correlations

    # Correlation summary
    println("All $(length(glob("CORR_*/*/*/$path*"))) Inter-station Correlations computed in $(T20D+T20O + T1D + T1O) seconds")
    try; rm("data/continuous_waveforms/", recursive=true); catch e; println(e); end
end

arg = 0
startdate = Date(2004)+Month(arg)
enddate = startdate+Month(1)-Day(1)
days = startdate:Day(1):enddate
println("Processing download for: ", startdate, " to ",enddate)

Tstart = Dates.now()
Tcbig = @elapsed map(dd -> correlate_big(dd, startdate, params), days[3:4])


dd = days[1]
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
println("There are $(length(filelist_basin)) node files and $(length(scedc_files)) SCEDC files available for $path.")

# download scedc and seisbasin data
ec2download(aws, "scedc-pds", scedc_files, "~/data")
ec2download(aws, "seisbasin", filelist_basin, "~/data")
println("Download complete!")


# preprocess - broadbands and seismometers are processed separately
allf = glob("data/*/$yr/$path/*")
accelerometers = filter(x -> isfile(x), joinpath.("data/continuous_waveforms/$yr/$path/", node_query))
broadbands = setdiff(allf, accelerometers) # get the rest of the files 

T_b = @elapsed map(f -> preprocess2(f, false, params, path), broadbands[1:3])
if length(accelerometers) != 0 # save time on pmap allocation overhead
    T_a = @elapsed pmap(f -> preprocess2(f, true, params, path), accelerometers) # integrate accelerometers
    println("Preprocessing Completed in $(T_a+T_b) seconds.")
else
    println("Preprocessing Completed in $T_b seconds.")
end

# autocorrelations    
fft_list_100 = glob("ffts/$path/100/*.jld2")
println("FFT list 100 is $(length(fft_list_100)) stations long")
fft_100_stations = unique([join(split(elt, ".")[1:2],".") for elt in fft_list_100]) # get stations to iterate over
pmap(x -> autocorrelate(x, fft_list_100, params), fft_100_stations)


fft_paths_20 = sort(glob("ffts/$path/20/*")) # alphabetic sort for all stations so that diagonal files are made correctly
fft_paths_1 = sort(glob("ffts/$path/1/*")) # alphabetic sort for all stations so that diagonal files are made correctly
println("There are $(length(fft_paths_20)) fft datas for the 20 HZ correlations")
# run 20 HZ correlations
chunks_20HZ, off_chunk_names_20HZ = get_blocks(fft_paths_20, params);

# correlate diagonal chunks
T20D = @elapsed pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_20HZ", true, params), chunks_20HZ)

# correlate off-diagonal chunks
T20O = @elapsed pmap(chunk -> offdiag_chunks(chunk, "CORR_20HZ", true, params), off_chunk_names_20HZ) # run mega correlations



# run low freq correlations 
chunks_1HZ, off_chunk_names_1HZ = get_blocks(fft_paths_1, params);

# correlate diagonal chunks
T1D = @elapsed pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_1HZ", false, params), chunks_1HZ)

# correlate off-diagonal chunks
T1O = @elapsed pmap(chunk -> offdiag_chunks(chunk, "CORR_1HZ", false, params), off_chunk_names_1HZ) # run mega correlations

# Correlation summary
println("All $(length(glob("CORR_*/*/*/$path*"))) Inter-station Correlations computed in $(T20D+T20O + T1D + T1O) seconds")








autocorr_names = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("AUTOCORR/*")]
Sauto = @elapsed pmap(x -> stack_auto(x, startdate), autocorr_names)

corr_names = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("CORR_20HZ/*")]
S20 = @elapsed pmap(x -> stack_corr(x, startdate, "CORR_20HZ"), corr_names)

corr_names_lf = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("CORR_1HZ/*")]
S1 = @elapsed pmap(x -> stack_corr(x, startdate, "CORR_1HZ"), corr_names_lf)
println("Stacking Completed for autocorrs and interstation cross correlations in $(Sauto+S20+S1) seconds")

# transfer to s3 
autos = glob("autocorrelations/*/*")
bigcorrs = glob("correlations/*/*")

pmap(x -> s3_put(aws, "seisbasin", x, read(x), acl="bucket-owner-full-control"), autos)
pmap(x -> s3_put(aws, "seisbasin", x, read(x), acl="bucket-owner-full-control"), bigcorrs)
Tend = Dates.now()
totalt = Tend-Tstart
println("Finished in $totalt")

# scp -i my_aws_key.pem ubuntu@ec2-34-221-181-73.us-west-2.compute.amazonaws.com:correlations/AZ.BZN/2004_January_AZ.BZN.h5 ~/Downloads/CORR_20HZ_2004_January_AZ.BZN.h5


function preprocess2(file::String,  accelerometer::Bool=false, rootdir::String="", samp_rates::Array{Float64, 1}=[1., 20., 100.],
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