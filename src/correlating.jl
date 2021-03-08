export correlate_day, stack_h5, preprocess, correlate_block, LLE_geo, add_location, rfft_raw,
        cc_medianmute, cc_medianmute!, name_corr, save_named_corr, load_corrs, write_jld2,
        foldersize, get_dict_name, get_scedc_files


# main processing functions
function correlate_day(dd::Date, sources::Array{String,1}=sources, all_stations::DataFrame=all_stations)
    """ Wrapper function for daily correlations"""
    path = join([Dates.year(dd),lpad(Dates.dayofyear(dd),3,"0")],"_") # Yeilds "YEAR_JDY"
    @eval @everywhere path = $path
    # get filepaths for source stations - would be faster if they're available
    #s3_get_file(aws, "seisbasin", "scedc_path/$yr/$path.csv", "scedc_path/$yr/$path.csv")
    #scedc_files = DataFrame(CSV.File("scedc_path/$yr/$path.csv")).Path
    scedc_files = get_scedc_files(dd)
    filter!(x -> any(occursin.(sources, x)), scedc_files)

    # filepaths for nodes
    filelist_b = [f["Key"] for f in s3_list_objects(aws, "seisbasin", "continuous_waveforms/$(yr)/$(path)/")]
    filelist_basin = [convert(String, x) for x in filelist_b]
    println(typeof(filelist_basin))
    println("There are $(length(filelist_basin)) node files available for $path")
    # download scedc and seisbasin data
    ec2download(aws, scedc, scedc_files, OUTDIR)
    ec2download(aws, basin, filelist_basin, OUTDIR)
    println("Download complete!")
    # preprocess data
    allf = glob("data/continuous_waveforms/$yr/$path/*")
    broadbands = filter(x -> any(occursin.(sources, x)) && !occursin("Q0066",x), allf)
    accelerometers = [f for f in allf if !any(occursin.(f, broadbands))]

    T_b = @elapsed pmap(f -> preprocess(f), broadbands)
    T_a = @elapsed pmap(f -> preprocess(f, true), accelerometers)
    println("Preprocessing Completed")
    # get indices for block correlation
    fft_paths = glob("ffts/$path/*")
    sources = filter(f -> any(occursin.(sources, f)), fft_paths)
    recievers = filter(f -> !any(occursin.(sources, f)), fft_paths)
    reciever_blocks = collect(Iterators.partition(recievers, convert(Int64, ceil(length(recievers/nprocs()))))) # chunk into 25 recievers to streamline IO and speed computation 
    println("Now Correlating $(length(reciever_blocks)) correlation blocks!")
    Tcorrelate = @elapsed pmap(rec_files -> correlate_block(sources, collect(rec_files), maxlag), reciever_blocks)
    rm("data/continuous_waveforms", recursive=true) # cleanup raw data
end
function stack_h5(tf::String, postfix::String)
    """ Stack and reformat correlations from source `tf`. 
        Postfix is used for naming files (ex 2019, B1_nodes etc)
        Writes to H5 files"""
    from_s = glob("CORR/$tf*/*/*")
    found_receivers = unique([join(split(f,".")[4:5],".") for f in from_s])
    components = ["EE","EN","EZ","NN","NZ","ZZ"]
    fname = join([tf, postfix, ".h5"])
    filename = joinpath("nodestack/$yr", fname)
    if !isdir(dirname(filename))
        mkpath(dirname(filename))
    end
    h5open(filename, "cw") do file
        source_station = split(tf,".")[2]
        source_loc = LLE_geo(source_station, all_stations)
        samplef = glob("CORR/$tf*/ZZ/*")[1]
        C = load_corr(samplef, "ZZ")
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
                    # load correlations for this receiver by component 
                    files = glob("CORR/$tf*$rec*/$comp/*")
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

# subroutines
function preprocess(file::String,  accelerometer::Bool=false,
                        freqmin::Float64=freqmin, freqmax::Float64=freqmax, cc_step::Int64=cc_step, 
                        cc_len::Int64=cc_len, half_win::Int64=half_win, water_level::Float64=water_level, samp_rate::Float64=fs)
    """
        Load raw seisdata file and process, saving to fft
    """
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
                    root_fft = "ffts/$path/"
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
function correlate_block(src::Array{String,1}, rec::Array{String,1}, maxlag::Float64)
    """ 
        Correlation function for sources to receivers with stacking by day
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
                cc_medianmute!(C, 10.) # remove correlation windows with high noise
                stack!(C)
                pair, comp = name_corr(C), C.comp
                save_named_corr(C,"CORR/$pair/$comp")
            catch e
                println(e)
            end
        end
    end
end

# helper functions - simplify larger coding blocks above 
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
    """ Adds locations to array of seisdata from a dataframe """
    name = split(s.name[1],".")[2]
    geo = LLE_geo(name, df)
    if !isnothing(geo)
        s.loc[1] = geo
    else
        println("Station $name can't be found in the dataframe")
    end
end
function rfft_raw(R::RawData,dims::Int=1)
    FFT = rfft(R.x,dims)
    FFT ./= rfftfreq(size(R.x,1), R.fs) .* 1im .* 2π # Integrate the accelerometers!
    FFT[1,:] .=0
    return FFTData(R.name, R.id,R.loc, R.fs, R.gain, R.freqmin, R.freqmax,
                    R.cc_len, R.cc_step, R.whitened, R.time_norm, R.resp,
                    R.misc, R.notes, R.t, FFT)
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
function foldersize(dir=".")
    """ returns total size of folder in GB """
    size = 0
    for (root, dirs, files) in walkdir(dir)
        size += sum(map(filesize, joinpath.(root, files)))
    end
    return size*10e-10
end
function get_dict_name(file::String)
    """ Helper function for get_scedc_files."""
    station = convert(String, split(split(file,"/")[end],"_")[1])
    component = split(file,"/")[end][10:12]
    return string(station, "_", component)
end
function get_scedc_files(dd::Date)
    ar_filelist = pmap(x -> s3query(aws, dd, enddate = dd, network=network, channel=x),["BH?", "HH?"])
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