# Script for docker container to run correlations on AWS Batch
# This file runs on batch through the correlation step. Stacking has yet to be done
using Distributed, Parallelism


addprocs()
# add Constants and locations dataframe
@everywhere begin 
    using Dates, CSV, DataFrames, Statistics, AbstractFFTs, Serialization, Glob
    using  SeisIO, SeisNoise, SCEDC, AWS, AWSS3, HDF5, JLD2
    aws = AWS.AWSConfig(region="us-west-2")
    num_procs = nprocs()

    # batch filepathing
    rootdir = "/scratch" # BATCH: for docker ecs image we have added 2 TB to docker scratch container: "/scratch"
    arg = ENV["AWS_BATCH_JOB_ARRAY_INDEX"] # BATCH
    all_stations = DataFrame(CSV.File("root/files/CAstations.csv")) # BATCH - locations file filepathing

    # map arg to find job which hasn't been done 
    # succeeded = [2,8,10,12,24,25,29,33,36,42,43,47,49,50,51,57,58,59] # jobs which have succeeded
    # all = collect(0:59)
    # all = [elt for elt in all if elt ∉ succeeded]
    all = [12,14,15,16,17,18,19,20]
    arg = all[parse(Int64, arg)+1] # add 1 because julia 1 indexes

    #ec2 filepathing
    # arg = 0#ARGS[1] # EC2 - accept month from commandline
    # rootdir = "/home/ubuntu" # EC2 - write files here
    # all_stations = DataFrame(CSV.File("files/CAstations.csv")) # EC2 - filepathing directly from home dir

    XMLDIR = joinpath(rootdir, "XML")

    # Dates for computation 
    startdate = Date(2000)+Month(arg)
    enddate = startdate+Month(1)-Day(1)
    days = startdate:Day(1):enddate
    month, yr = Dates.monthname(startdate), Dates.year(startdate)
    
    # preprocessing and correlation coefficients
    cc_step, cc_len = 3600, 3600
    maxlag, fs = 1200., 20. # maximum lag time in correlation, sampling frequency
    freqmin, freqmax = 0.01, 9.9
    half_win, water_level = 30, 0.01
    samp_rates = [1., 20., 100.] # for processing

    params = Dict("aws" => aws, "cc_step" => cc_step, "cc_len" => cc_len, "maxlag" => maxlag,
            "fs" => fs, "half_win" => half_win, "water_level" => water_level, "month" => month,
            "all_stations" => all_stations, "samp_rates" => samp_rates, "rootdir" => rootdir, "yr" => yr,
            "num_procs"=> num_procs, "freqmin" => freqmin, "freqmax" => freqmax, "XMLDIR" => XMLDIR)
end

@everywhere begin
    # correlation helper functions
    function correlate_pair(src::FFTData, rec::FFTData, pref::String="CORR", params::Dict=params)
        try
            C = correlate(src, rec, params["maxlag"])
            cc_medianmute!(C, 10.) # remove correlation windows with high noise
            stack!(C)
            name = join(split(name_corr(C), ".")[1:2],".")
            save_named_corr(C,"$(params["rootdir"])/$pref/$(yr)_$(params["month"])/$name/$(C.comp)")
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
        # Filter part of extraneous correlations (still some loss without checking station names)
        filter!(x -> (x[1]-3 < x[2]), pairs) # filter upper diagonal for these components. -3 to get all comps for autos. 
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
    function process_raw2(S::SeisData, samp_rate::Real, file::String, params::Dict; ϕshift::Bool=true)
        merge!(S)
        ungap!(S)
        detrend!(S)         # remove mean & trend from channel
        SeisIO.taper!(S, t_max=100)                # taper channel ends: use SeisIO's function, t_max = 100
        if fs ∉ S.fs && samp_rate <= S.fs[1] # if we're downsampling
            filtfilt!(S,fh=Float64(samp_rate/2),rt="Lowpass")    # lowpass filter before downsampling
        end
        resample!(S,fs=Float64(samp_rate)) # resample based on product
        # extend tapering
        SeisIO.taper!(S, t_max=20)
        phase_shift!(S, ϕshift=ϕshift) # timing offset from sampling period
        if split(S.name[1], ".")[1] == "CI" # Check if in CI network (SCEDC Data)
            # add instrument response
            add_response(S, file, params["XMLDIR"])
        end
        return S
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
                    names = Array{Union{Nothing,String}, 1}(nothing,3)
                    for (ind, samp_rate) in enumerate(params["samp_rates"])
                        S = SeisData()
                        try
                            S = process_raw2(data, samp_rate, file, params)
                        catch e# trying to sample non 100Hz data at 100 Hz - skip resampling 
                            println(e)
                            S = process_raw2(data, data.fs[1], file, params)
                        end
                        highpass!(S, 0.01, corners=2)
                        # remove instrument response
                        remove_resp!(S)

                        R = RawData(S,params["cc_len"],params["cc_step"])
                        SeisNoise.detrend!(R)
                        # bandpass below nyquist or 9.9 Hz (not interested in high frequencies)
                        bandpass!(R,params["freqmin"], minimum([R.fs[1]/2-0.01, params["freqmax"]]), zerophase=true)
                        # SMOOTHING LINE GOES HERE
                        SeisNoise.taper!(R)
                        FFT = compute_fft(R)
                        coherence!(FFT,params["half_win"], params["water_level"])
                        try # save fft 
                            root_fft = joinpath(params["rootdir"], "ffts/")
                            #save_fft(FFT, joinpath(params["rootdir"], root_fft))
                            serial_fft = "$(yr)_$(params["month"])_$(Int(samp_rate))_$(FFT.name)"
                            serialize(joinpath(root_fft, serial_fft), FFT)
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
    function sum_corrs(corrs::Array{CorrData,1})
        """ Implement sum function on array of corrs before stacking"""
        try
            sum_corr = Array{Float64, 2}(undef, size(corrs[1].corr)[1], length(corrs))
            for (ind, corr) in enumerate(corrs)
                sum_corr[:, ind] = corr.corr[:]
            end
            corr = deepcopy(corrs[1]) # keep metadata from first correlation
            corr.corr = sum_corr # update array of correlations
            return corr
        catch e 
            println("Error combining correlations: ", c)
            return nothing
        end
    end  
    function stack_auto(name::String, startdate::Date=startdate)
        try
            month = Dates.monthname(startdate) # get month for filename
            yr = Dates.year(startdate)
            # stack autocorrelations 

            CORROUT = expanduser(joinpath(params["rootdir"],"autocorrelations/$name/"))
            if !isdir(CORROUT)
                mkpath(CORROUT)
            end
            filename = joinpath(CORROUT,"$(yr)_$(month)_$name.h5") # get output filename

            # Get list of files to save
            autocorr_list = glob("AUTOCORR/$(yr)_$(params["month"])/$name*/*/*.jld2", params["rootdir"])
            if length(autocorr_list) == 0; return 1; end
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
                    try
                        comp_files = glob("AUTOCORR/$(yr)_$(params["month"])/$name*/$comp/*.jld2", params["rootdir"]) # uni-component files
                        if length(comp_files) == 0; continue; end # check if no files found - then continue

                        # load and clean autocorrelations
                        autocorrs = [load_corr(f, comp) for f in comp_files]
                        autocorrs = [c for c in autocorrs if typeof(c) == CorrData]

                        # Put correlations in single object, medianmute bad day correlations
                        corr_sum = sum_corrs(autocorrs) 
                        if isnothing(corr_sum); continue; end # check if nothing
                        cc_medianmute!(corr_sum, 3., false)

                        autocorr_mean = SeisNoise.stack(corr_sum, allstack=true, stacktype=mean)
                        autocorr_pws = SeisNoise.stack(corr_sum, allstack=true, stacktype=pws)
                        # save into file 
                        write(file, "$comp/linear", autocorr_mean.corr[:])
                        write(file, "$comp/pws", autocorr_pws.corr[:])
                    catch e 
                        println("Error stacking $name $comp: ", e)
                    end
                end
            end
        catch e
            println("Could not catch suberror. Cannot stack $name $startdate: ",e)
        end
    end
    function stack_corr(name::String, startdate::Date=startdate, prefix::String = "CORR_20HZ")
        try
            month = Dates.monthname(startdate) # get month for filename
            yr = Dates.year(startdate)
            # stack autocorrelations 
        
            CORROUT = expanduser(joinpath(params["rootdir"], "correlations/$name/"))
            if !isdir(CORROUT)
                mkpath(CORROUT)
            end
            filename = joinpath(CORROUT,"$(yr)_$(month)_$name.h5") # get output filename
            components = ["EE","EN","EZ", "NE", "NN","NZ", "ZE", "ZN", "ZZ"]
        
            receivers = glob("$prefix/$(yr)_$(params["month"])/$name/*/*", params["rootdir"])
            receivers = Set([join(split(x, ".")[end-2:end-1], ".") for x in receivers]) # just get reciever names
            corr_list = glob("$prefix/$(yr)_$(params["month"])/$name*/*/*", params["rootdir"])
            if length(corr_list) == 0; return 1; end

            println("Processing $(length(corr_list)) 20HZ corelations")
            C = load_corr(corr_list[1], convert(String, split(corr_list[1],"/")[end-1]))
            source_loc = C.loc
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
                        # add reciever location data if available in CAstations.csv
                        try
                            rec_network, rec_station = split(rec, ".")[1], split(rec, ".")[2]
                            rec_loc = LLE_geo(rec_network, rec_station, all_stations)
                            if !haskey(read(file), "$rec/meta")
                                write(file, "$rec/meta/lon", rec_loc.lon)
                                write(file, "$rec/meta/lat", rec_loc.lat)
                                write(file, "$rec/meta/el", rec_loc.el)
                                write(file, "$rec/meta/dist", get_dist(source_loc, rec_loc))
                                write(file, "$rec/meta/azi", get_azi(source_loc, rec_loc))
                                write(file, "$rec/meta/baz", get_baz(source_loc, rec_loc))
                            end
                        catch e 
                            println("Cannot get location information for $rec in CAstations.csv:", e)
                        end
                        for comp in components
                            try
                                # load correlations for this receiver by component 
                                files = glob("$prefix/$(yr)_$(params["month"])/$name/$comp/*$name..$rec.jld2", params["rootdir"])
                                if length(files) == 0; continue; end # check if no files found - then continue

                                corrs = [load_corr(f, comp) for f in files]
                                corrs = [c for c in corrs if typeof(c) == CorrData]

                                if comp == "EE"
                                    try
                                        if !haskey(read(file), "$rec/meta")
                                            write(file, "$rec/meta/lon", corrs[1].loc.lon)
                                            write(file, "$rec/meta/lat", corrs[1].loc.lat)
                                            write(file, "$rec/meta/el", corrs[1].loc.el)
                                            write(file, "$rec/meta/dist", get_dist(source_loc, corrs[1].loc))
                                            write(file, "$rec/meta/azi", get_azi(source_loc, corrs[1].loc))
                                            write(file, "$rec/meta/baz", get_baz(source_loc, corrs[1].loc))
                                        end
                                    catch e 
                                        println("Cannot add metadata:", e)
                                    end
                                end

                                # put daily correlations in single file, and remove noisy days
                                corr_sum = sum_corrs(corrs)
                                if isnothing(corr_sum); continue; end
                                cc_medianmute!(corr_sum, 3., false)

                                # implement various stacktypes
                                corr_mean = SeisNoise.stack(corr_sum, allstack=true, stacktype=mean)
                                corr_pws = SeisNoise.stack(corr_sum, allstack=true, stacktype=pws)
                                # save into file 
                                write(file, "$rec/$comp/linear", corr_mean.corr[:])
                                write(file, "$rec/$comp/pws", corr_pws.corr[:])
                            catch e 
                                println("Error stacking $name $comp: ", e)
                            end
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
        catch e 
            println("Failed to catch suberror while stacking $name $startdate. Cannot complete stack!")
        end
    end
    function stack_all(params::Dict=params)
        """ We stack each product, with robust error handling"""
        Sauto, S20, S1 = nothing, nothing, nothing
        try # AUTOCORRELATIONS
            autocorr_names = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("AUTOCORR/$(yr)_$(params["month"])/*", params["rootdir"])]
            Sauto = @elapsed robust_pmap(x -> stack_auto(x, startdate), autocorr_names)
        catch e
            println("Failed to catch error for autocorrelations stack for $(params["month"]).")
        end

        try # 20 HZ CORRELATIONS
            corr_names = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("CORR_20HZ/$(yr)_$(params["month"])/*", params["rootdir"])]
            S20 = @elapsed robust_pmap(x -> stack_corr(x, startdate), corr_names)
        catch e 
            println("Failed to catch error for 20 Hz correlations stacking for $(params["month"]).")
        end
        
        try # 1 HZ CORRELATIONS
            corr_names_lf = [join(split(split(auto, "/")[end], ".")[1:2],".") for auto in glob("CORR_1HZ/$(yr)_$(params["month"])/*", params["rootdir"])]
            S1 = @elapsed robust_pmap(x -> stack_corr(x, startdate, "CORR_1HZ"), corr_names_lf)
        catch e 
            println("Failed to catch error for 20 Hz correlations stacking for $(params["month"]).")
        end

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
    function cc_medianmute!(C::CorrData, cc_medianmute_α::Float64 = 10.0, bool::Bool = true)
        C.corr, inds = cc_medianmute(C.corr, cc_medianmute_α)
        if bool
            C.t = remove_medianmute(C, inds)
        end
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
    # add response for CI network data
    function yyyyjjj2date(yearday::String)
        @assert occursin(r"[1-2][0-9][0-9][0-9][0-3][0-9][0-9]",yearday)
        yint = parse(Int,yearday[1:4])
        dint = parse(Int,yearday[5:end])
        @assert dint <= 366 "Input day must be less than or equal to 366"
        return DateTime(yint) + Day(dint-1)
    end
    function read_resp(file::String,XMLDIR::String)
        s = yyyyjjj2date(file[end-9:end-3])
        t = s + Day(1)
        s = Dates.format(s, "yyyy-mm-dd HH:MM:SS")
        t = Dates.format(t, "yyyy-mm-dd HH:MM:SS")
        net = basename(file)[1:2]
        sta = split(basename(file),"_")[1][3:end]
        instpath = joinpath(XMLDIR,net * '_' * sta * ".xml" )
        return read_meta("sxml",instpath,s=s,t=t)
    end
    function add_response(S::SeisData, file::String, XMLDIR::String)
        resp = read_resp(file, joinpath(XMLDIR, "FDSNstationXML/CI")) 
        found_responses = resp[findfirst(x -> x == S.id[1], resp.id)] # filter response
        S.resp[1] = found_responses.resp # add response
        S.loc[1] = found_responses.loc
        S.gain[1] = found_responses.gain
    end
    function s3_hurl(file::String, params::Dict, split_id::String = "correlations")
        """ Robust s3_put with file repathing """
        try # try on each object in case some are corrupted (although s3_put should handle this)
            # adjust filepathing for cloud storage (removes /scratch/... or whatever prefix is)
            s3_name = join([split_id, convert(String, split(file, split_id)[2])]) 
            s3_put(params["aws"], "seisbasin", s3_name, read(file), acl="bucket-owner-full-control") # transfer to bucket
        catch e 
            try
                s3_name = join([split_id, convert(String, split(file, split_id)[2])]) # adjust fpath
                s3_put(params["aws"], "seisbasin", s3_name, read(file)) # transfer to bucket without acl (sometimes this bugs)
                println("Error uploading $file ($split_id) to seisbasin: ", e)
            catch e
                println("Error uploading $file ($split_id) to seisbasin (second failure): ", e)
            end
        end
    end
end
function XML_download(aws,XMLDIR)
    if !isdir(XMLDIR)
        mkpath(XMLDIR)
    end
    req = collect(s3_list_objects(aws,"scedc-pds","FDSNstationXML/CI/"))
    xmlin = [r["Key"] for r in req]
    xmlout = joinpath.(XMLDIR,basename.(xmlin))
    ec2download(aws, "scedc-pds", xmlin, XMLDIR)
    return nothing
end
function get_scedc_files(dd::Date, aws)
    try
        ar_filelist = pmap(x -> s3query(aws, dd, enddate = dd, network="CI", channel=x),["BH?", "HH?"])
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
    catch e 
        println("Cannot get SCEDC File list for $dd: ", e)
        return nothing
    end
end
function correlate_big(dd::Date, startdate::Date = startdate, params::Dict = params)
    """ Wrapper Function: Computes autocorrelations and correlations for a specific day"""
    try # Last resort failsafe
        println("Starting correlations for $dd.")
        aws, rootdir = params["aws"], params["rootdir"]
        yr = Dates.year(dd)
        path = join([Dates.year(dd),lpad(Dates.dayofyear(dd),3,"0")],"_") # Yeilds "YEAR_JDY"
        #@eval @everywhere path = $path

        ###### GET SCEDC FILES #######
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
        if isnothing(scedc_files); return 1; end # not worth computing if no SCEDC files - this should never happen
        ##############################

        # filepaths for iris and ncedc data
        iris_path  = S3Path("s3://seisbasin/iris_waveforms/$yr/$path/", config=aws)
        ncedc_path = S3Path("s3://seisbasin/ncedc_waveforms/$yr/$path/", config=aws)

        # format file strings for download 
        iris_query, ncedc_query = convert.(String, readdir(iris_path)), convert.(String, readdir(ncedc_path))
        filelist_basin = vcat(joinpath.("iris_waveforms/$yr/$path/", iris_query), joinpath.("ncedc_waveforms/$yr/$path/", ncedc_query))

        # print pre-download sumamry
        println("There are $(length(filelist_basin)) node files and $(length(scedc_files)) SCEDC files available for $path.")

        # download scedc and seisbasin (iris/ncedc) data
        ec2download(aws, "scedc-pds", scedc_files, joinpath(rootdir,"data"))
        ec2download(aws, "seisbasin", filelist_basin, joinpath(rootdir,"data"))
        println("Download complete!")


        # preprocess - broadbands and seismometers are processed separately
        raw_waveforms = glob("data/*/$yr/$path/*", rootdir)# params["rootdir"])
        println("$(length(raw_waveforms)) waveforms available for processing on $dd.")

        T_b = @elapsed serial_names = robust_pmap(f -> preprocess(f, params), raw_waveforms)
        println(serial_names[1:3])
        println(readdir())
        println(glob("ffts/*", "/scratch"))
        println("Preprocessing completed in $T_b seconds. $(length(serial_names)) ffts computed for $dd.")

        serial_names = collect(Iterators.flatten(filter(x -> !isnothing(x), serial_names))) # filter NaN arrays
        serial_names = convert(Array{String, 1}, filter(x -> !isnothing(x), serial_names)) # filter NaN elements
        serial_names = joinpath.(joinpath(params["rootdir"], "ffts/"), serial_names)

        # autocorrelations    
        fft_list_100 = filter(name -> occursin("100_", name), serial_names) # find 100Hz data
        fft_list_100 = convert(Array{String, 1}, fft_list_100) # convert to string arrray
        println("FFT list 100 is $(length(fft_list_100)) stations long for $dd.") # print summary
        fft_100_stations = unique([join(split(split(elt, "_")[end],".")[1:2] ,".") for elt in fft_list_100]) # get stations to iterate over
        println(fft_100_stations[1:5]) # print sample (CHECK)

        pmap(x -> autocorrelate(x, fft_list_100, params), convert.(String, fft_100_stations)) # Run autocorrelations 


        fft_paths_20 = sort(filter(name -> occursin("20_", name), serial_names)) # alphabetic sort for all stations so that diagonal files are made correctly
        fft_paths_1 = sort(filter(name -> occursin("1_", name), serial_names)) # alphabetic sort for all stations so that diagonal files are made correctly
        println("There are $(length(fft_paths_20)) fft datas for the 20 HZ correlations on $dd.") # print 20 Hz summary
        # run 20 HZ correlations
        chunks_20HZ, off_chunk_names_20HZ = get_blocks(fft_paths_20, params); # get chunking to distribute to workers

        # correlate diagonal chunks
        T20D = @elapsed robust_pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_20HZ", true, params), chunks_20HZ)

        # correlate off-diagonal chunks
        T20O = @elapsed robust_pmap(chunk -> offdiag_chunks(chunk, "CORR_20HZ", true, params), off_chunk_names_20HZ) # run mega correlations


        # run 1 Hz correlations 
        chunks_1HZ, off_chunk_names_1HZ = get_blocks(fft_paths_1, params); # chunking

        # correlate diagonal chunks
        T1D = @elapsed robust_pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_1HZ", false, params), chunks_1HZ)

        # correlate off-diagonal chunks
        T1O = @elapsed robust_pmap(chunk -> offdiag_chunks(chunk, "CORR_1HZ", false, params), off_chunk_names_1HZ) # run mega correlations

        # Correlation summary
        println("All $(length(glob("CORR_*/$(yr)_$(params["month"])/*/*/$path*", params["rootdir"]))) Inter-station Correlations computed in $(T20D + T20O + T1D + T1O) seconds")
        try; rm("$(params["rootdir"])/data/continuous_waveforms/", recursive=true); catch e; println(e); end
        GC.gc()
    catch e 
        println("This is embarrasing... Failed to compute correlations for $dd. Error:", e)
    end
end


println("Begin processing correlations for: ", startdate, " to ",enddate)

# Download instrument responses for CI data
T_download_XML = @elapsed XML_download(aws, params["XMLDIR"])
println("Metadata downloaded in $T_download_XML seconds.")

# Make directory for saving ffts into 
if !isdir("/scratch/ffts/")
    mkpath("/scratch/ffts/")
end
# Big wrapper function for processing and correlating everything

Tcbig = @elapsed map(dd -> correlate_big(dd, startdate, params), days)
println("All correlations completed in $Tcbig seconds.")

# Stack each of the 100, 20, and 1 Hz data
stack_all()


# transfer to s3
@everywhere begin
    autos = glob("autocorrelations/*/$(yr)_$(params["month"])*.h5", params["rootdir"])
    bigcorrs = glob("correlations/*/$(yr)_$(params["month"])*.h5", params["rootdir"])
end 

# robust transfer to s3
robust_pmap(x -> s3_hurl(x, params, "autocorrelations"), autos)
robust_pmap(x -> s3_hurl(x, params, "correlations"), bigcorrs)

println("Finished Correlations and upload for: ", startdate, " to ",enddate) # summary 


# Cleanup 
#old_corrs = glob("CORR_*/$(yr)_$(params["month"])/*/*/*.jld2", params["rootdir"])
# single day correlations
rm(joinpath(params["rootdir"], "CORR_20HZ/$(yr)_$(params["month"])/"), recursive=true)
rm(joinpath(params["rootdir"], "CORR_1HZ/$(yr)_$(params["month"])/"), recursive=true)
rm(joinpath(params["rootdir"], "AUTOCORR/$(yr)_$(params["month"])/"), recursive=true)
# month stacks 
rm.(autos)
rm.(bigcorrs)
println("Cleanup Complete!")