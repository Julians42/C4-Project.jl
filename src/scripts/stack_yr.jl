T = @elapsed using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays, AWSS3, ColorSchemes, Plots.PlotUtils, HDF5

#Add procs to access multiple cores
addprocs()
@everywhere using SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWSCore, StructArrays, AWSS3, Statistics, JLD2, Glob, HDF5

#coeffs
cc_step, cc_len = 3600, 3600
maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9
half_win, water_level = 30, 0.01

@everywhere begin
    using SCEDC, AWSCore, Dates, DataFrames, AWSS3
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    startdate = "$(ARGS[1])-01-01" # Select Start Date
    enddate = "$(ARGS[1])-12-31" # Select End Date
    days = Date(startdate):Day(1):Date(enddate)
end


function safe_download(aws::AWSConfig,seisbucket::String,file_list::Array{String},OUTDIR::String;
                        v::Int=0,XML::Bool=false)

    # check being run on AWS
    tstart = now()
    !localhost_is_ec2() && @warn("Running locally. Run on EC2 for maximum performance.")

    println("Starting Download...      $(now())")
    println("Using $(nworkers()) cores...")


    # create output files
    OUTDIR = expanduser(OUTDIR)
    outfiles = [joinpath(OUTDIR,f) for f in file_list]
    filedir = unique([dirname(f) for f in outfiles])
    for ii = 1:length(filedir)
        if !isdir(filedir[ii])
            mkpath(filedir[ii])
        end
    end

    # get XMLmesse
    if XML
        XMLDIR = joinpath(OUTDIR,"XML")
        getXML(aws,seisbucket,file_list,XMLDIR,v=v)
    end

    # send outfiles everywhere
    @eval @everywhere outfiles=$outfiles
    # do transfer to ec2
    startsize = diskusage(OUTDIR)
    if v > 0
        pmap(
            s3_file_map,
            fill(aws,length(outfiles)),
            fill(seisbucket,length(outfiles)),
            file_list,
            outfiles
        )
    else
        pmap(
            s3_get_file,
            fill(aws,length(outfiles)),
            fill(seisbucket,length(outfiles)),
            file_list,
            outfiles,
        )
    end

    println("Download Complete!        $(now())          ")
    tend = now()
    # check data in directory
    endsize = diskusage(OUTDIR)
    downloadsize = endsize - startsize
    downloadseconds = (tend - tstart).value / 1000
    println("Download took $(Dates.canonicalize(Dates.CompoundPeriod(tend - tstart)))")
    println("Download size $(formatbytes(downloadsize))")
    println("Download rate $(formatbytes(Int(round(downloadsize / downloadseconds))))/s")
    return nothing
end

function s3_file_map(aws::AWSConfig,seisbucket::String,filein::String,fileout::String)
    try
        s3_get_file(aws, seisbucket, filein, fileout)
        println("Downloading file: $filein       \r")
        return nothing
    catch e
        println("Unable to download file $filein")
    end
end

function s3_get_seed(
    aws::AWSConfig,seisbucket::String,
    filein::String,
    demean::Bool,
    detrend::Bool,
    msr::Bool,
    prune::Bool,
    rr::Bool,
    taper::Bool,
    ungap::Bool,
    unscale::Bool,
    resample::Bool,
    fs::Real,)
    f = s3_get(aws, seisbucket, filein)
    S = parseseed(f)

    # remove empty channels
    if prune == true
        prune!(S)
    end

    # Get list of channels with sane instrument codes
    CC = get_seis_channels(S)

    if msr == true
        @warn("Getting response not implemented yet.")
    end

    # unscale
    if unscale == true
        unscale!(S, chans=CC)
    end

    # Demean
    if demean == true
    demean!(S, chans=CC)
    end

    # Taper
    if taper == true
    taper!(S, chans=CC)
    end

    # Ungap
    if ungap == true
    ungap!(S, chans=CC)
    end

    # resample data
    if resample == true && fs != 0
        resample!(S, chans=CC, fs=fs)
    end

    # Remove response
    # need to implement attaching response
    if rr == true
    @warn("Removing response not implemented yet.")
    end

    return S
end

function formatbytes(bytes::Real, digits::Int=1)
    units = ["B", "KB", "MB", "GB", "TB","PB"]
    bytes = max(bytes,0)
    pow = Int(floor((bytes > 0 ? log(bytes) : 0) / log(1024)))
    pow = min(pow,length(units))
    powind = pow < length(units) ? pow + 1 : pow
    return string(round(bytes / 1024 ^ pow,digits=digits)) * units[powind]
end

function diskusage(dir)
    s = read(`du -k $(dir)`, String)
    kb = parse(Int, split(s)[1])
    return 1024 * kb
end

"""
parseseed(f)
Convert uint8 data to SeisData.
"""
function parseseed(f::AbstractArray)
    S = SeisData()
    SeisIO.SEED.parsemseed!(S,IOBuffer(f),SeisIO.KW.nx_new,SeisIO.KW.nx_add,false,0)
    return S
end


############################# Index Download #####################################

# Get list of CSVs to transfer corr_index/2017/2017_026_correlation_index.csv
month_index_unq = unique([(Dates.year(d), Dates.monthname(d)) for d in days])
month_fnames = ["month_index/$(ind[1])/$(ind[2]).csv" for ind in month_index_unq]

#Download CSVs
@eval @everywhere month_fnames = $month_fnames # Retain SCEDC download functionality, parallel download not actually needed
safe_download(aws, bucket2, month_fnames, "~/") # Retain seisbasin filepathing

stations = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV"
            ,"VINE", "ROPE", "ARNO", "LUCI", "ROUF", "KUZD", "ALLI", "CHN", "USB", "Q0048"]
df = DataFrame(CSV.File("files/full_socal.csv"))
@eval @everywhere df = $df
job_name = "linear/2018"
@eval @everywhere job_name = $job_name
############################# Get files for unique station pair ####################

to_download = unique(Iterators.flatten([DataFrame(CSV.File(file)).paths for file in month_fnames]))

safe_download(aws, bucket2, to_download, "~/")
println("Download Complete")



function csv_merge(large_index = Array{String, 1})
    """Returns dataframe with unique station pairs and csvs containing that station pair"""
    all_pair_paths = DataFrame(source = String[], receiver = String[], pair =String[], files = Array[])
    unq_pairs = unique(Iterators.flatten([DataFrame(CSV.File(file)).Files for file in large_index]))
    stations = [convert(String, strip(join(split(pair,".")[1:3],"."),'.')) for pair in unq_pairs]
    stations_recievers = [convert(String, strip(join(split(pair,".")[4:end],"."),'.')) for pair in unq_pairs]
    for (ind, pair) in enumerate(unq_pairs)
        pair_paths = Array{String, 1}(undef, 0)
        for csv in large_index
            try
                df = DataFrame(CSV.File(csv))
                path = df[(findall(x -> x==pair, df.Files)),:].paths[1]
                push!(pair_paths, path)
            catch 
                # That CSV doesn't have that pair - not really a problem!
            end
        end
        #print(pair_paths)
        push!(all_pair_paths, [stations[ind], stations_recievers[ind], pair, pair_paths])
    end
    return all_pair_paths
end

pair_paths_df = csv_merge(month_fnames)

sources = unique(pair_paths_df.source)


components =["EE", "EN", "EZ", "NE", "NN", "NZ", "ZE", "ZN", "ZZ"]
@eval @everywhere components = $components
@eval @everywhere pair_paths_df = $pair_paths_df
@eval @everywhere sources = $sources
@everywhere begin
    function corr_load(corr_large, key)
        try
            jld = jldopen(corr_large, "r")
            files = keys(jld[key])
            corrs_comp = [jld["$key/$file"] for file in files]
            close(jld)
            return corrs_comp
        catch e
            println(e)
            println("Likely $corr_large does not contain some components.")
        end
    end
    function LLE_geo(network, station, df)
        """ Find station matching location and return geoloc object"""
        try
            row = filter(row -> (row.network==network) & (row.station==station), df)[1,:]
            return GeoLoc(lat = float(row.latitude), lon = float(row.longitude), el = float(row.elevation))
        catch 
            return nothing # when station is not found in dataframe
        end
    end
    function postprocess_corrs(source, pair_df, stations) 
        """ Stacks monthly jld2 files into h5, saving comprehensive metadata
            - Arguments: 
                        source - string for source station
                        pair_df - dataframe containing file pathing information
                        stations - dataframe with location information
            - Returns (nothing):
                        h5 file for source station with stacked yearly correlations
                        for each receiver station 
                        contains receiver and source metadata
        """
        # filter receivers, throwing out bad GATR station
        df_source = filter(row -> (row.source ==source) & (!occursin("GATR", row.receiver)), pair_df) # get receivers
        list_of_receivers = df_source.receiver  # list of receiver names

        # get example correlation to extract metadata
        C = corr_load(df_source.files[1][1],"ZZ")[1]
        # file naming YEAR_NET.STA.h5
        yr, p_name = Dates.year(Date(C.id)), strip(join(split(C.name,".")[1:3],"."),'.')
        name_source = join([yr,p_name],"_")
        CORROUT = expanduser("source_processed/$yr/")
        if !isdir(CORROUT)
            mkpath(CORROUT)
        end
        filename = joinpath(CORROUT,"$(name_source).h5") # get output filename

        # get source location information
        src_net, src_sta = convert(String,split(source,".")[1]), convert(String,split(source,".")[2])
        src_geo = LLE_geo(src_net, src_sta, stations)
    
        # write metadata and add stacktypes
        h5open(filename,"cw") do file
            if !haskey(read(file),"meta") # if metadata isn't already added, add it
                write(file, "meta/corr_type", C.corr_type)
                write(file, "meta/cc_len", C.cc_len)
                write(file, "meta/cc_step", C.cc_step)
                write(file, "meta/whitened", C.whitened)
                write(file, "meta/time_norm", C.time_norm)
                write(file, "meta/notes", C.notes)
                write(file, "meta/maxlag", C.maxlag)
                write(file, "meta/starttime", Dates.format(u2d(C.t[1]), "yyyy-mm-dd HH:MM:SS"))
                write(file, "meta/samp_freq", fs)
                if !isnothing(src_geo) # add source location
                    write(file, "meta/lat", src_geo.lat)
                    write(file, "meta/lon", src_geo.lon)
                    write(file, "meta/el", src_geo.el)
                    write(file, "meta/dep", src_geo.dep)
                    write(file, "meta/az", src_geo.az)
                    write(file, "meta/inc", src_geo.inc)
                end
            end
            count=0
            for (ind, receiver) in enumerate(list_of_receivers)
                try
                    if !haskey(read(file), receiver)
                        # add location metadata
                        rec_net, rec_sta = split(receiver,".")[1], split(receiver,".")[2]
                        rec_geo = LLE_geo(rec_net, rec_sta, stations)
                        if !isnothing(src_geo) & !isnothing(rec_geo)
                            dist, azi, baz = get_dist(src_geo, rec_geo), get_azi(src_geo, rec_geo), get_baz(src_geo, rec_geo)
                            write(file, "$receiver/meta/dist", dist)
                            write(file, "$receiver/meta/azi", azi)
                            write(file, "$receiver/meta/baz", baz)
                            write(file, "$receiver/meta/lat", rec_geo.lat)
                            write(file, "$receiver/meta/lon", rec_geo.lon)
                            write(file, "$receiver/meta/el", rec_geo.el)
                            write(file, "$receiver/meta/dep", rec_geo.dep)
                            write(file, "$receiver/meta/az", rec_geo.az)
                            write(file, "$receiver/meta/inc", rec_geo.inc)
                        end
                        # stack and process reciever correlations
                        ar_corr_large = filter(row -> row.receiver == receiver, df_source).files[1]
                        try
                            comp_mean = Array{CorrData,1}(undef, 0)
                            comp_pws = Array{CorrData,1}(undef, 0)
                            comp_robust = Array{CorrData,1}(undef, 0)
                            # iterate and stack across components
                            for key in components
                                # Get all correlations for particular component
                                corr_singles = Iterators.flatten([corr_load(corr_large, key) for corr_large in ar_corr_large])
                                #filter!(x -> x! = nothing, corr_singles)
                                # stack and append to array
                                corr_mean = SeisNoise.stack(sum(corr_singles), allstack=true, stacktype=mean)
                                corr_pws = SeisNoise.stack(sum(corr_singles), allstack=true, stacktype=pws)
                                corr_robust = SeisNoise.stack(sum(corr_singles), allstack=true, stacktype=robuststack)
                                #println(length(corr_mean))
                                # save to disk here
                                write(file, "$receiver/$key/linear", corr_mean.corr[:])
                                write(file, "$receiver/$key/pws", corr_pws.corr[:])
                                write(file, "$receiver/$key/robust", corr_robust.corr[:])
                            end
                        catch e
                            println(e)
                            return nothing
                        end # end of trying :)
                    end
                    count +=1 # track progress
                    println("Iteration $count of $(length(list_of_receivers)) complete: processed receiver $receiver")
                catch e
                    println(e)
                end
            end
        end # file close
    end # end of function
end


T = @elapsed pmap(x-> postprocess_corrs(x, pair_paths_df,df), sources)
println("Year Stack Completed in $T seconds")

#postprocess_corrs(pair_paths_df.source[1],pair_paths_df[1:2,:],df)

# transfer to seisbasin
stacked_files = glob("source_processed/*/*.h5")
Transfer = @elapsed pmap(x -> s3_put(aws, "seisbasin", x, read(x)), stacked_files)
println("Transfer Completed in $Transfer seconds.")