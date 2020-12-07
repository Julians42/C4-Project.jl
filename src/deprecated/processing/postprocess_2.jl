T = @elapsed using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays, AWSS3, ColorSchemes, Plots.PlotUtils, HDF5

#coeffs
cc_step, cc_len = 3600, 3600
maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9
half_win, water_level = 30, 0.01

# Select plotting parameters
frequency_plots = [[0.1,0.5],[0.5,1.0],[0.1,1.]] #[[0.2,0.3],[0.3,0.4],[0.4,0.5]]  #[[0.5,1.],[1.,2.],[0.1,0.2]]
lw = 0.5 #Decrease line thickness by half from default

σs = 0.5:0.1:2
normal_x = -5:0.01:5
normal_y = [exp.(-normal_x.^2 / (2σ^2)) / (2π * σ^2) for σ in σs];
loadcolorscheme(:cm_maxamp,ColorSchemes.gist_heat.colors[end-50:-1:1], "maxamp color", "for waveform plot");

using Pkg 
ENV["GR"] = ""
Pkg.build("GR")

#Add procs to access multiple cores
addprocs()
@everywhere using SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWSCore, StructArrays, AWSS3, Statistics, JLD2, Glob, HDF5
#@everywhere begin # helper functions for safe correlation download
 
"""
ec2download(aws,bucket,filelist,OUTDIR)
Download files using pmap from S3 to EC2.
# Arguments
- `aws::AWSConfig`: AWSConfig configuration dictionary
- `bucket::String`: S3 bucket to download from.
- `filelist::Array{String}`: Filepaths to download in `bucket`.
- `OUTDIR::String`: The output directory on EC2 instance.
# Keywords
- `v::Int=0`: Verbosity level. Set v = 1 for download progress.
- `XML::Bool=false`: Download StationXML files for request. Downloads StationXML to
    `joinpath(OUTDIR,"XML")`.
"""
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
    s = read(`du -k $dir`, String)
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
#end

# plotting functions
function comp_corrs(corrs::Array{CorrData,1}, comp::String)
    corrs_comp = Array{CorrData,1}(undef,0)
    for (index, corr) in enumerate(corrs)
        if corr.comp == comp
            push!(corrs_comp, corr)
        end
    end
    return corrs_comp
end
#Returns B40XX node data only
function filter_nodes(corrs::Array{CorrData,1}, corr_group)
    corrs_nodes = Array{CorrData,1}(undef,0)
    for (index, corr) in enumerate(corrs)
        if occursin("NO.$corr_group", corr.name)
            push!(corrs_nodes, corr)
        end
    end
    return corrs_nodes
end
function vert_plot(corrs_passed::Array{CorrData,1}, sta::String, comp_iter::String, line::String, stack_type::String)
    # Add filtered correlations to appropriate plots
    T = collect(-corrs_passed[1].maxlag:1/corrs_passed[1].fs:corrs_passed[1].maxlag)
    y_max = length(corrs_passed)*5+9
    for j in 1:length(frequency_plots)
        #Select frequency
        fmin = frequency_plots[j][1]
        fmax = frequency_plots[j][2]
        # Define plots 
        ZZ_plot = plot(xlims = (0,100), ylims = (0,y_max), 
                        yticks=[],xlabel = "Time (s)", ylabel = "South (Bottom) to North (Top)", 
                        xtickfontsize=5,ytickfontsize=5,fontsize=5, xguidefontsize = 10, yguidefontsize = 10,
                        legendfontsize = 15)
        
        # Adjust scaling by comparison by first corr
        Cstack = bandpass(SeisNoise.stack(corrs_passed[1]), fmin, fmax)
        id_mid = round(Int, size(Cstack.corr, 1))
        scaled = 0.15*maximum(broadcast(abs, shorten(Cstack,200.).corr))
        max_amplitudes = Array{Float64, 1}(undef,0)
        plot_process = []
        # Stack over days and add to plot
        for i in 1:length(corrs_passed)
            Cstack = bandpass(SeisNoise.stack(corrs_passed[i]), fmin, fmax)
            push!(max_amplitudes, maximum(shorten(Cstack,200.).corr/scaled))
            push!(plot_process, Cstack.corr/scaled .+5*i)
        end
        title!("$(sta) S$line $(comp_iter) using $stack_type", fontsize=5)
        plot!(ZZ_plot, -T, plot_process, color=:cm_maxamp, colorbar_title="Normalized Maximum Amplitude", 
            line_z=max_amplitudes', fmt = :png, linewidth = lw, reuse = false, legend = false)
        plot!(size=(250,400),dpi=1000)
        filepath = "~/stack_plots/$stack_type/S$(line)_$(sta)_$(comp_iter)_$(fmin)to$(fmax).png"
        # ensure filepath is valid 
        DIR = dirname(filepath)
        if !isdir(DIR)
            mkpath(DIR)
        end
        png(ZZ_plot,"stack_plots/$stack_type/S$(line)_$(sta)_$(comp_iter)_$(fmin)to$(fmax).png")
    end
end

@everywhere begin
    using SCEDC, AWSCore, Dates, DataFrames, AWSS3
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    startdate = "2019-01-01" # Select Start Date
    enddate = "2019-12-31" # Select End Date
    days = Date(startdate):Day(1):Date(enddate)
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
job_name = "linear/2019"
@eval @everywhere job_name = $job_name
############################# Get files for unique station pair ####################

to_download = unique(Iterators.flatten([DataFrame(CSV.File(file)).paths for file in month_fnames]))

safe_download(aws, bucket2, to_download, "~/")
println("Download Complete")

#temporary


# @everywhere begin
#     function s3_file_map(aws::AWSConfig,seisbucket::String,filein::String,fileout::String)
#         try
#             s3_get_file(aws, seisbucket, filein, fileout)
#             println("Downloading file: $filein       \r")
#             return nothing
#         catch e
#             println("Unable to download file $filein")
#         end
#     end
# end
# temp_download = [x for x in to_download if occursin("February", x)==true]
# pmap(x-> s3_file_map(aws, bucket2, x, x), temp_download)



# safe_download(aws, bucket2, temp_download, "~/")






function csv_merge(large_index = Array{String, 1})
    """Returns dataframe with unique station pairs and csvs containing that station pair"""
    all_pair_paths = DataFrame(pair =String[], files = Array[])
    unq_pairs = unique(Iterators.flatten([DataFrame(CSV.File(file)).Files for file in large_index]))
    for pair in unq_pairs
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
        push!(all_pair_paths, [pair, pair_paths])
    end
    return all_pair_paths
end

pair_paths_df = csv_merge(month_fnames)

function csv_merge(large_index = Array{String, 1})
    """Returns dataframe with unique station pairs and csvs containing that station pair"""
    all_pair_paths = DataFrame(source = String[], pair =String[], files = Array[])
    unq_pairs = unique(Iterators.flatten([DataFrame(CSV.File(file)).Files for file in large_index]))
    stations = [convert(String, strip(join(split(pair,".")[1:3],"."),'.')) for pair in unq_pairs]
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
        push!(all_pair_paths, [stations[ind], pair, pair_paths])
    end
    return all_pair_paths
end

pair_paths_df = csv_merge(month_fnames)





# components =["EE", "EN", "EZ", "NE", "NN", "NZ", "ZE", "ZN", "ZZ"]
# # load, sum, stack and rotate correlations
# function postprocess_corrs(ar_corr_large)
#     try
#         if any(occursin.(".GATR.", ar_corr_large)) == true # crappy GATR station data - should catch earlier
#             return nothing
#         end
#         comp_mean = Array{CorrData,1}(undef, 0) # 9 element arrays for each component 
#         comp_pws = Array{CorrData,1}(undef, 0) # across stacktypes
#         comp_robust = Array{CorrData,1}(undef, 0)
#         # iterate and stack across components
#         for key in components
#             # Get all correlations for particular component
#             corr_singles = Iterators.flatten([corr_load(corr_large, key) for corr_large in ar_corr_large])
#             #filter!(x -> x! = nothing, corr_singles)
#             # stack and append to array
#             corr_mean = SeisNoise.stack(sum(corr_singles), allstack=true, stacktype=mean)
#             corr_pws = SeisNoise.stack(sum(corr_singles), allstack=true, stacktype=pws)
#             corr_robust = SeisNoise.stack(sum(corr_singles), allstack=true, stacktype=robuststack)
#             push!(comp_mean, corr_mean)
#             push!(comp_pws, corr_pws)
#             push!(comp_robust, corr_robust)
#         end
#         # Rotate for each stacktype based on locations 
#         for stack in [comp_mean, comp_pws, comp_robust]
#             #map(x -> add_corr_locations(x, df), stack)
#             SeisNoise.rotate!(stack, stack[1].azi, stack[1].baz) # rotate
#         # write to disk
#             #map(C-> save_processed(C, "processed/$job_name"), stack)
#         end
#         return [comp_mean, comp_pws, comp_robust] # array of arrays containing all 9 components by stack type. 
#     catch e
#         println(e)
#         return nothing
#     end
# end

# all_stacktypes = postprocess_corrs(pair_paths_df.files[1])
# #["RR","RT","RZ","TR","TT","TZ","ZR","ZT","ZZ"]
# # save files
# function save_hdf5(AC::Array{Array{CorrData,1}}, name::String, receiver::String)
#     println(name, receiver)
#     C = AC[1][1] # get metadata from first correlation
#     T = u2d(C.t[1]) # get starttime
#     D=Dict([ ("corr_type",C.corr_type), ("cc_len",C.cc_len),("cc_step",C.cc_len),
#     ("whitened",C.whitened), ("time_norm",C.time_norm), ("notes",C.notes),
#     ("maxlag",C.maxlag),("starttime",T)])
#     println(name)
#     # write metadata and add stacktypes
#     h5open(name,"cw") do file # read and write
#         if !haskey(read(file),"meta") # if metadata isn't already added, add it
#             write(file, "meta/corr_type", C.corr_type)
#             write(file, "meta/cc_len", C.cc_len)
#             write(file, "meta/cc_step", C.cc_step)
#             write(file, "meta/whitened", C.whitened)
#             write(file, "meta/time_norm", C.time_norm)
#             write(file, "meta/notes", C.notes)
#             write(file, "meta/maxlag", C.maxlag)
#             write(file, "meta/starttime", Dates.format(T, "yyyy-mm-dd HH:MM:SS"))
#         end
#         println("meta finished")
#         if !haskey(read(file), receiver)
#             println("test",receiver)
#             g=g_create(file,receiver)
#             # g[C.comp]=C.corr[:]
#             attrs(g)["Description"] = ["linear", "pws", "robust"]
#             #h=g_create(g,stacktype)
#             for i in 1:9
#                 comp = AC[1][i].comp
#                 c = g_create(g, comp)
#                 for (ind, stackt) in enumerate(["linear","pws","robust"])
#                     print(stackt)
#                     write(file, "$receiver/$comp/$stackt", AC[ind][i].corr[:])
#                 end
#             end
#         end
#     end
# end
# function save_processed(AC::Array{Array{CorrData,1}}, CORROUT::String)
#     # check if CORRDIR exists
#     CORROUT = expanduser(CORROUT)
#     if isdir(CORROUT) == false
#         mkpath(CORROUT)
#     end

#     #name is YEAR_PAIRNAME.h5 - we store all components and all stacktypes
#     name = AC[1][1].name
#     yr = Dates.year(Date(AC[1][1].id))
#     p_name = strip(join(split(name,".")[1:3],"."),'.')
#     receiver_name = strip(join(split(name,".")[end-3: end-1],"."),'.')
#     println(receiver_name)
#     name_source = join([yr,p_name],"_")

#     # create JLD2 file and save correlation
#     filename = joinpath(CORROUT,"$(name_source).h5")
#     save_hdf5(AC, filename, convert(String,receiver_name))
# end


# save_processed(all_stacktypes, "processed")

# indices = []

# h5open("processed/2019_CI.CHN.h5", "w") do file
#     g = g_create(file, "CI.ADO")
#     attrs(g)["Description"] = ["linear", "pws", "robust"]
#     println(g)
#     #h=g_create(g,stacktype)
#     for i in 1:9
#         comp = all_stacktypes[1][i].comp
#         println(comp)
#         c = g_create(g, comp)
#         for (ind, stackt) in enumerate(["linear","pws","robust"])
#             #stackt = "linear"

#             #g2 = g_create(c, stackt)
#             println("CI.ADO/$comp/$stackt")
#             write(file, "CI.ADO/$comp/$stackt", all_stacktypes[ind][i].corr[:])
#             push!(indices, [ind,i])
#             println(stackt, "done1")
#         end
#     end
# end

# fid = h5open("processed/2019_CI.CHN.h5","r")
# A = read(fid)
# collect(keys(A[collect(keys(A))[1]]))[1]
# A[collect(keys(A))[1]][1]




@everywhere begin 
    #using DataFrames, JLD2, SeisIO, SeisNoise, HDF5
    # get lat, lon, and elevation (LLE) for station 
    # gets array of corrs for particular file
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
    function save_hdf5(AC::Array{Array{CorrData,1}}, name::String, receiver::String)
        #println(name, receiver)
        C = AC[1][1] # get metadata from first correlation
        T = u2d(C.t[1]) # get starttime
        D=Dict([ ("corr_type",C.corr_type), ("cc_len",C.cc_len),("cc_step",C.cc_len),
        ("whitened",C.whitened), ("time_norm",C.time_norm), ("notes",C.notes),
        ("maxlag",C.maxlag),("starttime",T)])
        #println(name)
        # write metadata and add stacktypes
        h5open(name,"cw") do file # read and write
            if !haskey(read(file),"meta") # if metadata isn't already added, add it
                write(file, "meta/corr_type", C.corr_type)
                write(file, "meta/cc_len", C.cc_len)
                write(file, "meta/cc_step", C.cc_step)
                write(file, "meta/whitened", C.whitened)
                write(file, "meta/time_norm", C.time_norm)
                write(file, "meta/notes", C.notes)
                write(file, "meta/maxlag", C.maxlag)
                write(file, "meta/starttime", Dates.format(T, "yyyy-mm-dd HH:MM:SS"))
            end
            #println("meta finished")
            if !haskey(read(file), receiver)
                #println("test",receiver)
                g=g_create(file,receiver)
                # g[C.comp]=C.corr[:]
                attrs(g)["Description"] = ["linear", "pws", "robust"]
                #h=g_create(g,stacktype)
                for i in 1:9
                    comp = AC[1][i].comp
                    c = g_create(g, comp)
                    for (ind, stackt) in enumerate(["linear","pws","robust"])
                        #print(stackt)
                        write(file, "$receiver/$comp/$stackt", AC[ind][i].corr[:])
                    end
                end
            end
        end
    end
    function save_processed(AC::Array{Array{CorrData,1}}, CORROUT::String)
        # check if CORRDIR exists
        CORROUT = expanduser(CORROUT)
        if isdir(CORROUT) == false
            mkpath(CORROUT)
        end
    
        #name is YEAR_PAIRNAME.h5 - we store all components and all stacktypes
        name = AC[1][1].name
        yr = Dates.year(Date(AC[1][1].id))
        p_name = strip(join(split(name,".")[1:3],"."),'.')
        receiver_name = strip(join(split(name,".")[end-3: end-1],"."),'.')
        println(receiver_name)
        name_source = join([yr,p_name],"_")
    
        # create JLD2 file and save correlation
        filename = joinpath(CORROUT,"$(name_source).h5")
        save_hdf5(AC, filename, convert(String,receiver_name))
    end
    components =["EE", "EN", "EZ", "NE", "NN", "NZ", "ZE", "ZN", "ZZ"]
    # load, sum, stack and rotate correlations
    function postprocess_corrs(ar_corr_large)
        try
            if any(occursin.(".GATR.", ar_corr_large)) # crappy GATR station data - should catch earlier
                return nothing
            end
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
                push!(comp_mean, corr_mean)
                push!(comp_pws, corr_pws)
                push!(comp_robust, corr_robust)
            end
            # # Location patch - only add to 6 stacked corrs
            # for stack in [comp_mean, comp_pws, comp_robust]
            #     #map(x -> add_corr_locations(x, df), stack)
            #     SeisNoise.rotate!(stack, stack[1].azi, stack[1].baz) # rotate
            # # write to disk
            #     map(C-> save_processed(C, "processed/$job_name"), stack)
            # end
            save_processed([comp_mean, comp_pws, comp_robust], "processed")
            #return nothing #[comp_mean, comp_pws, comp_robust] # to check work - probably don't need to return unless plotting
        catch e
            println(e)
            return nothing
        end
    end
end

# process raw correlations 
T_process = @elapsed map(x -> postprocess_corrs(x), pair_paths_df.files)
println("Done")

#Transfer stacked correlations
stacked_files = joinpath.("processed/$job_name", readdir("processed/$job_name"))
@elapsed pmap(x -> s3_put(aws, bucket2, x, read(x)), stacked_files)
println("Done")













function postprocess_corrs(ar_corr_large)
    # let's write this function so that everything gets saved into one file
    # get/make the output dir (necessary?)
    CORROUT = expanduser("processed/$job_name")
    if isdir(CORROUT) == false
        mkpath(CORROUT)
    end
    # get the year
    yr = Dates.year(Date(ar_corr_large[1].id))
    # get the virtual source name
    s_name = strip(join(deleteat!(split(ar_corr_large[1].name,"."),[4,8]),"."),'.')
    # make list of receiver names
    list_of_receivers=#XX
    # name = join([yr,p_name,C.comp],"_")
    # get output filename
    filename = joinpath(CORROUT,"$(sname).h5")
    T = u2d(ar_corr_large[1].t[1]) # get starttime
    D=Dict([ ("corr_type",C.corr_type), ("cc_len",C.cc_len),("cc_step",C.cc_len),
    ("whitened",C.whitened), ("time_norm",C.time_norm), ("notes",C.notes),("dist",C.dist),
    ("azi",C.azi),("baz",C.baz),("maxlag",C.maxlag),("starttime",T)])
    # write metadata and add stacktypes
    h5open(filename,"w") do file
    m = g_create(file, "meta") #?
    m["meta"] = D #?
    #  g=g_create(file,"stack")
    #  # g[C.comp]=C.corr[:]
    #  # attrs(g)["Description"] = "linear stack"
    #  h=g_create(g,stacktype)
    #  h[C.comp]=C.corr[:]
    #  end
        for (ir,sta) in enumerate(list_of_receivers)
            iik = # list of indices of all files for a single receiver
            g=g_create(file,sta) # create group with receiver name
            try
                if any(occursin.(".GATR.", ar_corr_large[iik])) == true # crappy GATR station data - should catch earlier
                    return nothing
                end
                comp_mean = Array{CorrData,1}(undef, 0)
                comp_pws = Array{CorrData,1}(undef, 0)
                comp_robust = Array{CorrData,1}(undef, 0)
                # iterate and stack across components
                for key in components
                    # Get all correlations for particular component
                    corr_singles = Iterators.flatten([corr_load(corr_large, key) for corr_large in ar_corr_large[iik]])
                    #filter!(x -> x! = nothing, corr_singles)
                    # stack and append to array
                    corr_mean = SeisNoise.stack(sum(corr_singles), allstack=true, stacktype=mean)
                    corr_pws = SeisNoise.stack(sum(corr_singles), allstack=true, stacktype=pws)
                    corr_robust = SeisNoise.stack(sum(corr_singles), allstack=true, stacktype=robuststack)
                    # push!(comp_mean, corr_mean)
                    # push!(comp_pws, corr_pws)
                    # push!(comp_robust, corr_robust)
                    # save to disk here
                    g[comp_mean.comp]=corr_mean.corr[:]
                    g[comp_mean.comp]=corr_pws.corr[:]
                    g[comp_mean.comp]=corr_robust.corr[:]
                end
                # Location patch + write - only add to 6 stacked corrs
                # for stack in [comp_mean, comp_pws, comp_robust]
                    # map(x -> add_corr_locations(x, df), stack)
                    # SeisNoise.rotate!(stack, stack[1].azi, stack[1].baz) # rotate
                # write to disk
                    # map(C-> save_processed(C, "processed/$job_name"), stack)
                # end
                return  [comp_mean, comp_pws, comp_robust] # to check work - probably don't need to return unless plotting
            catch e
                println(e)
                return nothing
            end # end of trying
        end # end of loop over receivers
    end # file close
end # end of function







# Plotting functionality not yet available - coming soon 
################################### Make and Transfer Node Plots ########################################
println("Plotting nodal data now!")
types = ["mean", "robust", "pws"]
nodes = ["B4", "G1", "G2"]
components = ["ZZ", "TT","RR"]
for (ind, stacks) in enumerate([processed_mean, processed_robust])#, processed_pws])
    stack_type = types[ind]
    for (j,k) in Iterators.product(1:length(nodes), 1:length(components))
        node, component = nodes[j], components[k]
        filtered_nodes = filter_nodes(comp_corrs(stacks, component), node)
        vert_plot(filtered_nodes, source_station, component, node, stack_type)
    end
end

# Transfer node plots
stacked_plots = glob("*/*$(source_station)*.png", "stack_plots") # Transfer plots from this station only
@elapsed pmap(x -> s3_put(aws, bucket2, "stack_plots", read(x)), stacked_plots)