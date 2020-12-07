T = @elapsed using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays, AWSS3

#coeffs
cc_step, cc_len = 3600, 3600
maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9
half_win, water_level = 30, 0.01


using Pkg 
ENV["GR"] = ""
Pkg.build("GR")

#Add procs to access multiple cores
addprocs()


@everywhere using SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWSCore, StructArrays, AWSS3, Statistics, JLD2
@everywhere begin 
    function preprocess2(ar::SeisData, fs::Float64, freqmin::Float64, freqmax::Float64, cc_step::Int64, cc_len::Int64, half_win::Int64, water_level::Float64)
        try
            #println(ar.name)
            process_raw!(ar, fs) # downsample to 20 Hz 
            R = RawData(ar,cc_len,cc_step)
            SeisNoise.detrend!(R)
            SeisNoise.taper!(R)
            bandpass!(R,freqmin,freqmax,zerophase=true)
            FFT = compute_fft(R) # Compute the Fourier Transform
            coherence!(FFT,half_win, water_level)
            return(FFT)
        catch  a # Handle error catching
            println(a)
        end
    end
    function cc_medianmute!(C::CorrData, cc_medianmute_α::Float64 = 10.0)
        C.corr, inds = cc_medianmute(C.corr, cc_medianmute_α)
        C.t = remove_medianmute(C, inds)
        return nothing
    end
    function cc_medianmute(A::AbstractArray, cc_medianmute_α::Float64 = 10.0)
        #1. compute median of maximum amplitude of all corrs
        T, N = size(A)
        cc_maxamp = vec(maximum(abs.(A), dims=1))
        cc_medianmax = median(cc_maxamp)
        inds = findall(x-> x <= cc_medianmute_α*cc_medianmax,cc_maxamp)
        #NOTE: you cannot unbind entire array, so remove_nanandzerocol! is not used here.
        return A[:, inds], inds
    end
    remove_medianmute(C::CorrData, inds) = (return C.t[inds])
    function correlate5(ffts::Array{FFTData,1}, maxlag::Float64, dir::String)
        C = correlate(ffts[1],ffts[2],maxlag) 
        #clean_up!(C,freqmin,freqmax)
        cc_medianmute!(C, 10.)  #filter extraneous indices
        stack!(C)
        save_corr(C, dir) 
    end
    function correlatedir(ffts::Array{FFTData,1}, maxlag::Float64)
        # for j in 1:length(recievers) #i or 1 -- IMPORTANT
        C = correlate(ffts[1],ffts[2],maxlag)
        #Add second to filter only good indices to stack
        cc_medianmute!(C, 10.)
        stack!(C)
        pair = strip(join(deleteat!(split(C.name,"."),[4,8]),"."),'.')
        comp = C.comp
        starttime = u2d(C.t[1])
        save_corr2(C,"CORR/$pair/$comp")
        #try save_corr2(C, "CORR/$pair/$comp"); catch; println("$(C.name) alreadly correlated"); end
    end
    function save_corr2(C::CorrData, CORROUT::String)
        # check if CORRDIR exists
        CORROUT = expanduser(CORROUT)
        if isdir(CORROUT) == false
            mkpath(CORROUT)
        end

        #name is YEAR_JDY_PAIRNAME
        yr,j_day = Dates.year(Date(C.id)), lpad(Dates.dayofyear(Date(C.id)),3,"0")
        p_name = strip(join(deleteat!(split(C.name,"."),[4,8]),"."),'.')
        name = join([yr,j_day,p_name],"_")

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
end

#Load data based on path - little redundant to have two arguemnts currently
function load_data(OUTDIR::String, filelist::Array{String,1})
    ar = Array{SeisData,1}(undef,0)
    files = joinpath.(expanduser(OUTDIR),filelist)
    for i in 1:length(filelist)
        data = read_data("mseed",files[i])
        push!(ar, data)
    end
    return(ar)
end

#Add location from dataframe to array 
function add_loc_df(ar::Array{SeisData,1},df::DataFrame)
    for i in 1:length(ar)
        try
            station_name = split(ar[i].name[1],".")[2]
            for j in 1:size(df)[1]
                if station_name == df[j,2]
                    ar[i].loc[1].lat = df[j,3]
                    ar[i].loc[1].lon = df[j,4]
                    ar[i].loc[1].el = df[j,5]
                    break
                end
            end
        catch
            println("Station $(i) doesn't work!")
        end
    end
end

#Returns indices of source stations at index 1 and non-source stations at index 2
function indexer(ar::Array{String,1}, stations::Array{String,1})
    indices = Array{Int64,1}(undef,0)
    for i in 1:length(stations)
        for j in 1:length(ar)
            if occursin(stations[i],ar[j]) == true
                push!(indices, j)
            end
        end
    end
    not_indices = setdiff(1:length(ar), indices)
    return([indices,not_indices])
end

#Heatmap plot function - could be improved, some arguments are inbedded
function heatmapper(corrs::Array{CorrData,1}, maxbin::Float64, bin_size::Float64, minfreq::Float64, maxfreq::Float64)
    bins = collect(0.:bin_size:maxbin-bin_size)
    bins1 = collect(bin_size:bin_size:maxbin)
    counts = zeros(length(bins))
    # freqmin = 1.
    # freqmax = 2.
    C = corrs[1]
    Cmat = zeros(eltype(C.corr),size(C.corr,1),length(bins))
    for ii = 1:length(corrs)
        C = corrs[ii]
        ind = findall((C.dist .>= bins) .& (C.dist .<= bins1))[1]
        bandpass!(C, minfreq,maxfreq)
        clean_up!(C,minfreq,maxfreq)
        abs_max!(C)
        Cmat[:,ind] .+= C.corr[:]
        counts[ind] += 1
        #println("Reading file $ii")
    end

    maxdist = 600.
    ind = findall((bins .> 0.15) .& (bins .<= maxdist))
    Cmat = Cmat[:,ind]
    counts = counts[ind]
    bins = bins[ind]
    Cmat ./= counts'
    lags = -C.maxlag:1/C.fs:C.maxlag
    maxlag = 500.
    lagind = findall(abs.(lags) .<= maxlag)
    Cmat = Cmat[lagind,:]
    abs_max!(Cmat)
    # plot envolopes at 5 km/s and 1.3 km/s
    heatmap(
        lags[lagind],
        bins,
        Cmat',
        c=:balance,
        xlabel="Lag [s]",
        ylabel="Inter-station Distance [km]",
        legend = :none,
        dpi=500,
    )
    plot!([0.,80.],[0.1,400],line=(1.5,:dash),color=:black,label="")
    plot!([0.,-80.],[0.1,400],line=(1.5,:dash),color=:black,label="")
    plot!([0.,267.],[0.03,400],line=(1.5,:dash),color=:black,label="")
    plot!([0.,-267.],[0.03,400],line=(1.5,:dash),color=:black,label="")
    annotate!(2.5,17, Plots.text("5 km/s", 14, :dark, rotation = 78 ))
    annotate!(200,300, Plots.text("1.5 km/s", 14, :dark, rotation = 48 ))
    ylims!((minimum(bins),maximum(bins)))
    xlims!((-maxlag,maxlag))
    png(expanduser("plots/moveouts/heatplot_test.png"))
end

# Send bucket and date range for correlations to all cores
@everywhere begin
    using SCEDC, AWSCore, Dates 
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    startdate = "2017-01-01" # Select Start Date
    enddate = "2017-02-25" # Select End Date
    days = Date(startdate):Day(1):Date(enddate)
    network = "CI"
    channel = "BH?"
    OUTDIR = "~/data"
end
t1 = Dates.now()
#Loop through days for scedc-pds data download
for i in 1:length(days)
    filelist = s3query(aws, days[i], enddate = days[i], network=network, channel=channel)
    #println(filelist)
    @eval @everywhere filelist=$filelist
    #path = join([Dates.year(days[day]),lpad(Dates.dayofyear(days[day ]),3,"0")],"_") # Yeilds "YEAR_JDY"
    ec2download(aws,bucket, filelist, OUTDIR)
end

#Loop through days for seisbasin data download 
for i in 5:length(days)
    yr = Dates.year(days[i])
    path = join([Dates.year(days[i]),lpad(Dates.dayofyear(days[i]),3,"0")],"_") # Yeilds "YEAR_JDY"
    #S3 query call for data file information
    dict = collect(s3_list_objects(aws, "seisbasin", "continuous_waveforms/$(yr)/$(path)/", max_items=1000))
    filelist = Array{String,1}(undef,length(dict))
    #Index to filepath given by the "Key" element of the dictionary
    for i in 1:length(dict)
        filelist[i] = dict[i]["Key"]
    end
    @eval @everywhere filelist=$filelist
    ec2download(aws, bucket2, filelist, OUTDIR)
end

# Read in station locations and list source stations
df = CSV.File("files/all_locations_socal.csv") |> DataFrame! 
#sources = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV"]
sources = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV"
        ,"VINE", "ROPE", "ARNO", "LUCI", "ROUF", "KUZD", "ALLI"]

t_start = Dates.now()
#Preprocess and correlate each day
for i in 1:length(days)
    GC.gc() # Clean up Memory 
    #Get list of raw data file paths
    path = join([Dates.year(days[i]),lpad(Dates.dayofyear(days[i]),3,"0")],"_") # Yeilds "YEAR_JDY"
    fpaths = readdir("/home/ubuntu/data/continuous_waveforms/$(Dates.year(days[i]))/$(path)")
    files = joinpath.("/home/ubuntu/data/continuous_waveforms/$(Dates.year(days[i]))/$(path)",fpaths)


    #Load data into array 
    ar = Array{SeisData,1}(undef,0) # parallelize the reading of data into array
    files_clean = deepcopy(files); # Update files_clean to have only the clean data 
    for i in 1:length(files)
        data = read_data("mseed",files[i])
        # Check for gaps and sufficient data
        gaps = size(data.t[1])[1] # N-2 
        pts = size(data.x[1])[1]
        if gaps < 10 && pts > 2/3*fs*24*3600 #Ensure at least 2/3 of the data is present
            # If few gaps, data is clean
            push!(ar, data)
        else
            # Remove index from filelist
            deleteat!(files_clean, i)
        end
    end
    print
    add_loc_df(ar, df) # add location of array data from dataframe

    # Create DataFrame and coerce to CSV for seisbasin indexing
    sta_index = DataFrame(ms_filename = String[],network = String[],station = String[], location = String[], channel = String[],lat = String[], lon=String[], el = String[], sample_rate=String[])
    for i in 1:length(ar)
        f_name = split(files_clean[i],'/')[end]
        sta_name = split(ar[i].name[1],'.')
        push!(sta_index, Dict(:ms_filename => "$f_name", :network => "$(sta_name[1])", 
                              :station => "$(sta_name[2])", :location=> "$(sta_name[3])",
                              :channel => "$(sta_name[4])", :lat => "$(ar[i].loc[1].lat)", 
                              :lon => "$(ar[i].loc[1].lon)", :el => "$(ar[i].loc[1].el)", :sample_rate => "$fs"))
    end
    #Write DataFrame to csv - create path if necessary
    CORROUT = expanduser("files/index/$(Dates.year(days[i]))/")
    if isdir(CORROUT) == false
        mkpath(CORROUT)
    end
    CSV.write("files/index/$(Dates.year(days[i]))/$(path)_correlation_index.csv", sta_index)
    s3_put(aws, "seisbasin", "corr_index/$(Dates.year(days[i]))/$(path)_correlation_index.csv",read("files/index/$(Dates.year(days[i]))/$(path)_correlation_index.csv"))
#end # to get just the csv
    #Preprocess Data
    T1 = @elapsed ffts = pmap(x -> preprocess2(x, fs, freqmin, freqmax, cc_step, cc_len, half_win, water_level), ar)
    println("Data for $(days[i]) preprocessed in $(T1) seconds. Correlating now....")

    indices = indexer(files_clean, sources) # returns indices of sources
    pairs = Array{Array{Int64,1}}(undef, 0) #array of station pairs to correlate
    for i in 1:length(indices[1])
        for j in 1:length(indices[2])
            push!(pairs, [indices[1][i],indices[2][j]])
        end
    end

    #Correlate Station pairs
    T2 = @elapsed pmap(x -> correlatedir(x, maxlag), map(y -> ffts[y], pairs))
    println("Data for $(days[i]) correlated in $(T2) seconds!")


    # transfer that day of data 
    corr_paths = glob("*/*/$(path)*.jld2","CORR")
    Transfer = @elapsed pmap(x ->s3_put(aws, "seisbasin", "corr_data/$(join(deleteat!(split(x, "/"),[1]),"/"))", read(x)), corr_paths)
    println("$(length(corr_paths)) correlations transfered to $(bucket2) in $Transfer seconds!")
end
t2 = Dates.now()
t_diff = t2-t1
println("Correlations completed in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(t_diff)))).")
