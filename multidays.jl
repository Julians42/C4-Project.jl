# Install julia, add AIM role with s3 access, install packages.
using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays

@everywhere using SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWSCore, StructArrays
@everywhere begin 
    function preprocess2(ar::SeisData, fs::Float64, freqmin::Float64, freqmax::Float64, cc_step::Int64, cc_len::Int64, half_win::Int64, water_level::Float64)
        process_raw!(ar, fs) # downsample to 20 Hz 
        R = RawData(ar,cc_len,cc_step)
        detrend!(R)
        taper!(R)
        bandpass!(R,freqmin,freqmax,zerophase=true)
        FFT = compute_fft(R) 
        coherence!(FFT,half_win, water_level)
        return(FFT)
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
end

#Patch until SeisNoise removes Plots or updates it
using Pkg 
ENV["GR"] = ""
Pkg.build("GR")

cc_step, cc_len = 3600, 3600
maxlag, fs = 1500., 20. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9
half_win, water_level = 30, 0.01

# Functions 

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

#Choose startdate and enddate
startdate = "2018-07-01"
enddate = "2018-07-02"
days = Date(startdate):Day(1):Date(enddate)


addprocs()
@everywhere begin
    using SCEDC, AWSCore, Dates 
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "basin"
    startdate = "2018-07-01"
    enddate = "2018-07-02"
    days = Date(startdate):Day(1):Date(enddate)
    network = "CI"
    channel = "BH?"
    OUTDIR = "~/data"
end

#Loop through days for data download
for day in 1:length(days)
    filelist = s3query(aws, days[day], enddate = days[day], network=network, channel=channel)
    @eval @everywhere filelist=$filelist
    #path = join([Dates.year(days[day]),Dates.dayofyear(days[day ])],"_") # Yeilds "2018_182"
    ec2download(aws,bucket, filelist, OUTDIR)
end

# Read in station locations and list source stations
df = CSV.File("permanent_stations_socal.csv") |> DataFrame! 
sources = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV"]



# Now focus on 1 day at a time to reduce memory
for day in 1:length(days)

    #Get pathing to read in files (scedc-pds specific). Faster than more query calls but more verbose
    path = join([Dates.year(days[day]),Dates.dayofyear(days[day])],"_") # Yeilds "YEAR_JUL"
    names = readdir("/home/ubuntu/data/continuous_waveforms/$(Dates.year(days[day]))/$(path)")
    files = joinpath.("/home/ubuntu/data/continuous_waveforms/$(Dates.year(days[day]))/$(path)",names)
    
    #Load data into array 
    ar = Array{SeisData,1}(undef,0)
    for i in 1:length(files)
        data = read_data("mseed",files[i])
        push!(ar, data)
    end

    #Add location data and determine correlation indices
    add_loc_df(ar, df)
    indices = indexer(files, sources) # returns indices of sources
    pairs = Array{Array{Int64,1}}(undef, 0) #array of station pairs to correlate
    for i in 1:length(indices[1])
        for j in 1:length(indices[2])
            push!(pairs, [indices[1][i],indices[2][j]])
        end
    end
    println("Data for $(days[day]) loaded.")

    # Preprocess
    T1 = @elapsed ffts = pmap(x -> preprocess2(x, fs, freqmin, freqmax, cc_step, cc_len, half_win, water_level), ar)

    println("Data for $(days[day]) preprocessed in $(T1) seconds.")
    ans = nothing #reduce memory pressure
    GC.gc()
    #Correlate
    T2 = @elapsed pmap(x -> correlate5(x, maxlag, "~/CORR/$(path)"), map(y -> ffts[y], pairs))
    println("Data for $(days[day]) correlated in $(T2) seconds.")
end
