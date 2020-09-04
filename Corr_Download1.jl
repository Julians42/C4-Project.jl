T = @elapsed using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays, AWSS3

#coeffs
cc_step, cc_len = 3600, 3600
maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9
half_win, water_level = 30, 0.01

function indexer(ar::Array{String,1}, stations::Array{String,1})
    """
        Returns Indices of stations and recievers (others) in ar using occursin
    """
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

maxbin, bin_size = 300., 10. 
minfreq, maxfreq = 0.1, 0.2
#Heatmap plot function - could be improved, some arguments are inbedded
function heatmapper(corrs::Array{CorrData,1}, maxbin::Float64, bin_size::Float64, minfreq::Float64, maxfreq::Float64)
    """
        Function takes an array of correlation data, filters the data between minfreq and maxfreq

        maxbin is the cutoff for the maximum distance between correlations (in km), larger distances are discarded,
            the bin_size intuitively gives the range of correlation distances which are stacked together.

        A good starting place to check visibility of large-scale correlations is: 
            maxbin, bin_size = 600., 10. 
            minfreq, maxfreq = 0.1, 0.2
    """
    bins = collect(0.:bin_size:maxbin-bin_size)
    bins1 = collect(bin_size:bin_size:maxbin)
    counts = zeros(length(bins))
    # freqmin = 1.
    # freqmax = 2.
    C = corrs[1]
    Cmat = zeros(eltype(C.corr),size(C.corr,1),length(bins))
    for ii = 1:length(corrs)
        C = corrs[ii]
        ind = 0
        try
            ind = findall((C.dist .>= bins) .& (C.dist .<= bins1))[1]
        catch
            continue
        end
        bandpass!(C, minfreq,maxfreq)
        clean_up!(C,minfreq,maxfreq)
        abs_max!(C)
        Cmat[:,ind] .+= C.corr[:]
        counts[ind] += 1
        #println("Reading file $ii")
    end

    maxdist = 300.
    ind = findall((bins .> 0.15) .& (bins .<= maxdist))
    Cmat = Cmat[:,ind]
    counts = counts[ind]
    bins = bins[ind]
    Cmat ./= counts'
    lags = -C.maxlag:1/C.fs:C.maxlag
    maxlag = 300.
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
    CORROUT = expanduser("plots/moveouts/")
    if isdir(CORROUT) == false
        mkpath(CORROUT)
    end
    png(expanduser("plots/moveouts/heatplot_test.png"))
end
sources = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV"
        ,"VINE", "ROPE", "ARNO", "LUCI", "ROUF", "KUZD", "ALLI"]

# Make csv of keys 
@everywhere begin
    using SCEDC, AWSCore, Dates 
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    startdate = "2017-01-26" # Select Start Date
    enddate = "2017-02-02" # Select End Date
    days = Date(startdate):Day(1):Date(enddate)
    network = "CI"
    channel = "BH?"
    OUTDIR = "~/data"
end
dy = 1 # For debugging purposes working with 1 day
year = Dates.year(days[dy])
path = join([Dates.year(days[dy]),lpad(Dates.dayofyear(days[dy]),3,"0")],"_") # Yeilds "YEAR_JDY"
#Get CSV
corr_info = CSV.read("files/index/$year/$(path)_correlation_index.csv")

#Create station names from CSV
names = Array{String,1}(undef,0)
for i in 1:size(corr_info)[1]
    location = corr_info.location[i]
    if ismissing(location) == true #Deal with empty location information
        location = ""
    end
    name = join([corr_info.network[i],corr_info.station[i], location,corr_info.channel[i]],".")
    push!(names, name)
end

#Index for correlation pairs
source_index = indexer(names,sources)
len_src, len_rec = length(source_index[1]), length(source_index[2])

# Create DataFrame to store the names of the correlations with components 
corr_pairs = DataFrame(corr_name = String[], comp = String[]) #Pairs is length of sources * receivers
for i in 1:len_src
    for j in 1:len_rec
        pth = strip(join(vcat(split(names[source_index[1][i]],".")[1:3],split(names[source_index[2][j]],".")[1:3]),"."),'.')
        comp = join([names[source_index[1][i]][end],names[source_index[2][j]][end]],"")
        push!(corr_pairs, Dict(:corr_name => "$pth", :comp => "$comp"))
    end
end

# Make file paths by joining pth and comp from corr_pairs and file routing 
file_paths = Array{String, 1}(undef, length(corr_pairs.corr_name)); 
for i in 1:length(corr_pairs.corr_name)
    pth = "/home/ubuntu/CORR/$(corr_pairs.corr_name[i])/$(corr_pairs.comp[i])/$(path)_$(corr_pairs.corr_name[i]).jld2"
    file_paths[i] = pth
end

#Load all correlations into array
corrs = Array{CorrData,1}(undef,length(file_paths));
for i in 1:length(file_paths)
    try
        data = load_corr(file_paths[i],"$(corr_pairs.comp[i])")
        corrs[i] = data
    catch
        println("$(file_paths[i]) load unsuccessful!")
    end
end

#Make a moveout plot of all the corrs
heatmapper(corrs, 300., 10., 0.1, 0.2)

#Push moveout plot to s3 bucket (might want to update 3rd argument for file pathing)
s3_put(aws, "seisbasin","plots/$year/heatplot_$path.png", "plots/moveouts/heatplot_$path.png")