using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames,  SCEDC, AWSCore, Distributed, JLD2, Statistics

cc_step, cc_len = 3600, 3600
maxlag, fs = 1500., 40. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9

#Plot coefficients
frequency_plots = [[0.1,0.2],[0.2,0.5],[0.5,1.],[1.,2.],[2.,5.],[5.,10.],[10.,20.]]
scale, lw = 0.5, 0.5 #compress waveforms and thin line width


addprocs()

@everywhere begin
    using SCEDC, AWSCore, Dates 
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "basin"
    startdate = Date("2018-07-01")
    enddate = Date("2018-07-01")
    network = "CI"
    channel = "BH?"
    OUTDIR = "~/data"
end
filelist = s3query(aws, startdate, enddate=enddate, network=network,channel=channel)
@eval @everywhere filelist=$filelist
# do transfer
ar = ec2stream(aws,bucket,filelist);

# Gets station location information and add
df = CSV.File("permanent_stations_socal.csv") |> DataFrame! 

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
add_loc_df(ar,df)

#pick out sources by name (or any string)
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
# 2 element array of arrays storing source and reciever indices
indices = indexer(filelist, ["SVD","TA2"])

pairs = Array{Array{Int64,1}}(undef, 0)
for i in 1:length(indices[1])
    for j in 1:length(indices[2])
        push!(pairs, [indices[1][i],indices[2][j]])
    end
end


@everywhere begin
    using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames,SCEDC, AWSCore, Distributed
    function preprocess2(ar::SeisData, fs::Float64, freqmin::Float64, freqmax::Float64, cc_step::Int64, cc_len::Int64)
        process_raw!(ar, fs)
        R = RawData(ar,cc_len,cc_step)
        detrend!(R)
        taper!(R)
        bandpass!(R,freqmin,freqmax,zerophase=true)
        FFT = compute_fft(R) 
        R = nothing
        whiten!(FFT,freqmin,freqmax)
        return(FFT)
    end
end
ffts = pmap(x -> preprocess2(x, fs, freqmin, freqmax, cc_step, cc_len), ar)

#Delete ar 
ans = nothing
#ar = nothing
GC.gc()




# Correlate function with stacking
@everywhere begin 
    using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames,SCEDC, AWSCore, Distributed, Statistics
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
    function correlate5(ffts::Array{FFTData,1}, freqmin::Float64, freqmax::Float64, maxlag::Float64)
        C = correlate(ffts[1],ffts[2],maxlag) 
        #clean_up!(C,freqmin,freqmax)
        cc_medianmute!(C, 10.)  #filter extraneous indices
        stack!(C)
        save_corr(C, "~/CORR")
    end
end
T5 = @elapsed pmap(x -> correlate5(x, freqmin, freqmax, maxlag), map(y -> ffts[y], pairs))


corr_name = joinpath.("CORR/",readdir("CORR"))
corrs = Array{CorrData,1}(undef,0)
@everywhere begin
    using JLD2
    function testfunc(x)
        t = jldopen("$(x)")
        push!(corrs, t)
        close(t)
    end 
end
pmap(x-> testfunc(x), corr_name)

#load corrs into array
function testfunc(names::Array{String,1})
    corrs = Array{CorrData,1}(undef,0)
    for i in 1:length(names)
        comp  = "$(split(names[i],".")[4][end])$(split(names[i],".")[end-1][end])"
        t = load_corr(names[i], comp)
        push!(corrs, t)
    end
    return(corrs)
end 
corrs = testfunc(corr_name)
#map(x-> testfunc(x), corr_name)

#write into one file
fo = jldopen("july2018.jld2", "w") # name of the station pairs, most computationally expensive
for (index, value) in enumerate(corrs)
    # continue again if xcorr is empty
	isempty(value.corr) && continue
	pair = join(deleteat!(split(corrs[index].name,"."),[3,4,7,8]),".")
    comp = value.comp
    starttime = string(u2d(corrs[index].t[1]))
	groupname = joinpath(pair, comp, starttime)
	!haskey(fo, pair) && JLD2.Group(fo, pair) # if it doesn't have pair, make pair layer
	!haskey(fo[pair], comp) && JLD2.Group(fo[pair], comp) # component
	!haskey(fo, groupname) && (fo[groupname] = value) # save corr into layer 
end
close(fo)

#sudo apt install awscli
#Upload file to AWS bucket
s3_put(aws, "seisbasin", "2018/$(name)", read("july2018.jld2"))