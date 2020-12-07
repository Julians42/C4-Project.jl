# Install julia, add AIM role with s3 access, install packages.
using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics

#Get a timestamp for the start
t_start = Dates.now()

cc_step, cc_len = 3600, 3600
maxlag, fs = 1500., 40. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9

#Plot coefficients
frequency_plots = [[0.1,0.2],[0.2,0.5],[0.5,1.],[1.,2.],[2.,5.],[5.,10.],[10.,20.]]
scale, lw = 0.5, 0.5 #compress waveforms and thin line width

#select day range
# startdate = "2018-07-01"
# enddate = "2018-07-02"
# days = Date(startdate):Day(1):Date(enddate)

addprocs()

@everywhere begin
    using SCEDC, AWSCore, Dates 
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "basin"
    #select day range
    startdate = "2018-07-01"
    enddate = "2018-07-02"
    days = Date(startdate):Day(1):Date(enddate)
    startdate = days[2]
    enddate = days[2]
    network = "CI"
    channel = "BH?"
    OUTDIR = "~/data"
end
filelist = s3query(aws, startdate, enddate=enddate, network=network,channel=channel);
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
sources = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV"]
    #,"VINE", "ROPE", "ARNO", "LUCI", "ROUF", "KUZD", "ALLI"] #Node data
indices = indexer(filelist, sources)

pairs = Array{Array{Int64,1}}(undef, 0)
for i in 1:length(indices[1])
    for j in 1:length(indices[2])
        push!(pairs, [indices[1][i],indices[2][j]])
    end
end

#Patch until SeisNoise removes Plots or updates it
using Pkg 
ENV["GR"] = ""
Pkg.build("GR")

@everywhere begin
    using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames,SCEDC, AWSCore
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
    function correlate5(ffts::Array{FFTData,1}, maxlag::Float64)
        C = correlate(ffts[1],ffts[2],maxlag) 
        #clean_up!(C,freqmin,freqmax)
        cc_medianmute!(C, 10.)  #filter extraneous indices
        stack!(C)
        save_corr(C, "~/CORR")
    end
end
T5 = @elapsed pmap(x -> correlate5(x, maxlag), map(y -> ffts[y], pairs))


corr_name = joinpath.("CORR/",readdir("CORR")) #get array of name strings
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


#write into one file
fo = jldopen("2018-183.jld2", "a+") # name of the station pairs, most computationally expensive
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

t_end = Dates.now()
t_diff = t_end - t_start
println(t_diff)


scale = 5.
function plot_corrs2(corrs::Array{CorrData,1})
    # Add filtered correlations to appropriate plots
    T = collect(-corrs[1].maxlag:1/corrs[1].fs:corrs[1].maxlag)
    for j in 1:length(frequency_plots)
        #Select frequency
        fmin = frequency_plots[j][1]
        fmax = frequency_plots[j][2]
        
        #Determine the y scale axis
        distances = []
        for i in 1:length(corrs)
            push!(distances, corrs[i].dist)
        end
        #y_max = minimum(maximum(distances) +0.5, 300.)
        y_max = 300

        # Define plots 
        EE_plot = plot(title = "EE Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)
        EN_plot = plot(title = "EN Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)
        EZ_plot = plot(title = "EZ Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)
        NN_plot = plot(title = "NN Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)
        NZ_plot = plot(title = "NZ Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)
        ZZ_plot = plot(title = "ZZ Moveout Plot - Filtered $(fmin) to $(fmax) Hz", xlims = (-25,25), ylims = (-1,y_max), color = :black)

        for i in 1:length(corrs)
            if corrs[i].comp == "EE"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(EE_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif corrs[i].comp == "EN"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(EN_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif corrs[i].comp == "EZ"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(EZ_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif corrs[i].comp == "NN"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(NN_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif corrs[i].comp == "NZ"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(NZ_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            elseif corrs[i].comp == "ZZ"
                Cstack = bandpass(SeisNoise.stack(corrs[i]), fmin, fmax)
                plot!(ZZ_plot, T, Cstack.corr / maximum(broadcast(abs,Cstack.corr))*scale .+ corrs[i].dist, fmt = :png, color = :black, linewidth = lw, reuse = false, legend = false)
            end
        end
        png(EE_plot,"plots/EE_plot$(fmin)_$(fmax).png")
        png(EN_plot,"plots/EN_plot$(fmin)_$(fmax).png")
        png(EZ_plot,"plots/EZ_plot$(fmin)_$(fmax).png")
        png(NN_plot,"plots/NN_plot$(fmin)_$(fmax).png")
        png(NZ_plot,"plots/NZ_plot$(fmin)_$(fmax).png")
        png(ZZ_plot,"plots/ZZ_plot$(fmin)_$(fmax).png")
    end
end

# In command line run the following two lines (change based on day)
#sudo apt install awscli
#aws s3 cp 2018-182.jld2 s3://seisbasin/corr_data/2018/


####################### Does not work ###########################
#Upload file to AWS bucket
# name = "2018-182.jld2"
# s3_put(aws, "seisbasin", "2018/$(name)", read("2018-183.jld2"))