using SeisIO, SeisNoise, Distributed, Dates, AWSS3, AWSCore, JLD2, Statistics
addprocs()
@everywhere using SeisIO, SeisNoise, Dates, AWSS3, AWSCore

aws = aws_config(region="us-west-2")

# select days for download
startdate = "2019-02-01"
enddate = "2019-02-02"
days = Date(startdate):Day(1):Date(enddate)

# 2019.032.00.00.00.000.FDSNWS.SCEDC.mseed


# make the list 
# download data
# save into mseed with filename for S3
# upload to S3
# delete data
# loop through days
starttime = "2019-02-01"
endtime = "2019-02-02"
days = Date(startdate):Day(1):Date(enddate)

# make  a list day strings
aws = AWS.aws_config(region="us-west-2")
# # find the station list that day
# # get list from SCEDC
# lscedc = FDSNsta("*.*.*.HH*",src="SCEDC",s=starttime,t=endtime)
# println(lscedc.id)
# # get list from NCEDC
# lncedc = FDSNsta("*.*.*.HH*",src="NCEDC",s=starttime,t=endtime)
# println(lncedc.id)
# # get list from
# liris = FDSNsta("*.*.*.HH*",src="IRIS",s=starttime,t=endtime)
# println(liris.id)
#pmap(x-> download_seisdata(x),days)
function download_seisdata(days)
# this function will loop over
d=Dates.Date(days)
# find the station list that day for HH 100Hz data
# get list from SCEDC
lscedc_100 = FDSNsta("*.*.*.HH*",src="SCEDC",s=string(days[1]),t=string(days[1]+Dates.Day(1)))
println(lscedc_100.id)
lscedc_40 = FDSNsta("*.*.*.BH*",src="SCEDC",s=string(d),t=string(d+Dates.Day(1)))
println(lscedc_40.id)
# make list of station unique from the channels.
indices = [findall(x->x==i,lscedc_40.id) for i in lscedc_100.id]
println(lscedc_40.id[indices])
# get list from NCEDC
lncedc_100 = FDSNsta("*.*.*.HH*",src="NCEDC",s=string(d),t=string(d+Dates.Day(1)))
println(lncedc_100.id)
lncedc_40 = FDSNsta("*.*.*.BH*",src="NCEDC",s=string(d),t=string(d+Dates.Day(1)))
println(lncedc_40.id)
# get list from
liris_100 = FDSNsta("*.*.*.HH*",src="IRIS",s=string(d),t=string(d+Dates.Day(1)))
println(liris_100.id)
lncedc_40 = FDSNsta("*.*.*.BH*",src="NCEDC",s=string(d),t=string(d+Dates.Day(1)))
println(lncedc_40.id)
dir0="continuous/year/year_Doy"
# find the station list that day for BH 40Hz data (if no HH)
for sta in lncedc_100.id
    # make file name
    S = get_data("FDSN",sta,src="NCEDC",s=string(d),t=string(d+Dates.Day(1)),w=true)
    # find the default name and change.
    # change filename
    # upload to S3
    s3_put(aws, "seisbasin", filename,read(filename))
    # delete local file
    println(filename)
    return nothing
end

function download_seisdata(sta::String, d::Date)
    S = get_data("FDSN",sta,src="NCEDC",s=string(d),t=string(d+Dates.Day(1)),w=true)
end

S = get_data("FDSN",lscedc_100.id[1],src="NCEDC",s=string(d),t=string(d+Dates.Day(1)),w=true)

########################################################################################

@everywhere begin 
    cc_step, cc_len = 3600, 3600
    maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
    freqmin, freqmax = 0.05, 9.9
    half_win, water_level = 30, 0.01
    aws = aws_config(region="us-west-2")
    bucket = "scedc-pds"
    bucket2 = "seisbasin"
    network = "CI"
    channel1 = "BH?"
    channel2 = "HH?"
    OUTDIR = "~/data"
end

s3query(aws, d, enddate = d+Dates.Day(1), network="CI", channel="BH?")
function query_name(s::String)
    name = convert(String, split(s,"/")[end])
    net, cha, sta, comp = name[1:2], strip(convert(String, name[3:6]),'_'), 
                          strip(convert(String, name[11:12]),'_'), name[8:10]
    return join([net, cha, sta,""],".") #comp - add back in for type 
end
function get_dict_name(file::String)
    station = convert(String, split(split(file,"/")[end],"_")[1])
    component = split(file,"/")[end][10:12]
    return string(station, "_", component)
end
i=1
ar_filelist = pmap(x -> s3query(aws, days[i], enddate = days[i], network=network, channel=x),[channel1, channel2])
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

# Get SCEDC AWS files for comparison
scedc_brief = map(x -> query_name(x), filelist_scedc)


lscedc_100 = FDSNsta("*.*.*.HH*",src="SCEDC",s=string(days[1]),t=string(days[1]+Dates.Day(1)))
println(lscedc_100.id)
lscedc_40 = FDSNsta("*.*.*.BH*",src="SCEDC",s=string(d),t=string(d+Dates.Day(1)))
println(lscedc_40.id)
# make list of station unique from the channels.
indices = [findall(x->x==i,lscedc_40.id) for i in lscedc_100.id]
println(lscedc_40.id[indices])
# get list from NCEDC
lncedc_100 = FDSNsta("*.*.*.HH*",src="NCEDC",s=string(d),t=string(d+Dates.Day(1)))
println(lncedc_100.id)
lncedc_40 = FDSNsta("*.*.*.BH*",src="NCEDC",s=string(d),t=string(d+Dates.Day(1)))
println(lncedc_40.id)
# get list from
liris_100 = FDSNsta("*.*.*.HH*",src="IRIS",s=string(d),t=string(d+Dates.Day(1)))
println(liris_100.id)
liris_40 = FDSNsta("*.*.*.BH*",src="IRIS",s=string(d),t=string(d+Dates.Day(1)))
println(lncedc_40.id)
@everywhere begin
    function query_FDSN(file_str::String, src::String, dday::Date)
        println("Beginning metadata download for $src ($file_str).")
        if src in ["SCEDC","NCEDC"]
            return FDSNsta(file_str, src=src, s = string(dday), t = string(dday+Day(1))).id
        else # then data is IRIS so we filter to region
            return FDSNsta(file_str,src=src, s = string(dday), t = string(dday+Day(1)), reg=[31.,40.,-123.,-116.]).id
        end
    end
end
data_sources = [["*.*.*.HH*","NCEDC"],["*.*.*.BH*","NCEDC"],["*.*.*.HH*","IRIS"],["*.*.*.BH*","IRIS"]] #["*.*.*.HH*","SCEDC"],["*.*.*.BH*","SCEDC"],
T_list = @elapsed list_queries = pmap(x -> query_FDSN(x[1],x[2], days[i]), data_sources)
#scedc_100, scedc_40, ncedc_100, ncedc_40, iris_100, iris_40 = 
data_sources_ch = [[elt[1:end-3] for elt in array] for array in list_queries]    
data_2 = deepcopy(data_sources_ch)

# returns pairs which need to be downloaded
downloads = Array{Array{String,1}}(undef, 0)
for (ind, location) in enumerate(data_sources_ch)
    to_remove = Array{Int64,1}(undef, 0)
    for (rem, elt) in enumerate(location)
        if elt in scedc_brief # if elt in scedc (on aws) remove from download
            push!(to_remove, rem)
        end
        if (ind >1) && elt in Set(vcat(data_sources_ch[1:ind-1])) # check if in earlier array
            push!(to_remove, rem)
        end
    end
    indices = collect(1:length(location))
    filter!(e->e ∉ to_remove, indices) # filter out indicies 
    uniques = list_queries[ind][indices]
    #filter!(x -> x[end] ∉ ["E","N","Z"], uniques) # might have to figure this out later
    push!(downloads, uniques) # add station names to download 
end
println(downloads)



for elt in data_sources
    println(query_FDSN(elt[1],elt[2], days[1]))
end

FDSNsta(file_str, src=src, s = string(dday), t = string(dday+Day(1)), )


liris_100 = FDSNsta("*.*.*.HH*",src="IRIS", s = string(days[1]), t = string(days[1]+Day(1)), reg=[31.,40.,-123.,-116.])