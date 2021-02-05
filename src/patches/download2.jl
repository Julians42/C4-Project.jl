using SeisIO, SeisNoise, Distributed, Dates, AWSS3, AWSCore, JLD2, Statistics
addprocs()
@everywhere using SeisIO, SeisNoise, Dates, AWSS3, AWSCore

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

aws = aws_config(region="us-west-2")

# select days for download
startdate = "2019-02-01"
enddate = "2019-02-02"
days = Date(startdate):Day(1):Date(enddate)

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