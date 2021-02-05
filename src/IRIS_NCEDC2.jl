# Working on speedup for IRIS download
T = @elapsed using SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Glob, AWSS3

rootdir = ""
#rootdir = "/Users/julianschmitt/Downloads/Seismo/IRIS/" # set to "" if on EC2/ in container, but can run locally 
# Functions 
aws = aws_config(region="us-west-2")
# aws = aws_config(creds = AWSCredentials("KEY", "SECRET_KEY"), region="us-west-2") # if running locally - careful with AWS creds!
bucket = "scedc-pds"
bucket2 = "seisbasin"
network = "CI"
channel1 = "BH?"
channel2 = "HH?"
OUTDIR = "~/data"
function download_data(station, source::String, dday::Date)
    try
        # download data 
        data = get_data("FDSN", station, src=source, s=string(dday), t = string(dday+Day(1)))
        println("Downloaded data for $station")
        # save data
        for datum in [data[1], data[2], data[3]]
            try
                println(datum.id)
                if (source == "NCEDC") && (datum !=nothing)
                    println(length(datum.x))
                    if length(datum.x)>500000 # check second since NCEDC and IRIS send different formats
                        #println(joinpath(rootdir, "ncedc_waveforms/$yr/$path/$(datum.id)"))
                        wseis(joinpath(rootdir, "ncedc_waveforms/$yr/$path/$(datum.id)"), datum)
                    end
                else # data is from IRIS
                    if (datum !=nothing) && (length(SeisData(datum).x[1])>500000)
                        println("SeisData from IRIS, waveform length: $(length(SeisData(datum).x[1]))")
                        wseis(joinpath(rootdir, "iris_waveforms/$yr/$path/$(datum.id)"), SeisData(datum))
                    end
                end
            catch e
                println(e)
            end
        end
    catch e
        println(e)
    end
end

function query_name(s::String, full::Bool)
    name = convert(String, split(s,"/")[end])
    net, cha, sta, comp = name[1:2], strip(convert(String, name[3:6]),'_'), 
                          strip(convert(String, name[11:12]),'_'), name[10]
    if full
        return join([net, cha, sta, comp],".")
    else  
        return join([net, cha, comp],".") #comp - add back in for type 
    end
end
function get_dict_name(file::String)
    station = convert(String, split(split(file,"/")[end],"_")[1])
    component = split(file,"/")[end][10:12]
    return string(station, "_", component)
end
function query_FDSN(file_str::String, src::String, dday::Date)
    """ Query NCEDC and IRIS for metadata """
    println("Beginning metadata download for $src ($file_str).")
    if src in ["SCEDC","NCEDC"]
        return FDSNsta(file_str, src=src, s = string(dday), t = string(dday+Day(1))).id
    else # then data is IRIS so we filter to region
        return FDSNsta(file_str,src=src, s = string(dday), t = string(dday+Day(1)), reg=[31.,40.,-123.,-116.]).id
    end
end

arg = ENV["AWS_BATCH_JOB_ARRAY_INDEX"]
startdate = Date(2000)+Month(arg)
enddate = startdate+Month(1)-Day(1)
println(startdate, enddate)
days = startdate:Day(1):enddate

data_sources = [["*.*.*.HH*","NCEDC"],["*.*.*.BH*","NCEDC"],["*.*.*.HH*","IRIS"],["*.*.*.BH*","IRIS"]]
#list_queries, NCEDC_list, IRIS_list = [], [], []

for i in 1:length(days)
    yr = Dates.year(days[i])
    path = join([Dates.year(days[i]),lpad(Dates.dayofyear(days[i]),3,"0")],"_") # Yeilds "YEAR_JDY"
    @eval @everywhere yr, path, rootdir = $yr, $path, $rootdir
    println("Processing day $path...")
    ##################### Query SCEDC S3 for easy stations ########################
    ar_filelist = map(x -> s3query(aws, days[i], enddate = days[i], network=network, channel=x),[channel1, channel2])
    filelist_scedc_BH = ar_filelist[1]
    filelist_scedc_HH = ar_filelist[2]

    # create dictionary and overwrite HH keys with available BH data
    BH_keys = [get_dict_name(file) for file in filelist_scedc_BH]
    HH_keys = [get_dict_name(file) for file in filelist_scedc_HH]

    # Convert to dictionary 
    HH_dict = Dict([(name, file) for (name, file) in zip(HH_keys, filelist_scedc_HH)]) 
    BH_dict = Dict([(name, file) for (name, file) in zip(BH_keys, filelist_scedc_BH)]) 
    filelist_dict = merge(BH_dict, HH_dict) # BH dict overwrite HH_dict. This is essentually the union
    filelist_scedc = collect(values(filelist_dict)) # return values as array for download

    scedc_list = map(x -> query_name(x, true), filelist_scedc)
    scedc_brief = map(x -> query_name(x, false), filelist_scedc)
    scedc_sta = [split(elt, ".")[2] for elt in scedc_brief]

    save_sources = DataFrame(Bucket = vec(fill("SCEDC", (length(scedc_brief),1))), Source = scedc_list, Station = scedc_sta, ID = scedc_brief, Path = filelist_scedc)
    # write scedc name files immediately to bucket - donwload during correlate routine
    if !isdir(joinpath(rootdir, "scedc_path/$yr/"))
        mkpath(joinpath(rootdir, "scedc_path/$yr/"))
    end
    CSV.write(joinpath(rootdir, "scedc_path/$yr/$path.csv"), save_sources)
    s3_put(aws, "seisbasin", "scedc_path/$yr/$path.csv", read(joinpath(rootdir, "scedc_path/$yr/$path.csv")))
    println("There are $(length(save_sources.ID)) SCEDC stations for $path.")
    ####################### Query NCEDC and IRIS for extra missing stations ######################
    T_list = @elapsed list_queries = map(x -> query_FDSN(x[1],x[2], days[i]), data_sources) # pmap doesn't work bc server ....
    map(list -> filter!(x -> x[end] in ['E','N','Z'], list), list_queries) # Restrict to ENZ channels
    data_sources_ch = [[join(split(elt,".")[1:2],".") for elt in array] for array in list_queries] 
    data_sources_ch = [[join([split(elt,".")[1], split(elt,".")[2], elt[end]],".") for elt in array] for array in list_queries] 

    for (ind, source) in enumerate(data_sources_ch) # add new stations to dataframe
        for (j, station) in enumerate(source)
            if station ∉ save_sources.ID # check if already contained in list (lower priority stations won't be added)
                push!(save_sources, [data_sources[ind][2], list_queries[ind][j], split(station,".")[2], station, list_queries[ind][j]])
            end
        end
    end

    ####################### Download NCEDC and IRIS Data ######################
    # Write data and then upload to seisbasin
    if !isdir(joinpath(rootdir, "ncedc_waveforms/$yr/$path/"))
        mkpath(joinpath(rootdir, "ncedc_waveforms/$yr/$path/"))
    end
    if !isdir(joinpath(rootdir, "iris_waveforms/$yr/$path/"))
        mkpath(joinpath(rootdir, "iris_waveforms/$yr/$path/"))
    end

    NCEDC_list = filter(x -> x.Bucket =="NCEDC", save_sources).Path
    IRIS_list = filter(x -> x.Bucket =="IRIS", save_sources).Path
    println("Now trying download for $(length(NCEDC_list)) NCEDC stations and $(length(IRIS_list)) IRIS stations.")

    IRIS_stations = unique([string(x[1:end-1],"*") for x in IRIS_list])
    NCEDC_stations = unique([string(x[1:end-1],"*") for x in NCEDC_list])
    map(x -> download_data(x, "IRIS", days[i]),IRIS_stations)
    map(x -> download_data(x, "NCEDC", days[i]),NCEDC_stations)
    ####################### Filter and Upload data to seisbasin ######################
    # Can use glob on ec2 - doesn't like "/" on local
    iris_filepaths = joinpath.("iris_waveforms/$yr/$path/", readdir(joinpath(rootdir, "iris_waveforms/$yr/$path/"))) # duplicates contain "[NAME]" so we filter these
    ncedc_filepaths = joinpath.("ncedc_waveforms/$yr/$path/", readdir(joinpath(rootdir, "ncedc_waveforms/$yr/$path/")))
    filter!(x -> !occursin(".DS_Store", x), iris_filepaths)
    filter!(x -> !occursin(".DS_Store", x), ncedc_filepaths)

    map(x -> s3_put(aws, "seisbasin", x, read(joinpath(rootdir, x))), iris_filepaths);
    map(x -> s3_put(aws, "seisbasin", x, read(joinpath(rootdir, x))), ncedc_filepaths);

    println("Transfered $(length(iris_filepaths)+length(ncedc_filepaths)) raw seisdata objects from NCEDC and IRIS to seisbasin.")

    # Cleanup - doesn't delete directories because AWS batch gets mad
    rm.(joinpath.(rootdir, "iris_waveforms/$yr/$path", readdir(joinpath(rootdir, "iris_waveforms/$yr/$path"))), recursive=true)
    rm.(joinpath.(rootdir, "ncedc_waveforms/$yr/$path", readdir(joinpath(rootdir, "ncedc_waveforms/$yr/$path"))), recursive=true)
end
println("Finished", startdate, enddate)
