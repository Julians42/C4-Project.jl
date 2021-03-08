export download_data, get_scedc_files, query_name, get_dict_name, query_FDSN, evaluate_done, get_seisdata



function download_data(station, source::String, dday::Date, yr::Int64, path::String, rootdir::String="")
    try
        # download data 
        data = get_data("FDSN", station, src=source, s=string(dday), t = string(dday+Day(1)))
        println("Downloaded data for $station on $dday")
        # save data
        for ind in 1:3
            try
                datum = data[ind]
                if (source == "NCEDC") && (datum !=nothing)
                    if length(datum.x)>500000 # check second since NCEDC and IRIS send different formats
                        #println(joinpath(rootdir, "ncedc_waveforms/$yr/$path/$(datum.id)"))
                        wseis(joinpath(rootdir, "ncedc_waveforms/$yr/$path/$(datum.id)"), datum)
                    end
                else # data is from IRIS
                    if (datum !=nothing) && (length(SeisData(datum).x[1])>500000)
                        println("SeisData from IRIS, waveform length: $(length(SeisData(datum).x[1]))")
                        fpath = joinpath(rootdir, "iris_waveforms/$yr/$path/$(datum.id)")
                        wseis(fpath, SeisData(datum)) # write seismic data 
                        # upload to seisbasin and clean up
                        s3_put(aws, "seisbasin", fpath, read(fpath))
                        rm(fpath)
                    end
                end
            catch e
                println("Data Integrity Check Failed!")
            end
        end
    catch e
        println("Could not download data for $station.")
    end
end
function get_scedc_files(dd::Date)
    ar_filelist = pmap(x -> s3query(aws, dd, enddate = dd, network=network, channel=x),["BH?", "HH?"])
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
    return filelist_scedc
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
    """ Query NCEDC and IRIS for metadata. Try twice """
    println("Beginning metadata download for $src ($file_str).")
    try
        if src in ["SCEDC","NCEDC"]
            return FDSNsta(file_str, src=src, s = string(dday), t = string(dday+Day(1))).id
        else # then data is IRIS so we filter to region
            return FDSNsta(file_str,src=src, s = string(dday), t = string(dday+Day(1)), reg=[31.,40.,-123.,-116.]).id
        end
    catch e
        println(e)
        println("Could not download metadata for $file_str from $src on $dday.")
        try
            if src in ["SCEDC","NCEDC"]
                return FDSNsta(file_str, src=src, s = string(dday), t = string(dday+Day(1))).id
            else # then data is IRIS so we filter to region
                return FDSNsta(file_str,src=src, s = string(dday), t = string(dday+Day(1)), reg=[31.,40.,-123.,-116.]).id
            end
        catch
            println("Could not download metadata for $file_str from $src on $dday.")
            println("Error on both attempts.")
        end
    end
end
function evaluate_done(yr, path)
    try
        #iris_query = [elt["Key"] for elt in collect(s3_list_objects(aws, "seisbasin", "iris_waveforms/$yr/$path/", max_items=1000))]
        #ncedc_query = [elt["Key"] for elt in collect(s3_list_objects(aws, "seisbasin", "ncedc_waveforms/$yr/$path/", max_items=1000))]
        iris_path, ncedc_path = S3Path("s3://seisbasin/iris_waveforms/$yr/$path", aws), S3Path("s3://seisbasin/ncedc_waveforms/$yr/$path", aws)
        iris_query, ncedc_query = convert.(String, readdir(iris_path)), convert.(String, readdir(ncedc_path))
        println("$(length(iris_query)-1) iris waveforms and $(length(ncedc_query)-1) ncedc waveforms already downloaded for $path.")
        yr_num = convert(Int64, yr)
        if ((length(iris_query) < 20) && (yr_num <=2004)) || ((length(iris_query)<100) && (yr_num >=2005)) || ((length(ncedc_query) < 5) && yr_num <=  2003) || ((length(ncedc_query) < 80) && (yr_num >= 2004))
            println("Continuing to download phase for $path...")
            return false # we need to attempt to reprocess this day! 
        else
            println("Data is already in seisbasin for $path. Continuing to next day...")
            return true # day is already processed
        end
    catch e
        println(e)
        println("Error in evaluating day. Continuing to download phase")
    end
end

function get_seisdata(date::Date, data_sources::Array{Array{String,1},1}=[["*.*.*.HH*","NCEDC"],["*.*.*.BH*","NCEDC"],
            ["*.*.*.HH*","IRIS"],["*.*.*.BH*","IRIS"]], rootdir::String="")
    """ Downloads available seismic data from NCEDC and IRIS for given date. Uploads data to seisbasin """
    try 
        yr = Dates.year(date)
        path = join([Dates.year(date),lpad(Dates.dayofyear(date),3,"0")],"_") # Yeilds "YEAR_JDY"
        @eval @everywhere yr, path, rootdir = $yr, $path, $rootdir
        if evaluate_done(yr, path) == true
            return 0
        else
            ##################### Query SCEDC S3 for easy stations ########################
            filelist_scedc = get_scedc_files(date)

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
            T_list = @elapsed list_queries = map(x -> query_FDSN(x[1],x[2], date), data_sources) # pmap doesn't work bc server ....
            for list in list_queries
                try
                    filter!(x -> x[end] in ['E','N','Z'], list)
                catch e
                    println(e)
                    println("Likely no metadata from that source")
                end
            end
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
            map(x -> download_data(x, "IRIS", date, yr, path), IRIS_stations)
            map(x -> download_data(x, "NCEDC", date, yr, path), NCEDC_stations)
            ####################### Filter and Upload data to seisbasin ######################
            # Can use glob on ec2 - doesn't like "/" on local
            # iris_filepaths = joinpath.("iris_waveforms/$yr/$path/", readdir(joinpath(rootdir, "iris_waveforms/$yr/$path/"))) # duplicates contain "[NAME]" so we filter these
            # ncedc_filepaths = joinpath.("ncedc_waveforms/$yr/$path/", readdir(joinpath(rootdir, "ncedc_waveforms/$yr/$path/")))
            # filter!(x -> !occursin(".DS_Store", x), iris_filepaths)
            # filter!(x -> !occursin(".DS_Store", x), ncedc_filepaths)

            # map(x -> s3_put(aws, "seisbasin", x, read(joinpath(rootdir, x))), iris_filepaths);
            # map(x -> s3_put(aws, "seisbasin", x, read(joinpath(rootdir, x))), ncedc_filepaths);

            # println("Transfered $(length(iris_filepaths)+length(ncedc_filepaths)) raw seisdata objects from NCEDC and IRIS to seisbasin.")

            # # Cleanup - doesn't delete directories because AWS batch gets mad
            # rm.(joinpath.(rootdir, "iris_waveforms/$yr/$path", readdir(joinpath(rootdir, "iris_waveforms/$yr/$path"))), recursive=true)
            # rm.(joinpath.(rootdir, "ncedc_waveforms/$yr/$path", readdir(joinpath(rootdir, "ncedc_waveforms/$yr/$path"))), recursive=true)
            return 0
        end
    catch e
        println(e)
        println("Unable to process data for $path. Continuing ...")
        return 1
    end
end