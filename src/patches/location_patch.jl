# Patch Script for adding locations to processed h5 files
# Can only add locations once, script needs to be debugged to handle sparse data
# packages
T = @elapsed using SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Statistics, PyCall, Glob, StructArrays, AWSS3, HDF5
addprocs()
@everywhere begin
    using SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWSCore, StructArrays, AWSS3, Statistics, JLD2, Glob, HDF5
    function LLE_geo(network, station, df)
        """ Find station matching location and return geoloc object"""
        try
            station_indices, network_indices = findall(x -> x==station, df.station), findall(y -> y==network, df.network)
            row = df[collect(intersect(Set(station_indices), Set(network_indices)))[1],:]
            lat, lon, el = row.latitude[1], row.longitude[1], row.elevation[1]
            geo = GeoLoc(lat = float(lat), lon = float(lon), el = float(el))
            return geo
        catch 
            return nothing
        end
    end

    function hdf_addloc(file::String, df::DataFrame)
        """ add distance, azimuth, and back azimuth to processed h5 file """
        net_src, sta_src = split(split(file,"/")[end],".")[1], split(split(file,"/")[end],".")[2]
        geo_src = LLE_geo(net_src, sta_src, locations)
        h5open(file,"cw") do fid
            # get recievers
            reciever_stations = keys(read(fid))
            A = read(fid) # read file to access keys
            for rec in reciever_stations
                # Avoid adding to metafile and skip if all fields already filled 
                if rec == "meta"
                    continue
                end
                if haskey(A, "$rec/dist") && haskey(A, "$rec/azi") && haskey(A, "$rec/baz")
                    continue
                end
                # get reciever location
                println(rec) # optional print
                net_rec, sta_rec = split(rec,".")[1], split(rec,".")[2]
                geo_rec = LLE_geo(net_rec, sta_rec, df)
                dist, azi, baz = get_dist(geo_src, geo_rec), get_azi(geo_src, geo_rec), get_baz(geo_src, geo_rec)
                # write datafields
                if !haskey(A, "$rec/dist")
                    write(fid, "$rec/dist", dist)
                end
                if !haskey(A, "$rec/azi")
                    write(fid, "$rec/azi", azi)
                end
                if !haskey(A, "$rec/baz")
                    write(fid, "$rec/baz", baz)
                end
            end
        end
    end
end


# select files to redownload and process
dict = collect(s3_list_objects(aws, "seisbasin", "source_processed/2019/"))
filelist_basin = [fpath["Key"] for fpath in dict]

ec2download(aws, bucket2, filelist_basin, "~/")

# dataframe with all socal locations
locations = DataFrame(CSV.read("files/full_socal.csv"))

#Add location from dataframe to array 
T_loc = @elapsed pmap(x->hdf_addloc(x, locations), filelist_basin)

# reupload to seisbasin
Transfer = @elapsed pmap(x -> s3_put(aws, "seisbasin",x, read(x)), filelist_basin)



########################## EOS #############################
# print recievers as test
h5open(filelist_basin[1], "r") do fid
    recs = keys(read(fid))
    println(recs)
end

A = read(filelist_basin[1])
A[collect(keys(collect(keys(A))[1]))[1]]
# example
A = read("sourcename.h5")
A["$(collect(keys(A)))"] # returns all sources including "meta"
A["$(collect(keys(A))[1])"]["NN"]["pws"] # to access data for 1st station
A["$(collect(keys(A))[1])"]["dist"] # to access distance, azimuth, and baz