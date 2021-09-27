# script for clipping waveforms 
# we've used this to download clips of waveforms for basin data around earthquakes

# add packages and Distributed
using Pkg, Distributed, Statistics

# packges
addprocs()
@everywhere using AWS, AWSS3, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, Glob, HDF5

@eval @everywhere aws = AWS.AWSConfig(region="us-west-2")

locations = DataFrame(CSV.File("files/CAstations.csv"))
@eval @everywhere locations = $locations

# Helper functions
@everywhere begin 
    function LLE_geo(s::SeisData, df)
        """ Find station matching location and return geoloc object"""
        network, station = split(s.id[1], ".")[1], split(s.id[1],".")[2]
        try
            row = df[findfirst(x -> x.Station==station && x.Network==network, eachrow(df)),:]
            lat, lon, el = row.Latitude[1], row.Longitude[1], row.Elevation[1]
            geo = GeoLoc(lat = float(lat), lon = float(lon), el = float(el))
            s.loc[1] = geo
        catch e
            println(e)
            try # try lowercase csv format (older csvs)
                row = df[findfirst(x -> x.station==station && x.network==network, eachrow(df)),:]
                lat, lon, el = row.latitude[1], row.longitude[1], row.elevation[1]
                geo = GeoLoc(lat = float(lat), lon = float(lon), el = float(el))
                s.loc[1]=geo
            catch
                return nothing
            end
        end
    end
    function tsplit(S::SeisData, s::DateTime, t::DateTime)
        """ Returns clipped waveform based on start and end times """
        LLE_geo(S, locations)
        sync!(S, s=s, t=t)
        detrend!(S)
        taper!(S)
        filtfilt!(S, fl=0.05, fh=10., np=1)
        return S
    end
    function get_event_data(row)
        """ Clip event data based on event information """

        # get event information and file names
        yr = split(row.julian, "_")[1]
        files = glob("continuous_waveforms/$yr/$(row.julian)/*")
        println(files)
        # clip window of 3 hours for teleseismic events
        startsplit, endsplit = row.datetime, row.datetime+Dates.Minute(180)

        # open file and write clipped waveforms 
        h5open("sample_events/$(row.julian).h5", "cw") do file
            if !haskey(read(file), "meta") # event metadata
                write(file, "meta/date", "$(row.date)")
                write(file, "meta/datetime", "$(row.datetime)")
                write(file, "meta/lat", row.lat)
                write(file, "meta/lon", row.lon)
                write(file, "meta/magnitude", row.magnitude)
                write(file, "meta/julian", row.julian)
            end
            for seis in files # loop through data
                try
                    # grab data, split, and add location
                    f = tsplit(read_data("mseed", seis), startsplit, endsplit) # split from t=0 to 10 minutes after
                    network, station = split(f.id[1], ".")[1], split(f.id[1],".")[2]
                    comp = string(f.id[1][end])
                    LLE_geo(f, locations) # add location 
                    dist = get_dist(f.loc[1], GeoLoc(lat=row.lat, lon=row.lon)) # get distance
                    #if !haskey(read(file), "$station/meta") && !haskey(read(file), "$station/meta/dist")
                    try
                        # write receiver metadata
                        write(file, "$station/meta/dist", dist)
                        write(file, "$station/meta/lat", f.loc[1].lat)
                        write(file, "$station/meta/lon", f.loc[1].lon)
                        write(file, "$station/meta/el", f.loc[1].el)
                        write(file, "$station/meta/id", f.id[1])
                    catch
                        println("Failed in metadata")
                    end
                    #end
                    if !haskey(read(file), "$station/$comp")
                        # write data
                        write(file, "$station/$comp", f.x[1])
                    end
                catch e
                    println(e)
                    println(seis)
                end
            end
        end
    end
end

################ Earthquake Event information ##############
dates  = [Date(2019, 11, 20), Date(2019, 11, 24)]
dt = [DateTime(2019, 11, 20, 4, 27, 5), DateTime(2019, 11, 24, 0, 54, 1)]
magnitudes = [6.3, 6.3]
lat = [13.886, 51.381]
lon = [-93.207, 175.512]
julian = ["2019_324", "2019_328"]

# create dataframe to make summarizing easier 
events = DataFrame(date = dates, datetime = dt, lat=lat, lon = lon, magnitude = magnitudes, julian = julian)
@eval @everywhere events = $events
############################################################



########## Download full waveforms for event days ###########
# download data 
paths = Array{String,1}(undef, 0)
for (ind, row) in enumerate(eachrow(events))
    try
        #path = join([yr,lpad("$(row.julian)", 3,"0")],"_")
        yr = split(row.julian, "_")[1]
        path = row.julian 
        println(path)
        push!(paths, path)
        node_path  = S3Path("s3://seisbasin/continuous_waveforms/$yr/$path/", config=aws)
        files = joinpath.("continuous_waveforms/$yr/$path/", 
                            convert.(String, readdir(node_path)))
        #files = [elt["Key"] for elt in collect(s3_list_objects(aws, "seisbasin", "continuous_waveforms/$(yr)/$(path)/", max_items=1000))]
        println(files)
        #files = [convert(String, elt) for elt in files]
        ec2download(aws, "seisbasin", files, "~/")
    catch e
        println(e)
    end
end
##########################################################


######### Clip waveforms, write to event file ############
if !isdir("sample_events")
    mkpath("sample_events")
end

T = @elapsed pmap(row -> get_event_data(row), eachrow(events))
##########################################################



################# Upload to S3 ###########################
event_files = glob("sample_events/*")

map(x -> s3_put(aws, "seisbasin", x, read(x), acl= "bucket_owner_full_control"), event_files)
println("Done")
##########################################################



####### Example file read ######## 
fid = h5open(event_files[1],"r")
A = read(fid)
keys(A)
keys(A["B1131"]["Z"]["waveform"])
##################################