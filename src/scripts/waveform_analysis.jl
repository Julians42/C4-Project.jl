using Pkg, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWS, Distributed, AWSS3, Glob, HDF5, Statistics

ENV["GR"] = ""
Pkg.build("GR")

# ColorScheme
# σs = 0.5:0.1:2
# normal_x = -5:0.01:5
# normal_y = [exp.(-normal_x.^2 / (2σ^2)) / (2π * σ^2) for σ in σs];
# loadcolorscheme(:cm_maxamp,ColorSchemes.gist_heat.colors[end-30:-1:1], "maxamp color", "for waveform plot");

addprocs()
@everywhere begin
    using Pkg; Pkg.activate("ubuntu")
    using AWS, AWSS3, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, Glob, HDF5
end
@everywhere begin 
    aws = AWS.AWSConfig(region="us-west-2")
end
@everywhere begin 
    function plot_waveforms2(row)
        # get files for event
        for elt in ["EHE","EHN","EHZ"]
            raw_seis = glob("continuous_waveforms/*/*/*")
            filter!(x -> occursin(elt, x), raw_seis)
            key = join(["2019","$(row.julian)"])
            samp_day = filter(x -> occursin(key,x), raw_seis)
            startsplit, endsplit = row.datetime, row.datetime+Dates.Minute(3)

            # load, clip and sort waveforms
            ar = [tsplit(read_data("mseed", data), startsplit, endsplit) for data in samp_day]
            distances = [get_dist(s.loc[1], GeoLoc(lat=row.lat, lon=row.lon)) for s in ar]
            p = sortperm(distances)
            ar, distances = ar[p], distances[p] # reorder
            println(length(distances))
            # plot parameters - nasty stuff
            lmin, lmax = minimum(distances), maximum(distances)
            y_labels = (round.(collect(range(0, 5*length(ar), length=10)),digits=3),round.(collect(range(lmin, lmax, length=10)),digits=3))
            NS_plot = Plots.plot(ylims = (-5,length(ar)*5+5), 
                        yticks=y_labels,xlabel = "Time (s)", ylabel = "Station Nodes by Distance from Event",  
                        xtickfontsize=5,ytickfontsize=5,fontsize=5, xguidefontsize = 8, yguidefontsize = 10,
                        legendfontsize = 8, title="$elt Waveforms for $(row.date) event (mag: $(row.magnitude))",
                        titlefontsize= 8)
            
            ar_list, max_amp = [], []
            println(typeof(ar))
            for (ind, s) in enumerate(ar)
                # get distance band
                if count(s.x[1] .>1)<100 || count(s.x[1] .== 0) < .1*length(s.x[1])
                    println(lmax, "  ", lmin)
                    println(distances[ind])
                    push!(max_amp, maximum(abs.(s.x[1])))
                    shift_factor = ((distances[ind]-lmin)/(lmax-lmin))*length(ar)*5 # moveout by distance
                    println("shift factor:", shift_factor)
                    #data = [x < 1 ? x : 0 for x in s.x[1]] .* 20 .+ shift_factor
                    m, st = mean(s.x[1]), std(s.x[1])
                    println(m, "  ",st)
                    #data = [(x < m +2*st) && (x > m-2*st) ? x : (m+2*st) for x in (s.x[1])] .* 20 .+ shift_factor
                    data = s.x[1] .* 20 .+ shift_factor
                    push!(ar_list, data)
                end
            end
            println("$(length(ar_list)) waveforms to be plotted")
            # plotting
            T = collect(0.:1/ar[1].fs[1]:180.)
            Plots.plot!(NS_plot, T, ar_list, color=:cm_maxamp, colorbar_title="Maximum Amplitude", 
                    line_z=max_amp', fmt = :png, linewidth = 0.5, reuse = false, legend = false)    
            Plots.plot!(size=(250,600),dpi=500)
            filepath = joinpath(rootdir,
                "waveform_analysis/B1_$(key)_$(elt)_0.1_2.0.png")
            # ensure filepath is valid 
            DIR = dirname(filepath)
            if !isdir(DIR)
                mkpath(DIR)
            end
            println(filepath)
            png(NS_plot,filepath)
        end
    end
    function LLE_geo(s::SeisData, df)
        """ Add location to SeisData object from dataframe"""
        network, station = split(s.id[1], ".")[1], split(s.id[1],".")[2]
        try
            row = filter(row -> (row.network==network) & (row.station==station), df)[1,:]
            s.loc[1] = GeoLoc(lat = float(row.latitude), lon = float(row.longitude), el = float(row.elevation))
        catch 
            return nothing # when station is not found in dataframe
        end
    end
    function tsplit(S::SeisData, s::DateTime, t::DateTime)
        LLE_geo(S, locations)
        sync!(S, s=s, t=t)
        detrend!(S)
        taper!(S)
        filtfilt!(S, fl=0.05, fh=10., np=1)
        return S
    end
end
# get locations of nodes etc
#locations = DataFrame(CSV.File("files/modified_nodal.csv"))
rootdir=""

# create event date catelog 
# dates = [Date(2019, 11,12), Date(2019, 11, 13), Date(2019, 12,11), Date(2019, 12, 5), Date(2019, 11,29)]
# dt = [DateTime(2019, 11,12, 2,13,52,180), DateTime(2019, 11, 13, 19,26,53,240), DateTime(2019, 12,11,2,48,19,960),
#         DateTime(2019,12,5,8,55,31,650),DateTime(2019,11,29,0,49,42,870)]
# magnitudes = [3.97, 3.92, 3.87, 3.79, 3.72]
# lat = [32.79133, 35.60517, 31.70350, 35.69667, 35.71350]
# lon = [-115.54733, -117.40750,-116.15600, -117.61167,-117.56867]
# julian = [convert(Int64, floor(datetime2julian(dt2)-datetime2julian(DateTime(2019))))+1 for dt2 in dt]
# println(dates)
# events = DataFrame(date = dates, datetime = dt, lat=lat, lon = lon, magnitude = magnitudes, julian = julian)

# dates = [Date(2017, 2,10)]
# dt = [DateTime(2017, 2,10, 10,10,32)]
# magnitudes = [3.4]
# lat = [34.131]
# lon = [-117.046]
dates = [Date(2017, 2, 24), Date(2018, 8, 19), Date(2019, 5, 26), Date(2019, 12, 15)]
dt = [DateTime(2017, 2,24, 17, 38, 44), DateTime(2018, 8, 19, 0, 19, 40), DateTime(2019, 5, 26, 7,41,15), 
        DateTime(2019, 12, 15, 6, 11, 51)]
#julian = [join([Dates.year(dt2),lpad(convert(Int64, floor(datetime2julian(dt2)-datetime2julian(DateTime(Dates.year(dt2))))+1], "_"), 3, "0") for dt2 in dt]
julian = ["2017_055", "2018_231","2019_146", "2019_349"]
println(dates)
lat = [-23.259, -18.133,-5.812, 6.697]
lon = [-178.804, -178.153, -75.270, 125.174]
magnitudes = [6.9, 8.2, 8.0, 6.8]
events = DataFrame(date = dates, datetime = dt, lat=lat, lon = lon, magnitude = magnitudes, julian = julian)
# download data
# yr = "2017"
# @eval @everywhere locations = $locations
@eval @everywhere events = $events




########### Plotting Script ############
 
# for (ind, row) in enumerate(eachrow(events))

#     path = join(["2019","$(row.julian)"],"_")
#     files = [elt["Key"] for elt in collect(s3_list_objects(aws, "seisbasin", "continuous_waveforms/$(yr)/$(path)/", max_items=1000))]
#     ec2download(aws, "seisbasin", files, "~/")
#     plot_waveforms2(row)
# end

# for (ind, row) in enumerate(eachrow(events))
#     row1 = events[1,:]
#     path = join(["2019","$(row1.julian)"],"_")
#     #files = [elt["Key"] for elt in collect(s3_list_objects(aws, "seisbasin", "continuous_waveforms/$(yr)/$(path)/", max_items=1000))]
#     #ec2download(aws, "seisbasin", files, "~/")
#     println(row1)
#     plot_waveforms2(row1)
# end
# pmap(x -> plot_waveforms2(x), eachrow(events))

# file transfer 
#scp -i my_aws_key.pem /Users/julianschmitt/Documents/Schoolwork/Seismology/SeisCore.jl/docs/modified_nodal.csv ubuntu@ec2-34-223-59-140.us-west-2.compute.amazonaws.com:files/
# scp -i my_aws_key.pem ubuntu@ec2-34-223-59-140.us-west-2.compute.amazonaws.com:waveform_analysis/ /Users/julianschmitt/Downloads/waveform_analysis/ --recursive
####################################



############ Extract event raw waveforms and upload to aws for Laura! #############
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
locations = DataFrame(CSV.File("files/updated_sources.csv"))

@eval @everywhere locations = $locations
# read files, split and write to h5
@everywhere begin
    function get_event_data(row)
        yr = split(row.julian, "_")[1]
        files = glob("continuous_waveforms/$yr/$(row.julian)/*")
        println(files)
        startsplit, endsplit = row.datetime, row.datetime+Dates.Minute(180)
        h5open("sample_events/$(row.julian).h5", "cw") do file
            if !haskey(read(file), "meta")
                write(file, "meta/date", "$(row.date)")
                write(file, "meta/datetime", "$(row.datetime)")
                write(file, "meta/lat", row.lat)
                write(file, "meta/lon", row.lon)
                write(file, "meta/magnitude", row.magnitude)
                write(file, "meta/julian", row.julian)
            end
            for seis in files
                try
                    # grab data, split, and add location
                    f = tsplit(read_data("mseed", seis), startsplit, endsplit) # split from t=0 to 10 minutes after
                    network, station = split(f.id[1], ".")[1], split(f.id[1],".")[2]
                    comp = string(f.id[1][end])
                    LLE_geo(f, locations) # add location 
                    dist = get_dist(f.loc[1], GeoLoc(lat=row.lat, lon=row.lon)) # get distance
                    #if !haskey(read(file), "$station/meta") && !haskey(read(file), "$station/meta/dist")
                    try
                        # write metadata
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
if !isdir("sample_events")
    mkpath("sample_events")
end

T = @elapsed pmap(row -> get_event_data(row), eachrow(events))


# upload to s3 
event_files = glob("sample_events/*")

map(x -> s3_put(aws, "seisbasin", x, read(x), acl= "bucket_owner_full_access"), event_files)
println("Done")


# test file
fid = h5open(event_files[1],"r")
A = read(fid)
keys(A)
keys(A["B1131"]["Z"]["waveform"])