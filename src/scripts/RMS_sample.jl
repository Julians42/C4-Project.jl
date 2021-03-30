using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, JLD2, Statistics, Glob, 
            HDF5, AWSS3, AWS, Distributed, SCEDC, AWSCore

# multiple cores
addprocs()
@everywhere using SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, 
                    AWSS3, Statistics, JLD2, Glob, HDF5, AWS, AWSCore

# scp -i my_aws_key.pem /Users/julianschmitt/Downloads/updated_sources.csv ubuntu@ec2-54-202-138-141.us-west-2.compute.amazonaws.com:files/

all_stations = DataFrame(CSV.File("files/updated_sources.csv"))
@eval @everywhere all_stations = $all_stations
@everywhere begin
    function RMS_waveform(file::String)#, filt::Array{Float64, 1}=[0.05,1.])
        """ Compute RMS Waveform Metric """
        S = read_data("mseed", file)
        add_location(S, all_stations)
        demean!(S)
        a = S.x[1]
        SeisIO.taper!(S, t_max = 10., α= 0.05, N_min = 10)
        half = bandpass(S.x[1], 0.05, .5, 100.)
        one = bandpass(S.x[1], 0.05, 1., 100.)
        urban = bandpass(S.x[1], 1.5, 3., 100.)
        half_sum = sum(half.^2)
        one_sum = sum(one.^2)
        urban_sum = sum(urban.^2)
        st = std(S.x[1])
        npts = length(S.x[1])
        d = Dict("name"=>S.id[1], "network"=>split(S.id[1], ".")[1], "station" => split(S.id[1], ".")[2],
            "channel"=>"$(S.id[1][end])", "low_freq" => half_sum, "med_freq" => one_sum, 
                "urban" => urban_sum, "npts"=>npts, "datetime"=>DateTime(u2d(S.t[1][1,2] *10^(-6))), "location"=> S.loc[1], "std_raw"=>st)
        return d
    end

    # helper functions
    function LLE_geo(station, df)
        """ Find station matching location and return geoloc object"""
        try
            row = df[(findfirst(x -> x==station, df.station)),:]
            lat, lon, el = row.latitude[1], row.longitude[1], row.elevation[1]
            geo = GeoLoc(lat = float(lat), lon = float(lon), el = float(el))
            return geo
        catch 
            return nothing
        end
    end
    function add_location(s::SeisData,df::DataFrame)
        """ Adds locations to array of seisdata from a dataframe """
        name = split(s.name[1],".")[2]
        geo = LLE_geo(name, df)
        if !isnothing(geo)
            s.loc[1] = geo
        else
            println("Station $name can't be found in the dataframe")
        end
    end
end
# Select days and download data
#days = ["2017-02-05", "2017-03-19","2018-08-05","2019-05-12","2019-12-08"]

ss1, e1 = Date("2017-01-29"), Date("2017-03-15")
r2017 = ss1:Day(7):e1
ss2, e2 = Date("2018-07-22"), Date("2018-08-27")
r2018 = ss2:Day(7):e2
ss3, e3 = Date("2019-05-12"), Date("2019-06-16")
ss4, e4 = Date("2019-11-10"), Date("2019-12-15")
r2019 = vcat(collect(ss3:Day(7):e3), collect(ss4:Day(7):e4))
days = vcat(r2017, r2018, r2019)


aws = AWS.AWSConfig(region="us-west-2")
@eval @everywhere aws = $aws
# Notes
# 2017-02-05 had mag 3.3 just south of mexican border
# 2017-03-19 mag 2.6 near SLO
# 2019-05-12 nice and quiet - only mag ~2.6 ish off oregon coast
# put dictionaries in CSV to easily save
df = DataFrame(name= String[], network=String[], station=String[], channel=String[], low_freq=Float64[],
                med_freq=Float64[], urban=Float64[], npts=Int64[], datetime=DateTime[], location=GeoLoc[], std_raw=Float64[])
aws2 = AWSCore.aws_config(region="us-west-2") # having trouble with AWSS3 and SCEDC packages
@eval @everywhere aws2 = $aws2
for day in days
    d = Date(day)
    yr = Dates.year(d)
    path = join([yr,lpad(Dates.dayofyear(d),3,"0")],"_")

    # read files and download
    raw_waveforms = readdir(S3Path("s3://seisbasin/continuous_waveforms/$yr/$path/", config = aws))
    raw_waveforms = joinpath.("continuous_waveforms/$yr/$path/", convert.(String, raw_waveforms))
    ec2download(aws2, "seisbasin", raw_waveforms, "~/data")

    # get RMS for all files - returns an array of dictionaries
    println("Starting Processing of $(length(raw_waveforms)) files...")
    dicts = pmap(f -> RMS_waveform(f), joinpath.("data", raw_waveforms))
    
    map(d -> push!(df, d), dicts)
    println("Processing Completed!")

    # sadly have to clean up
    rm("data/continuous_waveforms/", recursive=true)
end


# save and upload to seisbasin
CSV.write("/home/ubuntu/RMS_values2.csv", df)
s3_put(aws2, "seisbasin", "RMS_values2.csv", read("RMS_values2.csv"))
println("Done")

# transfer file back to local
# scp -i my_aws_key.pem  ubuntu@ec2-172-31-16-183.us-west-2.compute.amazonaws.com:RMS_values2.csv /Users/julianschmitt/Downloads/RMS_values2.csv



ec2download(aws,"seisbasin", raw_waveforms[1:3], "/Users/julianschmitt/Downloads/")


