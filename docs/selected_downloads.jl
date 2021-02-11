# Install packages, use multiple cores
using SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, JLD2, Statistics, Glob, AWSS3, AWSCore, SCEDC, Distributed
aws = aws_config(region="us-west-2");

OUTDIR = "~/Downloads" # CHANGE THIS 
#optional multithreading - advised if you're doing a lot of files! 
addprocs() # add processes
@everywhere begin
    using SCEDC, AWSS3, AWSCore
    aws = aws_config(region="us-west-2");
end

# params
# these are the stations we want and years we have data for
recievers = ["CI.CPP","CI.PSR","CI.WLT","CI.RIO","CI.PASC","CI.GSA","CI.CAC","CI.CLT","CI.PDU",
             "CI.PSR", "CI.PASA", "CI.MLS","CI.CLT","CI.PDU","CI.CHN"]
months = ["January","February","March","April","May","June","July","August",
            "September","October","November","December"]
years = [2017,2018, 2019]


# lets just list all the correlation files on s3 and filter the ones we want afterwards
all_seisbasin = Array{Array{String,1}, 1}(undef, 0)
for year in years
    for month in months
        files = [elt["Key"] for elt in collect(s3_list_objects(aws, "seisbasin", "corr_large/2017/April/"))]
        push!(all_seisbasin, files)
    end
end
all_seisbasin


# collapse and filter
new_seisbasin = Array{String,1}(undef, 0)
for list in all_seisbasin
    for elt in list
        push!(new_seisbasin, elt)
    end
end
filter!(x-> any(occursin.(recievers, x)), new_seisbasin)

@eval @everywhere new_seisbasin = $new_seisbasin
ec2download(aws, "seisbasin", new_seisbasin, OUTDIR)






