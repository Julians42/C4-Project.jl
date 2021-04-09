# Script for docker container to run correlations on AWS Batch
include("SeisCore.jl")
using .SeisCore, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWS, Distributed, JLD2, Glob, AWSS3, HDF5, Statistics


addprocs()
# add Constants and locations dataframe
@everywhere begin 
    using .SeisCore, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWS, JLD2, Glob, AWSS3, HDF5, Statistics
    aws = AWS.AWSConfig(region="us-west-2")
    rootdir = ""
    network = "CI"
    channel1 = "BH?"
    channel2 = "HH?"
    OUTDIR = "~/data"
    cc_step, cc_len = 3600, 3600
    maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
    freqmin, freqmax = 0.05, 9.9
    half_win, water_level = 30, 0.01
    all_stations = DataFrame(CSV.File("locations.csv"))
end


arg = ENV["AWS_BATCH_JOB_ARRAY_INDEX"]
startdate = Date(2000)+Month(arg)
enddate = startdate+Month(1)-Day(1)
days = startdate:Day(1):enddate
println("Processing download for: ", startdate, " to ",enddate)

# Big wrapper function for processing and correlating everything
Tcbig = @elapsed map(dd -> correlate_big(dd, startdate), days)


# Stack each of the 100, 20, and 1 Hz data
stack_all()


# transfer to s3 
autos = glob("autocorrelations/*/*")
bigcorrs = glob("correlations/*/*")

pmap(x -> s3_put(aws, "seisbasin", x, read(x)), autos)
pmap(x -> s3_put(aws, "seisbasin", x, read(x)), bigcorrs)

println("Finished Correlations and upload for: ", startdate, " to ",enddate)
