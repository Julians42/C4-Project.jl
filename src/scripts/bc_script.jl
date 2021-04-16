# Script for docker container to run correlations on AWS Batch
using SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWS, Distributed, JLD2, Glob, AWSS3, HDF5, Statistics, AbstractFFTs


addprocs()
@everywhere include("/root/SeisCore.jl")
@everywhere using .SeisCore
# add Constants and locations dataframe
#@everywhere include("/home/ubuntu/SeisCore.jl/SeisCore.jl")
@everywhere begin 
    using SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWS, JLD2, Glob, AWSS3, HDF5, Statistics, AbstractFFTs
    aws = AWS.AWSConfig(region="us-west-2")
    rootdir = ""
    network = "CI"
    channel1 = "BH?"
    channel2 = "HH?"
    OUTDIR = "~/data"
    num_procs = nprocs()
    cc_step, cc_len = 3600, 3600
    maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
    freqmin, freqmax = 0.05, 9.9
    half_win, water_level = 30, 0.01
    samp_rates = [1., 20., 100.] # for processing
    #all_stations = DataFrame(CSV.File("/home/ubuntu/SeisCore.jl/docs/updated_sources.csv"))
    all_stations = DataFrame(CSV.File("/root/updated_sources.csv"))
    params = Dict("aws" => aws, "cc_step" => cc_step, "cc_len" => cc_len, "maxlag" => maxlag,
            "fs" => fs, "half_win" => half_win, "water_level" => water_level, 
            "all_stations" => all_stations, "samp_rates" => samp_rates, "rootdir" => rootdir,
            "OUTDIR" => OUTDIR, "num_procs"=> num_procs, "freqmin" => freqmin, "freqmax" => freqmax)
end


arg = ENV["AWS_BATCH_JOB_ARRAY_INDEX"]
startdate = Date(2004)+Month(arg)
enddate = startdate+Month(1)-Day(1)
days = startdate:Day(1):enddate
println("Processing download for: ", startdate, " to ",enddate)

# Big wrapper function for processing and correlating everything
Tcbig = @elapsed map(dd -> correlate_big(dd, startdate, params), days[1:3])


# Stack each of the 100, 20, and 1 Hz data
stack_all()


# transfer to s3
@everywhere begin
    autos = glob("autocorrelations/*/*")
    bigcorrs = glob("correlations/*/*")
end 

pmap(x -> s3_put(aws, "seisbasin", x, read(x), acl="bucket-owner-full-control"), autos)
pmap(x -> s3_put(aws, "seisbasin", x, read(x), acl="bucket-owner-full-control"), bigcorrs)

println("Finished Correlations and upload for: ", startdate, " to ",enddate)
