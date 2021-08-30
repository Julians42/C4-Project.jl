# script for computing correlations for the BASIN project 

# you need to specify which year of correlations to run
# eg the command line to run this file should be something like:
# `julia basin_correlations.jl XXXX` where XXXX is the year, e.g. 2019

####################### Packages ########################
using Distributed, PyCall # add distributed computing

#Add all available cores
addprocs()

# send packages to all cores
@everywhere using SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWS, StructArrays, AWSS3, Statistics, JLD2, Glob, AbstractFFTs, HDF5

# add BASIN routines
@everywhere include("../src/BASIN.jl")
@everywhere using .BASIN
##########################################################



#################### Constants ##########################
cc_step, cc_len = 3600, 3600
maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9
half_win, water_level = 30, 0.01
aws = AWS.AWSConfig(region="us-west-2")
scedc = "scedc-pds"
basin = "seisbasin"
network = "CI"
channel1 = "BH?"
channel2 = "HH?"
OUTDIR = "/home/ubuntu/data"
rootdir = "/home/ubuntu/"
sources = ["TA2","LPC","CJM", "IPT", "SVD", "SNO", "DEV", "VINE", "ROPE", 
        "ARNO", "LUCI", "ROUF", "KUZD", "ALLI", "CHN", "USB", "Q0048", "Q0066",
        "PASC","RUS","FUL","SRN","BRE"]
all_stations = DataFrame(CSV.File("../../docs/updated_sources.csv"))

##########################################################



################## User selected job #####################
yr = parse(Int64, ARGS[1])
job_id = join(["_", yr]) # postfix for files 

# wrap all parameters in dictionary to pass collectively
params = Dict("aws" => aws, "cc_step" => cc_step, "cc_len" => cc_len, "maxlag" => maxlag, "yr" => yr,
        "fs" => fs, "half_win" => half_win, "water_level" => water_level, "sources" => sources,
        "all_stations" => all_stations, "fs" => fs, "rootdir" => rootdir, "network" => network,
        "OUTDIR" => OUTDIR, "freqmin" => freqmin, "freqmax" => freqmax)# Read in station locations and list source stations

# get date ranges based off user argument input
dates = nothing
if yr == 2017
    dates = collect(Date("2017-01-26"):Day(1):Date("2017-03-26"))
elseif yr == 2018
    dates = collect(Date("2018-07-20"):Day(1):Date("2018-08-27"))
elseif yr == 2019 # two separate data ranges for BASIN data in 2019
    dates = vcat(collect(Date("2019-05-11"):Day(1):Date("2019-06-16")), collect(Date("2019-11-10"):Day(1):Date("2019-12-15")))
else
    print("That is not a year when BASIN data was deployed. Try 2017, 2018, or 2019 as the input argument!")
    return 1
end
##########################################################



################## Correlate each day ####################
Tfull = @elapsed map(dayofyear -> correlate_day(dayofyear, params), dates)
##########################################################



############### Stack day correlations ###################
found_sources = unique([join(split(f,".")[1:2],".") for f in readdir("CORR/")]) # get files to save
T = @elapsed pmap(x -> stack_h5(x, job_id, params), found_sources)
##########################################################



############# Send files to S3 bucket ###################
stacked_files = glob("nodestack/*/*")
Transfer = @elapsed pmap(x -> s3_put(aws, "seisbasin", joinpath("BASIN", x), read(x)), stacked_files)
println("Finished")
##########################################################