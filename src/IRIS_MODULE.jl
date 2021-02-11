include("SeisCore.jl")
T=@elapsed using .SeisCore, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Glob, AWSS3
println(T)

# Constants
#rootdir = "/Users/julianschmitt/Downloads/Seismo/IRIS/" # set to "" if on EC2/ in container, but can run locally 
rootdir = ""
arg = ENV["AWS_BATCH_JOB_ARRAY_INDEX"] # get date range from environment variable
startdate = Date(2000)+Month(arg)
enddate = startdate+Month(1)-Day(1)
days = startdate:Day(1):enddate
println("Processing download for: ", startdate, " to ",enddate)

data_sources = [["*.*.*.HH*","NCEDC"],["*.*.*.BH*","NCEDC"],["*.*.*.HH*","IRIS"],["*.*.*.BH*","IRIS"]]

# get seismic data for all days
map(x -> SeisCore.get_seisdata(x, data_sources, rootdir), days)

println("Finished: ", startdate, " to ", enddate)