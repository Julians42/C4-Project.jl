include("SeisCore.jl")
T=@elapsed using .SeisCore, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, JLD2, Glob, AWSS3
println(T)

# Variables
#rootdir = "/Users/julianschmitt/Downloads/Seismo/IRIS/" # set to "" if on EC2/ in container, but can run locally 
#locations = DataFrame(CSV.File("../docs/updated_sources.csv")) # include locations dataframe 

aws = AWS.AWSConfig(region="us-west-2")
rootdir = ""
network = "CI"
channel1 = "BH?"
channel2 = "HH?"
OUTDIR = "~/data"

# get date range from environment variable
arg = ENV["AWS_BATCH_JOB_ARRAY_INDEX"]
startdate = Date(2004)+Month(arg)
enddate = startdate+Month(1)-Day(1)
days = startdate:Day(1):enddate
println("Processing download for: ", startdate, " to ",enddate)

data_sources = [["*.*.*.HH*","NCEDC"],["*.*.*.BH*","NCEDC"],["*.*.*.HH*","IRIS"],["*.*.*.BH*","IRIS"]]

# get seismic data for all days
map(x -> SeisCore.get_seisdata(x, aws, data_sources, rootdir), days)

println("Finished: ", startdate, " to ", enddate)