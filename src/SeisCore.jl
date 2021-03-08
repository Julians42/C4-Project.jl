module SeisCore

# packages
using Base, Core, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, AWSS3, Glob,
        Statistics, JLD2, HDF5, AbstractFFTs, Plots

# variables
aws = aws_config(region="us-west-2")
network = "CI"
channel1 = "BH?"
channel2 = "HH?"
OUTDIR = "~/data"
rootdir = "" # manipulate for easy filepathing
#locations = DataFrame(CSV.File("../docs/updated_sources.csv")) # include locations dataframe 

# functions
include("download.jl")
include("process.jl")
include("correlate.jl")
end