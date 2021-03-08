module SeisCore

# packages
using Base, Core, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWSCore, Distributed, AWSS3, Glob, HDF5,
        Statistics, AbstractFFTs, JLD2, Plots

# variables
aws = aws_config(region="us-west-2")
network = "CI"
channel1 = "BH?"
channel2 = "HH?"
OUTDIR = "~/data"
rootdir = "" # manipulate for easy filepathing

# functions
include("download.jl")
include("processing.jl")
end