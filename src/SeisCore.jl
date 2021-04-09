module SeisCore

# packages
using Base, Core, SeisIO, SeisNoise, Dates, CSV, DataFrames, SCEDC, AWS, Distributed, AWSS3, Glob,
        Statistics, JLD2, HDF5, AbstractFFTs, Plots

# functions
include("download.jl")
include("process.jl")
include("correlate.jl")
include("big_correlate.jl")
end