module BASIN
""" Module for the BASIN Project """

using Distributed, SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWSCore, StructArrays, AWSS3, Statistics, JLD2, Glob, AbstractFFTs, HDF5
include("utils.jl")
include("wrapper_functions.jl")

end 