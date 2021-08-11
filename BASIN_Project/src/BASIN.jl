module BASIN
""" Module for the BASIN Project """

using Distributed, SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWSCore, StructArrays, AWSS3, Statistics, JLD2, Glob, AbstractFFTs, HDF5
include("helper_routines.jl")
include("wrapper_functions.jl")

end 