using Pkg

Pkg.add(["AWS", "AWSS3", "AbstractFFTs", "CSV", "DataFrames", "Dates", "Distributed",
            "Glob", "HDF5", "JLD2", "Parallelism", "Plots", "PyCall", "SeisIO",
            "SeisNoise", "Serialization", "Statistics", "StructArrays"])

# add custom       
using Pkg; Pkg.add(PackageSpec(url="https://github.com/tclements/SCEDC.jl", rev="master"))