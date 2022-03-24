module BASIN_local
""" Module for the BASIN Project """

using Distributed, SeisIO, SeisNoise, Dates, CSV, DataFrames, StructArrays, Statistics, JLD2, Glob, AbstractFFTs, HDF5
include("utils.jl")
include("wrapper_functions.jl")

end 