# script for computing correlations for the BASIN project 

# you need to specify which year of correlations to run
# eg the command line to run this file should be something like:
# `julia basin_correlations.jl XXXX` where XXXX is the year, e.g. 2019

####################### Packages ########################
# using Distributed, PyCall # add distributed computing

#Add all available cores
# addprocs()

# send packages to all cores
# @everywhere using SeisIO, SeisNoise, Dates, CSV, DataFrames,SCEDC, AWS, StructArrays, AWSS3, Statistics, JLD2, Glob, AbstractFFTs, HDF5
using SeisIO, SeisNoise, Dates, CSV, DataFrames, StructArrays, Statistics, JLD2, Glob, AbstractFFTs, HDF5, PyCall
include("../src/wrapper_functions.jl")
include("../src/utils.jl")
##########################################################



#################### Constants ##########################
cc_step, cc_len = 3600, 3600
maxlag, fs = 300., 20. # maximum lag time in correlation, sampling frequency
freqmin, freqmax = 0.05, 9.9
half_win, water_level = 30, 0.01
network = "CI"
channel1 = "BH?"
channel2 = "HH?"
OUTDIR = "/mnt/DATA0/BASIN"   # where temporary data gets read
DATADIR = "/mnt/DATA3/BASIN"  # where data archive sit.
@assert OUTDIR != DATADIR "OUTDIR and DATADIR must be different"    # must be different so that outdir can be safely deleted at the end.
rootdir = "~/C4-Project.jl/"   # where the codes sit
sources = ["SVD"] 
all_stations = DataFrame(CSV.File("../../docs/updated_sources.csv"))

##########################################################



################## User selected job #####################
yr=2017
# wrap all parameters in dictionary to pass collectively
params = Dict("aws" => "local", "cc_step" => cc_step, "cc_len" => cc_len, "maxlag" => maxlag, "yr" => yr,
        "fs" => fs, "half_win" => half_win, "water_level" => water_level, "sources" => sources,
        "all_stations" => all_stations, "fs" => fs, "rootdir" => rootdir, "network" => network,
        "freqmin" => freqmin, "freqmax" => freqmax,"datadir"=>DATADIR,"outdir"=>OUTDIR)# Read in station locations and list source stations
dates = collect(Date("2017-01-26"):Day(1):Date("2017-03-26"))
##########################################################
## copied correlation_day 

dd=dates[1]
path = join([Dates.year(dd),lpad(Dates.dayofyear(dd),3,"0")],"_") # Yields "YEAR_JDY"
sources, maxlag, OUTDIR, yr, aws = params["sources"], params["maxlag"], params["outdir"], params["yr"], params["aws"]
filepath = "$(params["datadir"])/continuous_waveforms/$(yr)/$(path)/" 
println(filepath)
newdir = "$(params["outdir"])/continuous_waveforms/$(yr)/$(path)/"
allf = glob("continuous_waveforms/$yr/$path/*", "$OUTDIR") # if on local
broadbands = filter(x -> any(occursin.(sources, x)) && !occursin("Q0066",x), allf)
accelerometers = [f for f in allf if !any(occursin.(f, broadbands))]
## preprocessing
T_b = @elapsed map(f -> preprocess(f, false, params, path), broadbands)
T_a = @elapsed map(f -> preprocess(f, true, params, path), accelerometers)
println("Preprocessing Completed in $(T_b + T_a) seconds.")


##########
println(path)
fft_paths = glob("ffts/$path/*", "$OUTDIR")
sources = filter(f -> any(occursin.(sources, f)), fft_paths)
println(sources)
recievers = filter(f -> !any(occursin.(sources, f)), fft_paths)
print("There are $(length(recievers)) available for correlation.")
reciever_blocks = collect(Iterators.partition(recievers, convert(Int64, ceil(length(recievers)/nprocs())))) 
println("Now Correlating $(length(reciever_blocks)) correlation blocks!")
Tcorrelate = @elapsed map(rec_files -> correlate_block(sources, collect(rec_files), maxlag,params), reciever_blocks)
println("$(length(reciever_blocks)) blocks correlated in $Tcorrelate seconds for $path.")

################## Correlate each day ####################
# Tfull = @elapsed map(dayofyear -> correlate_day(dayofyear, params), dates)
##########################################################



############### Stack day correlations ###################
found_sources = unique([join(split(f,".")[1:2],".") for f in readdir("$OUTDIR/CORR/")]) # get files to save
T = @elapsed pmap(x -> stack_h5(x, job_id, params), found_sources)
##########################################################



############# Send files to S3 bucket ###################
stacked_files = glob("nodestack2/*/*","$OUTDIR")
println(stacked_files)
# Transfer = @elapsed pmap(x -> s3_put(aws, "seisbasin", joinpath("BASIN", x), read(x)), stacked_files)
println("Finished")
##########################################################
