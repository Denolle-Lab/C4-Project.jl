export preprocess, correlate_block, correlate_day, stack_h5

####################### Preprocessing routine ##########################
function preprocess(file::String,  accelerometer::Bool=false, params::Dict=params, path::String=path)
    """
        Load raw seisdata file and process, saving to fft per frequency
    this is exactly the same as preprocess2 from SeisCore.jl
    """
    # unpack needed params
    rootdir, OUTDIR, samp_rate, all_stations = params["rootdir"], params["outdir"], params["fs"], params["all_stations"]
    freqmin, freqmax, cc_step, cc_len = params["freqmin"], params["freqmax"], params["cc_step"], params["cc_len"]
    half_win, water_level = params["half_win"], params["water_level"]

    try
        data = SeisData()
        try # assume mseed
            read_data!(data, "mseed", file)
        catch
            read_data!(data, "seisio", file)
        end
        gaps = size(data.t[1])[1] # N-2 gaps (eg gaps = 12 tests for 10 gaps)
        pts = size(data.x[1])[1]
        fs_temp = data.fs[1]
        windows = u2d.(SeisIO.t_win(data.t[1], data.fs[1]) .* 1e-6)
        startend = windows[:,2] .- windows[:,1] .> Second(cc_len)
        if gaps < 25 && any(startend) # If few gaps and sufficient (2/3) data present, data is clean
            try
                add_location(data, all_stations)
                S = SeisData()
                try
                    S = process_raw(data, samp_rate)
                catch # trying to sample non 100Hz data at 100 Hz - skip resampling 
                    S = process_raw(data, data.fs[1])
                end

                # truncate the data to night time (first half of the UTC time for LA stations)
                # this starts at 5pm. Move to 8pm.!!!! Here i do not change the start time!!!!
                sync!(S,s=u2d(S.t[1][1,2]*1e-6)+Dates.Hour(3),t=u2d(S.t[1][1,2]*1e-6)+Dates.Hour(12))
                ### Comment out above to keep all of the data
                R = RawData(S,cc_len,cc_step)
                SeisNoise.detrend!(R)
                bandpass!(R,freqmin,freqmax,zerophase=true)
                SeisNoise.taper!(R)
                FFT = nothing
                if accelerometer # accelerometer - so we integrate
                    FFT = rfft_raw(R,1)
                else # then its a seismometer
                    FFT = compute_fft(R)
                end
                coherence!(FFT,half_win, water_level) # try no coherence ??
                try # save fft 
                    crap=joinpath(OUTDIR, "ffts")
                    if isdir(crap)==false
                        mkpath(crap)
                    end
                    root_fft = joinpath(OUTDIR, "ffts/$path/")
                    save_fft(FFT, root_fft)
                catch e
                    println(e)
                end
                return data.id[1]
                println("Successfully processed $(data.id[1])")
            catch e
                println(e)
            end
        end
    catch e 
        println(e)
    end
end

############## Correlate two blocks of day waveforms ###################
function correlate_block(src::Array{String,1}, rec::Array{String,1}, maxlag::Float64, params::Dict=params)
    """ 
        Correlation function for pair of fft data
        - noise filter and stacking
        - saves directly to disk: ensure directory is correct if not on AWS
    """
    # load sources and recievers
    crap=joinpath("$(params["outdir"])/CORR")
    if isdir(crap)==false
        mkpath(crap)
    end
    sources = map(x-> load_fft(x, string(x[end-7:end-5])), src)
    receivers = map(x-> load_fft(x, string(x[end-7:end-5])), rec)
    for source in sources
        for receiver in receivers
            try
                C = correlate(source,receiver,maxlag)
                # cc_medianmute!(C, 10.) # remove correlation windows with high noise 
                stack!(C)
                pair, comp = name_corr(C), C.comp
                save_named_corr(C,"$(params["outdir"])/CORR/$pair/$comp")
            catch e
                # add counter - iterate num bad print
                println(e)
            end
        end
    end
end

#################### Downloads and Correlates day of data ###########################
function correlate_day(dd::Date, params::Dict=params)
    """ Wrapper function for daily correlations 
    """
    path = join([Dates.year(dd),lpad(Dates.dayofyear(dd),3,"0")],"_") # Yields "YEAR_JDY"
    @eval @everywhere path = $path

    # unpack needed params
    sources, maxlag, OUTDIR, yr, aws = params["sources"], params["maxlag"], params["outdir"], params["yr"], params["aws"]

    # here download the data from S3
    if aws!="local"
        # get filepaths for SCEDC source stations 
        scedc_files = get_scedc_files(dd, params)
        filter!(x -> any(occursin.(sources, x)), scedc_files)

        # filepaths for nodes
        filelist_b = S3Path("s3://seisbasin/continuous_waveforms/$(yr)/$(path)/", config=aws) # use new AWS functions
        filelist_basin = joinpath.("continuous_waveforms/$(yr)/$(path)/", 
                                    convert.(String, readdir(filelist_b))) # add directory 
        println("There are $(length(filelist_basin)) node files available for $path")
        # download scedc and seisbasin data
        try
            ec2download(aws, "scedc-pds", scedc_files, OUTDIR)
            ec2download(aws, "seisbasin", filelist_basin, OUTDIR)
            println("Download complete!")
        catch e
            println(e)
            println("Failed to process for $path. Continuing to next day.")
            return 1
        end
    else  # Local / workstation use

          # filepaths for nodes
        filepath = "$(params["datadir"])/continuous_waveforms/$(yr)/$(path)/"
        # convert.(String, readdir("$(params["datadir"])/continuous_waveforms/$(yr)/$(path)/"))) # add directory 
        println(filepath)
        newdir = "$(params["outdir"])/continuous_waveforms/$(yr)/$(path)/"
        # THIS IS NOT GREAT
        # chown("ffts",39101,39123) # PLEASE EDIT HERE!!! First integer is the user ID, the second is the Denolle lab group number
        println(newdir, )
        # cp(filepath, newdir,force=true)
        filelist_basin = readdir(newdir)
        println(filelist_basin)

    end
    # preprocess data
    if params["aws"]!="local"
    allf = glob("continuous_waveforms/$yr/$path/*", "/home/ubuntu/data/") # if on aws
    else
    allf = glob("continuous_waveforms/$yr/$path/*", "$OUTDIR") # if on local
    end
    broadbands = filter(x -> !occursin("Q0066",x), allf)
    println("printing the filenames of the velocitmeters")
    println(broadbands)
    accelerometers = filter(x -> any(occursin("Q0066",x)), allf)

    T_b = @elapsed pmap(f -> preprocess(f, false, params, path), broadbands)
    T_a = @elapsed pmap(f -> preprocess(f, true, params, path), accelerometers)
    println("Preprocessing Completed in $(T_b + T_a) seconds.")
    println(path)

    # get indices for block correlation
    if params["aws"]!="local"
    fft_paths = glob("ffts/$path/*", "/home/ubuntu/")
    else
    fft_paths = glob("ffts/$path/*", "$OUTDIR")
    end
    sources = filter(f -> any(occursin.(sources, f)), fft_paths)
    receivers = filter(f -> !any(occursin.(sources, f)), fft_paths)
    print("There are $(length(receivers)) available for correlation.")
    if length(receivers) == 0 # if no ffts available for that day 
        return 1
    else # then data available: Correlate!
        receiver_blocks = collect(Iterators.partition(receivers, convert(Int64, ceil(length(receivers)/nprocs())))) 
        # println("Now Correlating $(length(receiver_blocks)) correlation blocks!")
        Tcorrelate = @elapsed pmap(rec_files -> correlate_block(sources, collect(rec_files), maxlag,params), receiver_blocks)
        println("$(length(receiver_blocks)) blocks correlated in $Tcorrelate seconds for $path.")
    end
    if params["aws"]!="local" # if on aws, remove data from EC2 instance
    rm("/home/ubuntu/data/continuous_waveforms", recursive=true) # cleanup raw data
    else
    rm("$OUTDIR/ffts/$path", recursive=true) # cleanup raw data on the SSD working disk
    # rm("$OUTDIR/continuous_waveforms/$yr/$path", recursive=true) # cleanup raw data on the SSD working disk
    end
end

######################## Stacking Routine ############################
function stack_h5(tf::String, postfix::String, params::Dict=params)
    """Stacks all files for source-receiver pair, saves to h5 file by source"""
    
    # extract needed parameters from dict
    yr, all_stations,OUTDIR = params["yr"], params["all_stations"],params["outdir"]
    # scrape correlation filenames - get unique sources and receiver combinations
    from_s = glob("CORR/$tf*/*/*","$OUTDIR")
    found_receivers = unique([join(split(f,".")[4:5],".") for f in from_s])
    components = ["EE","EN","EZ", "NE", "NN","NZ", "ZE", "ZN", "ZZ"]
    fname = join([tf, postfix, ".h5"])
    filename = joinpath("$OUTDIR/nodestack/$yr", fname)
    if !isdir(dirname(filename))
        mkpath(dirname(filename))
    end
    # create file and open
    h5open(filename, "cw") do file
        source_station = split(tf,".")[2]
        source_loc = LLE_geo(source_station, all_stations)
        samplef = glob("CORR/$tf*/ZZ/*","$OUTDIR")[1]
        C = load_corr(samplef, "ZZ")
        # add metadata information about correlation processing and source location 
        if !haskey(read(file), "meta")
            write(file, "meta/corr_type", C.corr_type)
            write(file, "meta/cc_len", C.cc_len)
            write(file, "meta/cc_step", C.cc_step)
            write(file, "meta/whitened", C.whitened)
            write(file, "meta/time_norm", C.time_norm)
            write(file, "meta/notes", C.notes)
            write(file, "meta/maxlag", C.maxlag)
            write(file, "meta/starttime", Dates.format(u2d(C.t[1]), "yyyy-mm-dd HH:MM:SS"))
            write(file, "meta/samp_freq", C.fs)
            println(source_loc.lon)
            if !isnothing(source_loc)
                write(file, "meta/lon", source_loc.lon)
                write(file, "meta/lat", source_loc.lat)
                write(file, "meta/el", source_loc.el)
            end
        end
        # iterate through node receivers 
        for rec in found_receivers
            try
                # sample_r = glob("CORR/$tf*$rec*/ZZ/*")[1]
                # Cr = load_corr(sample_r, "ZZ")
                rec_station = split(rec, ".")[2]
                rec_loc = LLE_geo(rec_station, all_stations)
                if !haskey(read(file), "$rec/meta")
                    write(file, "$rec/meta/lon", rec_loc.lon)
                    write(file, "$rec/meta/lat", rec_loc.lat)
                    write(file, "$rec/meta/el", rec_loc.el)
                    write(file, "$rec/meta/dist", get_dist(source_loc, rec_loc))
                    write(file, "$rec/meta/azi", get_azi(source_loc, rec_loc))
                    write(file, "$rec/meta/baz", get_baz(source_loc, rec_loc))
                end
                for comp in components
                    try
                        # load correlations for this receiver by component 
                        files = glob("CORR/$tf*$rec*/$comp/*.jld2","$OUTDIR")
                        if length(files) == 0; continue; end
                        corrs = [load_corr(f, comp) for f in files]
                        corrs = [c for c in corrs if typeof(c) == CorrData]
        
                        # implement linear and phase-weighted stack methods. 
                        corr_sum = sum_corrs(corrs) # combine correlations into single CorrData object 
                        if isnothing(corr_sum); continue; end # if corr_sum fails to return, skip component

                        dims = size(corr_sum.corr)
                        # cc_medianmute2!(corr_sum, 3., false) # skip metadata copy which errors on sum_corrs
                        println("$(length(corrs)) day correlations. Dims: $dims. ",
                                "Threw out $(dims[2]-size(corr_sum.corr)[2]) day windows.")

                        # then stack
                        corr_mean = SeisNoise.stack(corr_sum, allstack=true, stacktype=mean)
                        corr_pws = SeisNoise.stack(corr_sum, allstack=true, stacktype=pws)

                        # save to file 
                        write(file, "$rec/$comp/linear", corr_mean.corr[:])
                        write(file, "$rec/$comp/pws", corr_pws.corr[:])
                    catch e 
                        println("Error in stacking $rec $comp: ", e)
                    end
                end
            catch e
                println(e)
            end
        end
    end
end
