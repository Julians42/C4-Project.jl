function correlate_big(dd::Date, startdate::Date = startdate, params::Dict = params)
    """ Wrapper Function: Computes autocorrelations and correlations for a specific day"""
    try # Last resort failsafe
        println("Starting correlations for $dd.")
        aws, rootdir = params["aws"], params["rootdir"]
        yr = Dates.year(dd)
        path = join([Dates.year(dd),lpad(Dates.dayofyear(dd),3,"0")],"_") # Yeilds "YEAR_JDY"
        #@eval @everywhere path = $path

        ###### GET SCEDC FILES #######
        # get filepaths for source stations - would be faster if they're available
        scedc_files = nothing

        if !isdir(joinpath(rootdir, "scedc_path/$yr/")); mkpath(joinpath(rootdir, "scedc_path/$yr/")); end
        try # most filelists are stored on seisbasin
            s3_get_file(aws, "seisbasin", "scedc_path/$yr/$path.csv", joinpath(rootdir, "scedc_path/$yr/$path.csv"))
            scedc_files = DataFrame(CSV.File(joinpath(rootdir, "scedc_path/$yr/$path.csv"))).Path
        catch e # in case file not found/present we use SCEDC lookup
            println(e)
            #scedc_files = get_scedc_files(dd, aws)
        end
        if isnothing(scedc_files); return 1; end # not worth computing if no SCEDC files - this should never happen
        ##############################
        #scedc_files = get_scedc_files(dd, aws)

        # filepaths for iris and ncedc data
        iris_path  = S3Path("s3://seisbasin/iris_waveforms/$yr/$path/", config=aws)
        ncedc_path = S3Path("s3://seisbasin/ncedc_waveforms/$yr/$path/", config=aws)

        # format file strings for download 
        iris_query, ncedc_query = convert.(String, readdir(iris_path)), convert.(String, readdir(ncedc_path))
        filelist_basin = vcat(joinpath.("iris_waveforms/$yr/$path/", iris_query), joinpath.("ncedc_waveforms/$yr/$path/", ncedc_query))

        # print pre-download sumamry
        println("There are $(length(filelist_basin)) node files and $(length(scedc_files)) SCEDC files available for $path.")


        # download scedc and seisbasin (iris/ncedc) data
        ec2download(aws, "scedc-pds", scedc_files, joinpath(rootdir,"data"))
        ec2download(aws, "seisbasin", filelist_basin, joinpath(rootdir,"data"))
        println("Download complete!")

        # preprocess - broadbands and seismometers are processed separately
        raw_waveforms = glob("data/*/$yr/$path/*", rootdir)# params["rootdir"])
        println("$(length(raw_waveforms)) waveforms available for processing on $dd.")

        T_b = @elapsed serial_names = robust_pmap(f -> preprocess(f, params), raw_waveforms)
        println(serial_names[1:3])
        println(readdir())
        println(glob("ffts/*", "/scratch"))
        println("Preprocessing completed in $T_b seconds. $(length(serial_names)) ffts computed for $dd.")

        serial_names = collect(Iterators.flatten(filter(x -> !isnothing(x), serial_names))) # filter NaN arrays
        serial_names = convert(Array{String, 1}, filter(x -> !isnothing(x), serial_names)) # filter NaN elements
        serial_names = joinpath.(joinpath(params["rootdir"], "ffts/"), serial_names)

        # autocorrelations    
        fft_list_100 = filter(name -> occursin("100_", name), serial_names) # find 100Hz data
        fft_list_100 = convert(Array{String, 1}, fft_list_100) # convert to string arrray
        println("FFT list 100 is $(length(fft_list_100)) stations long for $dd.") # print summary
        fft_100_stations = unique([join(split(split(elt, "_")[end],".")[1:2] ,".") for elt in fft_list_100]) # get stations to iterate over
        println(fft_100_stations[1:5]) # print sample (CHECK)

        pmap(x -> autocorrelate(x, fft_list_100, params), convert.(String, fft_100_stations)) # Run autocorrelations 


        fft_paths_20 = sort(filter(name -> occursin("20_", name), serial_names)) # alphabetic sort for all stations so that diagonal files are made correctly
        fft_paths_1 = sort(filter(name -> occursin("1_", name), serial_names)) # alphabetic sort for all stations so that diagonal files are made correctly
        println("There are $(length(fft_paths_20)) fft datas for the 20 HZ correlations on $dd.") # print 20 Hz summary
        # run 20 HZ correlations
        chunks_20HZ, off_chunk_names_20HZ = get_blocks(fft_paths_20, params); # get chunking to distribute to workers

        # correlate diagonal chunks
        T20D = @elapsed robust_pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_20HZ", true, params), chunks_20HZ)

        # correlate off-diagonal chunks
        T20O = @elapsed robust_pmap(chunk -> offdiag_chunks(chunk, "CORR_20HZ", true, params), off_chunk_names_20HZ) # run mega correlations


        # run 1 Hz correlations 
        chunks_1HZ, off_chunk_names_1HZ = get_blocks(fft_paths_1, params); # chunking

        # correlate diagonal chunks
        T1D = @elapsed robust_pmap(chunk -> diag_chunks(convert(Array,chunk), "CORR_1HZ", false, params), chunks_1HZ)

        # correlate off-diagonal chunks
        T1O = @elapsed robust_pmap(chunk -> offdiag_chunks(chunk, "CORR_1HZ", false, params), off_chunk_names_1HZ) # run mega correlations

        # Correlation summary
        println("All $(length(glob("CORR_*/$(yr)_$(params["month"])/*/*/$path*", params["rootdir"]))) Inter-station Correlations computed in $(T20D + T20O + T1D + T1O) seconds")
        try; rm("$(params["rootdir"])/data/continuous_waveforms/", recursive=true); catch e; println(e); end
        GC.gc()
    catch e 
        println("This is embarrasing... Failed to compute correlations for $dd. Error:", e)
    end
end
dd = Date("2001-03-02")
correlate_big(dd)