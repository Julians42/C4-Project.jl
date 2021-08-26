# Helper Routines for basin_correlations.jl

export get_dict_name, get_scedc_files, LLE_geo, add_location, rfft_raw # preprocessing routines
export cc_medianmute, cc_medianmute!, cc_medianmute2!, remove_medianmute # medianmute 
export save_named_corr, write_jld2, load_corrs, sum_corrs, name_corr, foldersize # post-processing routines


#################### Get file names from SCEDC (uses Tim's SCEDC.jl Pkg) #############
function get_dict_name(file::String)
    """ Helper function for get_scedc_files."""
    station = convert(String, split(split(file,"/")[end],"_")[1])
    component = split(file,"/")[end][10:12]
    return string(station, "_", component)
end
function get_scedc_files(dd::Date, params::Dict=params)
    aws, network = params["aws"], params["network"]
    ar_filelist = pmap(x -> s3query(aws, dd, enddate = dd, network=network, channel=x),["BH?", "HH?"])
    filelist_scedc_BH = ar_filelist[1]
    filelist_scedc_HH = ar_filelist[2]
    # create dictionary and overwrite HH keys with available BH data

    BH_keys = [get_dict_name(file) for file in filelist_scedc_BH]
    HH_keys = [get_dict_name(file) for file in filelist_scedc_HH]

    # Convert to dictionary 
    HH_dict = Dict([(name, file) for (name, file) in zip(HH_keys, filelist_scedc_HH)]) 
    BH_dict = Dict([(name, file) for (name, file) in zip(BH_keys, filelist_scedc_BH)]) 
    filelist_dict = merge(HH_dict, BH_dict) # BH dict overwrite HH_dict. This is essentually the union
    filelist_scedc = collect(values(filelist_dict)) # return values as array for download
    return filelist_scedc
end

################ Adds locations to SCEDC files ##################
function LLE_geo(station, df)
    """ Find station matching location and return geoloc object"""
    try
        row = df[(findfirst(x -> x==station, df.station)),:]
        lat, lon, el = row.latitude[1], row.longitude[1], row.elevation[1]
        geo = GeoLoc(lat = float(lat), lon = float(lon), el = float(el))
        return geo
    catch 
        return nothing
    end
end
function add_location(s::SeisData,df::DataFrame)
    """ Adds locations to array of seisdata from a dataframe """
    name = split(s.name[1],".")[2]
    geo = LLE_geo(name, df)
    if !isnothing(geo)
        s.loc[1] = geo
    else
        println("Station $name can't be found in the dataframe")
    end
end


############ raw waveform integral via FFT ##############
function rfft_raw(R::RawData,dims::Int=1)
    FFT = rfft(R.x,dims)
    FFT ./= rfftfreq(size(R.x,1), R.fs) .* 1im .* 2π # Integrate the accelerometers!
    FFT[1,:] .=0
    return FFTData(R.name, R.id,R.loc, R.fs, R.gain, R.freqmin, R.freqmax,
                    R.cc_len, R.cc_step, R.whitened, R.time_norm, R.resp,
                    R.misc, R.notes, R.t, FFT)
end

################### Signal Medianmute ###################
function cc_medianmute(A::AbstractArray, cc_medianmute_α::Float64 = 10.0)
    """
        Remove noisy correlation windows before stacking
        - Remove if average noise is greater than 10x the average
    """
    T, N = size(A)
    cc_maxamp = vec(maximum(abs.(A), dims=1))
    cc_medianmax = median(cc_maxamp)
    inds = findall(x-> x <= cc_medianmute_α*cc_medianmax,cc_maxamp)
    return A[:, inds], inds
end
remove_medianmute(C::CorrData, inds) = (return C.t[inds])
function cc_medianmute!(C::CorrData, cc_medianmute_α::Float64 = 10.0)
    C.corr, inds = cc_medianmute(C.corr, cc_medianmute_α)
    C.t = remove_medianmute(C, inds)
    return nothing
end
function cc_medianmute2!(C::CorrData, cc_medianmute_α::Float64 = 10.0, bool::Bool = true)
    """ Median mute function which ignore metadata when bool=false"""
    C.corr, inds = cc_medianmute(C.corr, cc_medianmute_α)
    if bool
        C.t = remove_medianmute(C, inds)
    end
    return nothing
end


############# Save single day correlations #################
function save_named_corr(C::CorrData, CORROUT::String)
    """ Implements custom naming scheme for project """
    CORROUT = expanduser(CORROUT) # ensure file directory exists
    if isdir(CORROUT) == false
        mkpath(CORROUT)
    end

    yr,j_day = Dates.year(Date(C.id)), lpad(Dates.dayofyear(Date(C.id)),3,"0") # YEAR, JULIAN_DAY
    p_name = name_corr(C) 
    name = join([yr,j_day,p_name],"_") #YEAR_JDY_CORRNAME

    # create JLD2 file and save correlation
    filename = joinpath(CORROUT,"$(name).jld2")
    file = jldopen(filename, "a+")
    if !(C.comp in keys(file))
        group = JLD2.Group(file, C.comp)
        group[C.id] = C
    else
        file[C.comp][C.id] = C
    end
    close(file)
end

function write_jld2(corr_folder::String, file_dir::String)
    pair_corr_names = glob("$corr_folder/*/*.jld2","CORR")
    pair_corrs = load_corrs(pair_corr_names)
    if !isdir("$file_dir") # ensure filepathing
        mkpath("$file_dir")
    end
    fo = jldopen("$(file_dir)/$(corr_folder).jld2", "w") # name of the station pairs, most computationally expensive
    for (index, value) in enumerate(pair_corrs)
        # continue if no data in corr 
        isempty(value.corr) && continue
        comp = value.comp
        starttime = string(Date(u2d(value.t[1]))) # eg 2017-01-01
        groupname = joinpath(comp, starttime)
        !haskey(fo, comp) && JLD2.Group(fo, comp) # component layer
        !haskey(fo, groupname) && (fo[groupname] = value) # save corr into layer 
    end
    close(fo)
end



############# Load list of correlations #####################
function load_corrs(file_list::Array{String,1})
    """ Load list of correlations, extract component"""
    corrs_in_pair = Array{CorrData,1}(undef, length(file_list))
    for (index, name) in enumerate(file_list)
        comp = string(split(name, "/")[end-1]) # get component
        corrs_in_pair[index] = load_corr(name, comp)
    end
    return corrs_in_pair
end


########### Turns array of SeisData's into single SeisData object ##########
function sum_corrs(corrs::Array{CorrData,1})
    """ Implement sum function on array of corrs before stacking"""
    try
        sum_corr = Array{Float64, 2}(undef, size(corrs[1].corr)[1], length(corrs))
        for (ind, corr) in enumerate(corrs)
            sum_corr[:, ind] = corr.corr[:]
        end
        corr = deepcopy(corrs[1]) # keep metadata from first correlation
        corr.corr = sum_corr # update array of correlations
        return corr
    catch e 
        println("Error combining correlations: ", c)
        return nothing
    end
end  


###################### Miscelaneous helpers ####################
function name_corr(C::CorrData)
    """ Returns corr name string: CH1.STA1.LOC1.CH2.STA2.LOC2 """
    return strip(join(deleteat!(split(C.name,"."),[4,8]),"."),'.')
end

function foldersize(dir=".")
    """ returns total size of folder in GB """
    size = 0
    for (root, dirs, files) in walkdir(dir)
        size += sum(map(filesize, joinpath.(root, files)))
    end
    return size*10e-10
end
