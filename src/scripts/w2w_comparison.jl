include("/Users/julianschmitt/Documents/Schoolwork/Seismology/SeisCore.jl/src/SeisCore.jl")
using .SeisCore, SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, JLD2, Statistics, Glob, ColorSchemes, 
            Plots.PlotUtils, HDF5, Images, Statistics, AbstractFFTs, MAT, DSP

# open nan's velocity file and noise file
fmat = matopen("/Volumes/T7/seis_data/W2W/5_SVD_velocity.mat") # total time 240 seconds
locations = DataFrame(CSV.File("/Volumes/T7/seis_data/W2W/fullsocal_inmesh_fixed.csv"));
locations.index = collect(1:919) # add index to match with source_df
fid = HDF5.h5open("/Volumes/T7/seis_data/nodestack/2019/CI.SVD_2019_B1.h5", "r")


# Select Reciever and find position in nan's matlab file
reciever = "B1260"
SVD = filter(x -> occursin(reciever, x.station), locations)
rec_mat = SVD.index[1]
println("$reciever ZZ simulation is located at index: ", rec_mat)

# functions - create CorrData object from mat file and resample CorrData object
function create_matcorr(file::MAT.MAT_v5.Matlabv5File, station::String, df::DataFrame, source::String= "SVD",
        sim_frequency::Float64=500., component::String="V_SzUz")
    """ Returns correlation with metadata from MATLAB file """
    try
    rec_df = filter(x -> occursin(station, x.station), df) 
    mat_index = rec_df.index[1] # get index for matfile
    corr = read(file, component)[:, mat_index]

    # write to corrdata object
    S = CorrData()
    S.corr = reshape(corr,size(corr,1),1)
    S.name = join([rec_df.network[1], rec_df.station[1],"", component[[4,6]]], ".")
    S.fs = sim_frequency # Based on the simulation
    geo = GeoLoc(lat = rec_df.latitude[1], lon = rec_df.longitude[1], el = convert(Float64,rec_df.elevation[1]))
    S.loc = geo

    # get distance - first find source location
    src = filter(x -> occursin(source, x.station), df)
    src_geo = GeoLoc(lat = src.latitude[1], lon = src.longitude[1], el=convert(Float64, src.elevation[1]))
    S.dist = get_dist(src_geo, geo)
    return S
    catch e
        println(e)
        return nothing
    end
end
function c_resample!(C::CorrData,fs::Real)
    @assert fs > 0 "New fs must be greater than 0!"
    if C.fs == fs 
        return nothing 
    end
    T = eltype(C.corr)
    C.corr = resample_kernel(C.corr, T(fs / C.fs))
    C.fs = Float64(fs)
    if C.freqmax >  fs / 2
        C.freqmax = Float64(fs / 2) 
    end
    return nothing
end
c_resample(C::CorrData,fs::Real) = (U = deepcopy(C); c_resample!(U,fs); return U)
function resample_kernel(x::AbstractMatrix,rate::Real)
    T = eltype(x)
    h = resample_filter(rate)
    self = FIRFilter(h, rate)
    # Get delay, in # of samples at the output rate, caused by filtering processes
    τ = timedelay(self)
    # Use setphase! to
    #   a) adjust the input samples to skip over before producing and output (integer part of τ)
    #   b) set the ϕ index of the PFB (fractional part of τ)
    setphase!(self, τ)
    # calculate number of zeros 
    Nrows, Ncols = size(x)
    outLen = ceil(Int, Nrows*rate)
    reqInlen = inputlength(self, outLen)
    reqZerosLen = reqInlen - Nrows
    xPadded     = [x; zeros(eltype(x), reqZerosLen, Ncols)]
    out = zeros(T,outLen,Ncols)
    @inbounds for ii = 1:Ncols
        out[:,ii] .= filt(deepcopy(self), xPadded[:,ii])
    end
    return out
end

function plot_w2w(matfile::MAT.MAT_v5.Matlabv5File, noisefile::HDF5File,
    receiver::String, locations_df::DataFrame, freqs::Array{Float64,1}=[0.2,0.35],
    stacktype::String="linear", component::String="ZZ", source::String="SVD")

    """ Plot waveform comparison between simulation and data. This function 
        assumes that the simulation data is in velocity, while noise data needs 
        a double integral to get to velocity. Flat scaling is done for visibility. 
        Returns L2 norm. """
    noise_corr = get_corrs(noisefile, stacktype, component, receiver)[1];
    sim_corr = create_matcorr(matfile, receiver, locations_df)
    c_resample!(noise_corr, 20.)
    c_resample!(sim_corr, 20.)

    # scale data to match simulation
    sim_data = bandpass(sim_corr,freqs[1], freqs[2]).corr[1:1000]
    noise_data = reverse(bandpass(fftderivative(fftderivative(noise_corr)),freqs[1],freqs[2]).corr[5000:6001])
    noise_data .*= maximum(sim_data)/maximum(noise_data) # scale 

    # plot 
    w2w = plot(sim_data, label="Simulation waveform", color = "blue",
            title="$receiver $component from $source Simulation and $stacktype Noise waveform comparison: $(freqs[1]) to $(freqs[2])")
    plot!(w2w, noise_data, color="red", label = "Noise waveform", dpi = 300)
    fpath = "/Users/julianschmitt/Desktop/SeisPlots2/w2w/$(receiver)_from_$(source)_$(component)_$(stacktype)_$(freqs[1])_to_$(freqs[2]).png"
    if !isdir(dirname(fpath))
        mkpath(dirname(fpath))
    end
    png(w2w, fpath)
end

# iterate through B1 stations
B1_stations = ["B1$(lpad(elt,3,'0'))" for elt in collect(1:260)]
for sta in B1_stations
    try
        plot_w2w(fmat, fid, sta, locations) 
    catch e
        println(sta)
        #println(e)
    end
end