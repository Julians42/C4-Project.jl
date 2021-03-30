include("/Users/julianschmitt/Documents/Schoolwork/Seismology/SeisCore.jl/src/SeisCore.jl")
using .SeisCore, SeisIO, SeisNoise, Plots, Dates, CSV, DataFrames, JLD2, Statistics, Glob, ColorSchemes, 
            Plots.PlotUtils, HDF5, Images, Statistics, AbstractFFTs, MAT, DSP

# open nan's velocity file and noise file
fmat = matopen("/Volumes/T7/seis_data/W2W/5_SVD_velocity.mat") # total time 240 seconds
locations = DataFrame(CSV.File("/Volumes/T7/seis_data/W2W/fullsocal_inmesh_fixed.csv"));
locations.index = collect(1:919) # add index to match with source_df
fid = HDF5.h5open("/Volumes/T7/seis_data/nodestack/2019/CI.SVD_2019_B1.h5", "r")

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
function fftnderivative(A::CorrData, n::Int64)
    """ Computes discrete derivative n times via FFT"""
    C = deepcopy(A)
    for ind in 1:n
        detrend!(C)
        taper!(C)
        F = fft(C.corr,1)
        F .*= AbstractFFTs.fftfreq(length(C.corr)).* 1im .* 2π
        C.corr = real.(ifft(F,1))/(C.fs^2)
    end
    return C
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

function cc_w2w(matfile::MAT.MAT_v5.Matlabv5File, noisefile::HDF5File,
    receiver::String, locations_df::DataFrame, nderivs::Int64=2, freqs::Array{Float64,1}=[0.2,0.35], 
    sim_comp::String="V_SzUz", stacktype::String="linear", component::String="ZZ", source::String="SVD")

    """ Plot waveform comparison between simulation and data. This function 
        assumes that the simulation data is in velocity, while noise data needs 
        a double integral to get to velocity. Flat scaling is done for visibility. 
        Returns L2 norm. """
    noise_corr = get_corrs(noisefile, stacktype, component, receiver)[1];
    sim_corr = create_matcorr(matfile, receiver, locations_df, "SVD", 500., sim_comp)
    c_resample!(noise_corr, 20.)
    c_resample!(sim_corr, 20.)

    # scale data to match simulation
    sim_data = bandpass(sim_corr,freqs[1], freqs[2]).corr[1:1500]
    sim_data ./= maximum(sim_data) # scale to 1
    noise_data = reverse(bandpass(fftnderivative(noise_corr,nderivs),freqs[1],freqs[2]).corr[4502:6001])
    noise_data .*= maximum(sim_data)/maximum(noise_data) # scale 
    
    #Calculate correlation coefficient
    ccoef = sum(noise_data .* sim_data)/(sqrt(sum(sim_data.^2))*(sqrt(sum(noise_data.^2))))
    return [receiver, ccoef]
end
function plot_w2w(matfile::MAT.MAT_v5.Matlabv5File, noisefile::HDF5File,
    receiver::String, locations_df::DataFrame, nderivs::Int64=2, freqs::Array{Float64,1}=[0.2,0.35], 
    sim_comp::String="V_SzUz", stacktype::String="linear", component::String="ZZ", source::String="SVD", 
    line::String="")

    """ Plot waveform comparison between simulation and data. This function 
        assumes that the simulation data is in velocity, while noise data needs 
        a double integral to get to velocity. Flat scaling is done for visibility. 
        Returns L2 norm. """
    noise_corr = get_corrs(noisefile, stacktype, component, receiver)[1];
    sim_corr = create_matcorr(matfile, receiver, locations_df, "SVD", 500., sim_comp)
    c_resample!(noise_corr, 20.)
    c_resample!(sim_corr, 20.)

    # scale data to match simulation
    sim_data = bandpass(sim_corr,freqs[1], freqs[2]).corr[1:1500]
    sim_data ./= maximum(sim_data) # scale to 1
    noise_data = reverse(bandpass(fftnderivative(noise_corr,nderivs),freqs[1],freqs[2]).corr[4502:6001])
    noise_data .*= maximum(sim_data)/maximum(noise_data) # scale 
    
    #Calculate correlation coefficient
    ccoef = sum(noise_data .* sim_data)/(sqrt(sum(sim_data.^2))*(sqrt(sum(noise_data.^2))))
    #L2norm_flipped = 1/2*sum((sim_data .+ noise_data).^2)
    # plot 
    w2w = plot(sim_data, label="Simulation", color = "blue",
            title="$receiver $component from $source Sim and Noise comparison: $(freqs[1]) to $(freqs[2])")
    plot!(w2w, noise_data, color="red", label = "Noise", dpi = 300)
    #plot!(w2w, -noise_data, color="green", label = "Noise Flipped")
    fpath = "/Users/julianschmitt/Desktop/SeisPlots2/w2w_$line/$(sim_comp)/$(receiver)_from_$(source)_$(component)_$(stacktype)_$(freqs[1])_to_$(freqs[2]).png"
    if !isdir(dirname(fpath))
        mkpath(dirname(fpath))
    end
    png(w2w, fpath)
    
    return [receiver, ccoef]
end

# file access across multiple years
fid_2017 = HDF5.h5open("/Volumes/T7/seis_data/nodestack/2017/CI.SVD_2017.h5", "r")
fid_2018 = HDF5.h5open("/Volumes/T7/seis_data/nodestack/2018/CI.SVD_2018.h5", "r")
fid_2019 = HDF5.h5open("/Volumes/T7/seis_data/nodestack/2019/CI.SVD_2019_B1.h5", "r")

# full list of stations
B1_stations = ["B1$(lpad(elt,3,'0'))" for elt in collect(1:260)]
B2_stations = ["B2$(lpad(elt,3,'0'))" for elt in collect(1:60)]
B3_stations = ["B3$(lpad(elt,3,'0'))" for elt in collect(1:86)]
B4_stations = ["B40$(lpad(elt,2,'0'))" for elt in collect(1:96)]
B5_stations = ["B50$(lpad(elt,2,'0'))" for elt in collect(1:61)]
B6_stations = ["B60$(lpad(elt,2,'0'))" for elt in collect(1:33)]
G1_stations = ["G10$(lpad(elt,2,'0'))" for elt in collect(1:60)]
G2_stations = ["G20$(lpad(elt,2,'0'))" for elt in collect(1:50)]
G3_stations = ["G30$(lpad(elt,2,'0'))" for elt in collect(1:33)]
G4_stations = ["G40$(lpad(elt,2,'0'))" for elt in collect(1:14)]

function w2w_wrapper(sta::String)
    line = convert(String, sta[1:2])
    fh5 = nothing
    if any(occursin.(line, ["SB4", "SG2", "SG1"])); fh5 = fid_2017
    elseif any(occursin.(line, ["SB3", "SB2", "SB6"])); fh5 = fid_2018
    else; fh5 = fid_2019
    end
    try
        # write plots
        cc_N = plot_w2w(fmat, fh5, sta, locations, 2, [0.2,0.35], "V_SyUy", "linear", "NN", "SVD", line) 
        cc_Z = plot_w2w(fmat, fh5, sta, locations, 2, [0.2,0.35], "V_SzUz", "linear", "ZZ", "SVD",line) 
        return [cc_N, cc_Z]
    catch e
        println(e)
    end
    println(sta)
end

# map plotting to all nodes
ccs = map(x-> w2w_wrapper(x), vcat(B1_stations, B2_stations, B3_stations, B4_stations, B5_stations,
                                    B6_stations, G1_stations, G2_stations, G3_stations, G4_stations)[collect(1:30:end)])
        

# save correlation coefficients to dataframe 
NN_df = DataFrame([[n[1][1] for n in ccs], [n[1][2] for n in ccs]])
ZZ_df = DataFrame([[n[2][1] for n in ccs], [n[2][2] for n in ccs]])
CSV.write("/Users/julianschmitt/Downloads/ZZ_correlation_coefficients.csv", ZZ_df)
CSV.write("/Users/julianschmitt/Downloads/NN_correlation_coefficients.csv", NN_df)