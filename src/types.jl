#!/usr/local/bin julia
# coding=utf-8


struct LTEBase
    seis_data :: Array{Float64}                    # seismic velocity data
    disp_data :: Union{Array{Float64}, Nothing}    # seismic displacement data for dsar computation
    channels  :: Array{String}                     # list of available channels (ej. "SHZ", "HHN", etc)
    fs        :: Integer                           # sampling rate
    freq      :: Array{Float64}                    # frequency data
    fq_band   :: Tuple{Float64, Float64}           # frequency band
    fqminmax  :: Tuple{Integer, Integer}           # position of fq band
    nfs       :: Integer                           # npts of freq between fq band
    nwin      :: Integer                           # number of windows to reduce
    lwin      :: Integer                           # length of each window
    nswin     :: Union{Integer, Nothing}            # number of subwindows for moving average
    lswin     :: Union{Integer, Nothing}            # length of subwindow for moving average
    nadv      :: Union{Float64, Nothing}            # number of points to advance a subwindow for moving average
    NW        :: Float64                           # number of tapers
    pad       :: Float64                           # padding frequency domain
    add_param :: Bool                              # true/false to compute aditional parameters
    polar     :: Bool                              # true/false to compute polarization parameters
    pe_order  :: Integer                           # permutation entropy order/dimension
    pe_delta  :: Integer                           # permutation entropy tau
    ap_twin   :: Float64                           # time window for aditional parameters
    ap_th     :: Float64                           # threshold for aditional parameters
end


struct Angles
    azimuth :: Float64
    elev    :: Float64
    phyhh   :: Float64
    phyvh   :: Float64
end


struct PParams
    freq    :: Union{Array{Float64}, Nothing}
    degree  :: Array{Float64}
    rect    :: Array{Float64}
    azimuth :: Array{Float64}
    elev    :: Array{Float64}
    phyhh   :: Array{Float64}
    phyvh   :: Array{Float64}
end


struct SParams
    freq      :: Array{Float64}  # frequency series
    S         :: Array{Float64}  # power spectrum density
    erg       :: Float64         # energy
    dominant  :: Float64         # dominant frequency
    centroid  :: Float64         # centroid frequency
end


struct OptParams
    vlf   :: Float64 
    lf    :: Float64
    vlar  :: Float64
    rsam  :: Float64
    lrar  :: Float64
    mf    :: Float64
    rmar  :: Float64
    hf    :: Float64
    dsar  :: Union{Float64, Nothing} 
end