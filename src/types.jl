#!/usr/local/bin julia
# coding=utf-8


struct LTEBase{T<:Real, J<:Int}
    seis_data :: Array{T}                    # seismic velocity data
    disp_data :: Union{Array{T}, Nothing}    # seismic displacement data for dsar computation
    channels  :: Vector{String}              # list of available channels (ej. "SHZ", "HHN", etc)
    fs        :: J                           # sampling rate
    freq      :: Array{T}                    # frequency data
    fq_band   :: Vector{T}                   # frequency band
    fqminmax  :: Tuple{J, J}                 # position of fq band
    nfs       :: J                           # npts of freq between fq band
    nwin      :: J                           # number of windows to reduce
    lwin      :: J                           # length of each window
    nswin     :: Union{J, Nothing}           # number of subwindows for moving average
    lswin     :: Union{J, Nothing}           # length of subwindow for moving average
    nadv      :: Union{T, Nothing}           # number of points to advance a subwindow for moving average
    NW        :: T                           # number of tapers
    pad       :: T                           # padding frequency domain
    add_param :: Bool                        # true/false to compute aditional parameters
    polar     :: Bool                        # true/false to compute polarization parameters
    pe_order  :: J                           # permutation entropy order/dimension
    pe_delta  :: J                           # permutation entropy tau
    ap_twin   :: T                           # time window for aditional parameters
    ap_th     :: T                           # threshold for aditional parameters
end


struct Angles{T<:Real}
    azimuth :: T
    elev    :: T
    phyhh   :: T
    phyvh   :: T
end


struct PParams{T<:Real}
    degree  :: Array{T}
    rect    :: Array{T}
    azimuth :: Array{T}
    elev    :: Array{T}
    phyhh   :: Array{T}
    phyvh   :: Array{T}
end


struct SParams{T<:Real}
    freq      :: Array{T}                    # frequency series
    S         :: Array{T}                    # power spectrum density
    erg       :: T                           # energy
    dominant  :: T                           # dominant frequency
    centroid  :: T                           # centroid frequency
end


struct OptParams{T<:Real}
    vlf   :: T 
    lf    :: T
    vlar  :: T
    rsam  :: T
    lrar  :: T
    mf    :: T
    rmar  :: T
    hf    :: T
    dsar  :: Union{T, Nothing} 
end