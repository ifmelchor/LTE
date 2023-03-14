#!/usr/local/bin julia
# coding=utf-8

# Utility functions for SAP.jl

# These functions were taken from SeisNoise julia package
# see more in github.com/tclements/SeisNoise.jl

# Ivan Melchor 2023


using NumericalIntegration
using Entropies


"""
    _PE(*args)

Perform permutation entropy
"""
function _PE(data::AbstractArray{T}, ord::J, delta::J, fq_band::Vector{T}, fs::J) where {T<:Real, J<:Int}

    filt_data = _filter(data, fq_band, fs)
    PE = Entropies.permentropy(filt_data, τ=ord, m=delta, base=2)

    return PE / log2(factorial(ord))
end


"""
    _dsar(*args)

Compute DSAR measure
"""
function _dsar(data::AbstractArray{T}, fs::J, twindow::T, threshold::T) where {T<:Real, J<:Int}

    # integrate seismic displacement
    ndata = length(data)
    data   = cumul_integrate(data, range(1, ndata))

    # filter in 4-8 and 8-16 Hz
    mfdata = abs.(_filter(data, [4.,8.], fs))
    hfdata = abs.(_filter(data, [8.,16.], fs))

    # do the ratio
    dsar   = mfdata./hfdata

    # remove outlier
    if threshold > 0
        dsar_f = _removeoutliers(dsar, convert(Int64, 15*fs), convert(Int64, twindow*fs), threshold)
    else
        dsar_f = mean(dsar)
    end
        
    return dsar_f
end


"""
    _optparams(*args)

Compute additional LTE parameters from seismic data
"""
function _optparams(data::AbstractArray{T}, fs::J, twindow::T, threshold::T) where {T<:Real, J<:Int}

    remove_outlier = x -> _removeoutliers(x, convert(Int64, 15*fs), convert(Int64, twindow*fs), threshold)

    vlf  = abs.(_filter(data, [.01,.1], fs))
    lf   = abs.(_filter(data, [.1,2.], fs))
    vlar = vlf./lf
    rsam = abs.(_filter(data, [2.,4.5], fs))
    lrar = lf./rsam
    mf   = abs.(_filter(data, [4.,8.], fs))
    rmar = rsam./mf
    hf   = abs.(_filter(data, [8.,16.], fs))

    if threshold > 0
        fparams = map(remove_outlier, [vlf, lf, vlar, rsam, lrar, mf, rmar, hf])
        vlf_f  = fparams[1]
        lf_f   = fparams[2]
        vlar_f = fparams[3]
        rsam_f = fparams[4]
        lrar_f = fparams[5]
        mf_f   = fparams[6]
        rmar_f = fparams[7]
        hf_f   = fparams[8]
    else
        vlf_f  = mean(vlf)
        lf_f   = mean(lf)
        vlar_f = mean(vlar)
        rsam_f = mean(rsam)
        lrar_f = mean(lrar)
        mf_f   = mean(mf)
        rmar_f = mean(rmar)
        hf_f   = mean(hf)
    end
    
    dsar_f = _dsar(data, fs, twindow, threshold)
    
    return OptParams(vlf_f, lf_f, vlar_f, rsam_f, lrar_f, mf_f, rmar_f, hf_f, dsar_f)
end

