#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

using ChaosTools
using Statistics
using DataFrames
using DSP


function _PE(data::Array{Float64,1}, ord::Integer, delta::Integer, fq_band::Tuple{Float64,Float64}, fs::Integer)

    filt_data = _filt(data, fq_band, fs)
    PE = ChaosTools.permentropy(data, τ=ord, m=delta, base=2)
    normalizedPE = PE / log2(factorial(ord))

    return PE
end


function _filt(data::Array{Float64,1}, fq_band::Tuple{Float64,Float64}, fs::Integer; ctr=true, buttorder=4)

    # filter data (with zero phase distortion) of buterworth of order 4

    if ctr
        data .-= mean(data)
    end

    responsetype = Bandpass(fq_band[1], fq_band[2], fs=fs)
    designmethod = Butterworth(buttorder)
    datafilt = filtfilt(digitalfilter(responsetype, designmethod), data)

    return datafilt
end


function _dsar(data::Array{Float64,1}, fs::Integer; mf_band=(4.,8.), hf_band=(8.,16.), twindow=2, threshold=5.)

    mfdata = abs.(_filt(data, mf_band, fs))
    hfdata = abs.(_filt(data, hf_band, fs))
    dsar = mfdata./hfdata

    # remove outlier
    dsar_f = _removeoutliers(dsar, fs, twindow, threshold)

    return dsar_f
end


function _optparams(data::Array{Float64,1}, fs::Integer; twindow=2, threshold=5.)

    remove_outlier = x -> _removeoutliers(x, fs, twindow, threshold)

    vlf = abs.(_filt(data, (.01,.1), fs))
    lf = abs.(_filt(data, (.1,2.), fs))
    vlar = vlf./lf
    rsam = abs.(_filt(data, (2.,4.5), fs))
    lrar = lf./rsam
    mf = abs.(_filt(data, (4.,8.), fs))
    rmar = rsam./mf
    hf = abs.(_filt(data, (8.,16.), fs))
    
    fparams = map(remove_outlier, [vlf, lf, vlar, rsam, lrar, mf, rmar, hf])
    vlf_f  = fparams[1]
    lf_f   = fparams[2]
    vlar_f = fparams[3]
    rsam_f = fparams[4]
    lrar_f = fparams[5]
    mf_f   = fparams[6]
    rmar_f = fparams[7]
    hf_f = fparams[8]

    return OptParams(vlf_f, lf_f, vlar_f, rsam_f, lrar_f, mf_f, rmar_f, hf_f)
end


function _standarize(data::Array{Union{Float64, Missing},1})

    u = mean(data)
    s = std(data)
    z = (data .- u) ./ s

    return z
end


function _removeoutliers(data::Array{Union{Float64, Missing},1}, fs::Integer, twindow::Integer, threshold::Float64)

    ndat = size(data)[1]
    z = _standarize(data)
    cond = z .> threshold

    for i in eachindex(cond)
        if convert(Bool,cond[i])
            nin = i - twindow*fs
            nfi = i + twindow*fs

            if nin <= 0
                nin = 1
            end

            if nfi >= ndat
                nfi = ndat
            end

            z[nin:nfi] = missing
        end
    end

    return mean(skipmissing(z))
end


# for a future
#function _removeresponse(data)
