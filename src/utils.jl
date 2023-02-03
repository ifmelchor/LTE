#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

using ChaosTools
using DataFrames


"""
    _fqbds(*args)

Return freq time series
"""
function _fqbds(ndata::Int64, fs::Int64, fq_band::Vector{Float64}; pad=1.0)
    
    _, fftleng, halffreq = Multitaper.output_len(range(1,ndata), pad)
    freq    = fs*range(0,1,length=fftleng+1)[1:halffreq]
    freqmin = fq_band[1]
    freqmax = fq_band[2]
    fleft   = findfirst(x -> x >= freqmin, freq)
    fright  = findfirst(freqmax .<= freq)
    #frmin = findmin(abs.(freq.-fq_band[1]))[2]
    #frmax = findmin(abs.(freq.-fq_band[2]))[2]

    return freq[fleft:fright], (fleft, fright)
end


"""
    _PE(*args)

Perform permutation entropy
"""
function _PE(data::Array{Float64,1}, ord::Int64, delta::Int64, fq_band::Vector{Float64}, fs::Int64)

    filt_data = _filt(data, fq_band, fs)
    PE = ChaosTools.permentropy(data, τ=ord, m=delta, base=2)
    normalizedPE = PE / log2(factorial(ord))

    return PE
end


"""
    _filt(*args)

Filter data between fq_band
"""
function _filt(data::Array{Float64,1}, fq_band::Vector{Float64}, fs::Int64; ctr=true, buttorder=4)

    # filter data (with zero phase distortion) of buterworth of order 4

    if ctr
        data .-= mean(data)
    end

    responsetype = Bandpass(fq_band[1], fq_band[2], fs=fs)
    designmethod = Butterworth(buttorder)
    datafilt = filtfilt(digitalfilter(responsetype, designmethod), data)

    return datafilt
end


"""
    _dsar(*args)

Compute DSAR measure from displacement data
"""
function _dsar(data::Array{Float64,1}, fs::Int64; mf_band=(4.,8.), hf_band=(8.,16.), twindow=2, threshold=5.)

    mfdata = abs.(_filt(data, mf_band, fs))
    hfdata = abs.(_filt(data, hf_band, fs))
    dsar = mfdata./hfdata

    # remove outlier
    dsar_f = _removeoutliers(dsar, fs, twindow, threshold)

    return dsar_f
end


"""
    _optparams(*args)

Compute additional LTE parameters from seismic data
"""
function _optparams(s_data::Array{Float64,1}, d_data::Union{Array{Float64,1},Nothing}, fs::Int64; twindow=2, threshold=5.)

    remove_outlier = x -> _removeoutliers(x, fs, twindow, threshold)

    vlf = abs.(_filt(s_data, (.01,.1), fs))
    lf = abs.(_filt(s_data, (.1,2.), fs))
    vlar = vlf./lf
    rsam = abs.(_filt(s_data, (2.,4.5), fs))
    lrar = lf./rsam
    mf = abs.(_filt(s_data, (4.,8.), fs))
    rmar = rsam./mf
    hf = abs.(_filt(s_data, (8.,16.), fs))
    
    fparams = map(remove_outlier, [vlf, lf, vlar, rsam, lrar, mf, rmar, hf])
    vlf_f  = fparams[1]
    lf_f   = fparams[2]
    vlar_f = fparams[3]
    rsam_f = fparams[4]
    lrar_f = fparams[5]
    mf_f   = fparams[6]
    rmar_f = fparams[7]
    hf_f   = fparams[8]

    if !isnothing(d_data)
        dsar_f = _dsar(d_data, fs, twindow=twindow, threshold=threshold)
    else
        dsar_f = nothing
    end
    
    return OptParams(vlf_f, lf_f, vlar_f, rsam_f, lrar_f, mf_f, rmar_f, hf_f, dsar_f)
end

function Base.:+(op1::OptParams, op2::OptParams)
    vlf  = op1.vlf  + op2.vlf
    lf   = op1.lf   + op2.lf
    vlar = op1.vlar + op2.vlar
    rsam = op1.rsam + op2.rsam
    lrar = op1.lrar + op2.lrar
    mf   = op1.mf   + op2.mf
    rmar = op1.rmar + op2.rmar
    hf   = op1.hf   + op2.hf

    if !isnothing(op1.dsar) & !isnothing(op2.dsar)
        dsar = op1.dsar + op2.dsar
    else
        dsar = nothing
    end

    return OptParams(vlf, lf, vlar, rsam, lrar, mf, rmar, hf, dsar)
end

function Base.:/(op::OptParams, n::Int64)
    vlf  = op.vlf  / n
    lf   = op.lf   / n
    vlar = op.vlar / n
    rsam = op.rsam / n
    lrar = op.lrar / n
    mf   = op.mf   / n
    rmar = op.rmar / n
    hf   = op.hf   / n

    if !isnothing(op.dsar)
        dsar = op.dsar / n
    else
        dsar = nothing
    end

    return OptParams(vlf, lf, vlar, rsam, lrar, mf, rmar, hf, dsar)
end


"""
    _standarize(*args)

Standarize a time series
"""
function _standarize(data::Array{Union{Float64, Missing},1})

    u = mean(data)
    s = std(data)
    z = (data .- u) ./ s

    return z
end


"""
    _removeoutliers(*args)

Search for spikes and replace by nan
"""
function _removeoutliers(data::Array{Union{Float64, Missing},1}, fs::Int64, twindow::Int64, threshold::Float64)

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
