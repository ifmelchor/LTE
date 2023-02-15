#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

using Entropies
using Statistics
using DSP


"""
    _fqbds(*args)

Return freq time series
"""
function _fqbds(ndata::J, fs::J, fq_band::Vector{T}; pad::T=1.0) where {T<:Real, J<:Int}
    
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
function _PE(data::AbstractArray{T}, ord::J, delta::J, fq_band::Vector{T}, fs::J) where {T<:Real, J<:Int}

    filt_data = _filt(data, (fq_band[1],fq_band[2]), fs)
    PE = Entropies.permentropy(data, τ=ord, m=delta, base=2)
    normalizedPE = PE / log2(factorial(ord))

    return PE
end


"""
    _filt(*args)

Filter data between fq_band
"""
function _filt(data::AbstractArray{T}, fq_band::Tuple{T, T}, fs::J; ctr::Bool=true, buttorder::J=4) where {T<:Real, J<:Int}

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
function _dsar(data::AbstractArray{T}, fs::J, twindow::T, threshold::T) where {T<:Real, J<:Int}

    mfdata = abs.(_filt(data, (4.,8.), fs))
    hfdata = abs.(_filt(data, (8.,16.), fs))
    dsar = mfdata./hfdata

    # remove outlier
    dsar_f = _removeoutliers(dsar, convert(Int64, 15*fs), convert(Int64, twindow*fs), threshold)

    return dsar_f
end


"""
    _optparams(*args)

Compute additional LTE parameters from seismic data
"""
function _optparams(s_data::AbstractArray{T}, d_data::Union{Array{T},Nothing}, fs::J, twindow::T, threshold::T) where {T<:Real, J<:Int}

    remove_outlier = x -> _removeoutliers(x, convert(Int64, 15*fs), convert(Int64, twindow*fs), threshold)

    vlf  = abs.(_filt(s_data, (.01,.1), fs))
    lf   = abs.(_filt(s_data, (.1,2.), fs))
    vlar = vlf./lf
    rsam = abs.(_filt(s_data, (2.,4.5), fs))
    lrar = lf./rsam
    mf   = abs.(_filt(s_data, (4.,8.), fs))
    rmar = rsam./mf
    hf   = abs.(_filt(s_data, (8.,16.), fs))
    
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
        dsar_f = _dsar(d_data, fs, twindow, threshold)
    else
        dsar_f = nothing
    end
    
    return OptParams(vlf_f, lf_f, vlar_f, rsam_f, lrar_f, mf_f, rmar_f, hf_f, dsar_f)
end

# function Base.:+(op1::OptParams, op2::OptParams)
#     vlf  = op1.vlf  + op2.vlf
#     lf   = op1.lf   + op2.lf
#     vlar = op1.vlar + op2.vlar
#     rsam = op1.rsam + op2.rsam
#     lrar = op1.lrar + op2.lrar
#     mf   = op1.mf   + op2.mf
#     rmar = op1.rmar + op2.rmar
#     hf   = op1.hf   + op2.hf

#     if !isnothing(op1.dsar) && !isnothing(op2.dsar)
#         dsar = op1.dsar + op2.dsar
#     else
#         dsar = nothing
#     end

#     return OptParams(vlf, lf, vlar, rsam, lrar, mf, rmar, hf, dsar)
# end

# function Base.:/(op::OptParams, n::Int64)
#     vlf  = op.vlf  / n
#     lf   = op.lf   / n
#     vlar = op.vlar / n
#     rsam = op.rsam / n
#     lrar = op.lrar / n
#     mf   = op.mf   / n
#     rmar = op.rmar / n
#     hf   = op.hf   / n

#     if !isnothing(op.dsar)
#         dsar = op.dsar / n
#     else
#         dsar = nothing
#     end

#     return OptParams(vlf, lf, vlar, rsam, lrar, mf, rmar, hf, dsar)
# end


"""
    _standarize(*args)

Standarize a time series
"""
function _standarize(data::AbstractArray{T}) where T<:Real

    u = mean(data)
    s = std(data)
    z = (data .- u) ./ s

    return z
end


"""
    _removeoutliers(*args)

Search for spikes and replace by nan
    compute mean
"""
function _removeoutliers(data::AbstractArray{T}, twindow_in::J, twindow::J, threshold::T; return_mean::Bool=true) where {T<:Real, J<:Int}

    ndat = size(data)[1]
    x = convert(Array{Union{Missing,Float64}}, data)
    z = convert(Array{Union{Missing,Float64}}, _standarize(data))
    cond = z .> threshold

    for i in eachindex(cond)
        if Bool(cond[i])
            nin = i - twindow_in
            nfi = nin + twindow

            if nin < 1
                nin = 1
            end

            if nfi > ndat
                nfi = ndat
            end

            x[nin:nfi] .= missing
        end
    end

    if return_mean
        return mean(skipmissing(x))
    else
        return x
    end
end


