#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

using DSP, FFTW


"""
    _fqbds(*args)

Return freq time series
"""
function _fqbds(ndata::J, fs::J, fq_band::Vector{T}; pad::T=1.0) where {T<:Real, J<:Int}
    
    _, fftleng, halffreq = Multitaper.output_len(range(1,ndata), pad)
    freq    = Array{Float64}(fs*range(0,1,length=fftleng+1)[1:halffreq])
    freqmin = fq_band[1]
    freqmax = fq_band[2]
    fleft   = findfirst(x -> x >= freqmin, freq)
    fright  = findfirst(freqmax .<= freq)
    #frmin = findmin(abs.(freq.-fq_band[1]))[2]
    #frmax = findmin(abs.(freq.-fq_band[2]))[2]

    return freq[fleft:fright], (fleft, fright)
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
    _spec_white(*args)

Spectral withening as Benson et al. 2007
"""
function _spec_white(A::AbstractArray{T}, fs::J, win_freq::T, fq_band::Tuple{T,T}) where {T<:Real, J<:Integer}

    npts = size(A,1)
    R = rfft(A)

    rnpts = size(R,1)
    R_w = Array{T}(undef, rnpts)

    deltaf = fs / npts  # frequency step

    # smoothing amplitude spectrum
    halfwindow = convert(Int64, round(win_freq / deltaf / 2.0))

    z = zeros(T, halfwindow)
    Apadded =  vcat(z, abs.(R), z)
    
    # total size of the averaging window
    npt = 2 * halfwindow + 1

    # do moving average
    for i in 1:npts
        ak = a_padded[i:i+npt] 
        w = sum(ak)/sum(ak>0)
        R_w[i] = R / w
    end

    A_w = irfft(R_w, npts)

    return _filt(A_w, fq_band, fs)
end


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


