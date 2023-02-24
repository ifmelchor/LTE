#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

using Entropies
using NumericalIntegration
using Statistics
using DSP

"""
   _empty_lte(*args)

Genera un dict vacio para llenar durante el procesado.
"""
function _empty_dict(channels::Vector{String}, add_param::Bool, polar::Bool, base::STABase)
    dict = Dict()
    
    for chan in channels
        dict[chan] = Dict()
        dict[chan]["specgram"] = Array{Float64}(undef, base.nwin, base.nfs)
        for attr in ("energy", "fq_dominant", "fq_centroid", "perm_entr")
            dict[chan][attr] = Array{Float64}(undef, base.nwin)
        end
    end

    if add_param
        dict["opt"] = Dict()
        for attr in ("vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf", "dsar")
            dict["opt"][attr] = Array{Float64}(undef, base.nwin)
        end
    end

    if length(channels)==3 && polar
        dict["polar"] = Dict()
        for attr in ("degree", "rect", "azimuth", "elev", "phyhh", "phyvh")
            dict["polar"][attr] = Array{Float64}(undef, base.nwin, base.nfs)
        end
    end

    return dict
end


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
    _PE(*args)

Perform permutation entropy
"""
function _PE(data::AbstractArray{T}, ord::J, delta::J, fq_band::Vector{T}, fs::J) where {T<:Real, J<:Int}

    filt_data = _filt(data, (fq_band[1],fq_band[2]), fs)
    PE = Entropies.permentropy(filt_data, τ=ord, m=delta, base=2)

    return PE / log2(factorial(ord))
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

Compute DSAR measure
"""
function _dsar(data::AbstractArray{T}, fs::J, twindow::T, threshold::T) where {T<:Real, J<:Int}

    # integrate seismic displacement
    ndata = length(data)
    data   = cumul_integrate(data, range(1, ndata))

    # filter in 4-8 and 8-16 Hz
    mfdata = abs.(_filt(data, (4.,8.), fs))
    hfdata = abs.(_filt(data, (8.,16.), fs))

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

    vlf  = abs.(_filt(data, (.01,.1), fs))
    lf   = abs.(_filt(data, (.1,2.), fs))
    vlar = vlf./lf
    rsam = abs.(_filt(data, (2.,4.5), fs))
    lrar = lf./rsam
    mf   = abs.(_filt(data, (4.,8.), fs))
    rmar = rsam./mf
    hf   = abs.(_filt(data, (8.,16.), fs))

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


