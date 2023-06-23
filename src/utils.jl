#!/usr/local/bin julia
# coding=utf-8

# Utility functions for SAP.jl

# Ivan Melchor

using DSP
using PointwiseKDEs

"""
    _fqbds(*args)

Return freq time series
"""
function _fqbds(ndata::J, fs::J, fq_band::Vector{T}; pad::T=1.0) where {T<:Real, J<:Integer}
    
    _, fftleng, halffreq = Multitaper.output_len(range(1,ndata), pad)
    freq    = Array{Float32}(fs*range(0,1,length=fftleng+1)[1:halffreq])
    freqmin = fq_band[1]
    freqmax = fq_band[2]
    fleft   = findfirst(x -> x >= freqmin, freq)
    fright  = findfirst(freqmax .<= freq)
    #frmin = findmin(abs.(freq.-fq_band[1]))[2]
    #frmax = findmin(abs.(freq.-fq_band[2]))[2]

    return freq[fleft:fright], (fleft, fright)
end


"""
    _filter(*args)

Filter a time series using a Bandpass Butterworth of second order
"""
function _filter(data::Array{T}, fs::J, fq_band::Vector{T}) where {T<:Real, J<:Real}
    # U = deepcopy(data)
    # _filter!(U, fs, fq_band)

    filter = digitalfilter(Bandpass(fq_band[1], fq_band[2], fs=fs), Butterworth(2))

    return filtfilt(filter, data)
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
    _scaler(*args)

MinMax scaler of a time series
"""
function _scaler(x::AbstractArray{T}; smin=0., smax=1.) where T<:Real
    x_std = (x .- minimum(x)) ./ (maximum(x) - minimum(x)) 
    x_scaled = x_std .* (smax - smin) .+ smin
    return x_scaled
end



"""
    _removeoutliers(*args)

Search for spikes and replace by nan
    compute mean
"""
function _removeoutliers(data::AbstractArray{T}, twindow_in::J, twindow::J, threshold::T; return_mean::Bool=true) where {T<:Real, J<:Integer}

    ndat = size(data)[1]
    x = convert(Array{Float32}, data)
    z = convert(Array{Float32}, _standarize(data))
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

            x[nin:nfi] .= NaN
        end
    end

    # exclude NaN
    x = x[map(!isnan,x)]

    if length(x) > 1
        if return_mean
            return mean(x)
        else
            return x
        end
    else
        return NaN
    end
end


"""
    _moving_avg(*args)

normalization moving average as Benson et al. 2007
"""
function _moving_avg!(data::AbstractArray{T}, halfwindow::J) where {T<:Any, J<:Integer}

    # pad data
    npts  = size(data, 1)
    z = zeros(Float32, halfwindow)
    data_pad = vcat(z, abs.(data), z)

    # get total size of the averaging window
    npt = 2 * halfwindow

    # do moving average
    for i in 1:npts
        ak = data_pad[i:i+npt] 
        w = sum(ak)/sum(ak.>0)
        data[i] /= w
    end

end


"""
    _circular_avg(*args)

circular mean
"""
function _circmean(data::Array{T}; high=pi, low=0.) where T<:Real

    data_to_rad = data .* pi/180
    sin_data = sin.((data_to_rad .- low)*2*pi ./ (high - low))
    cos_data = cos.((data_to_rad .- low)*2*pi ./ (high - low))
    res = atan(sum(sin_data), sum(cos_data))
    res += 2*pi
    circmean = res*(high-low)/2.0/pi + low
    circmean *= 180/pi

    if circmean > 180
        circmean -= 180
    end

    return circmean
end


"""
    _circular_std(*args)

circular standar deviation
"""
function _circstd(data::Array{T}; high=pi, low=0., normalize=false) where T<:Real

    data_to_rad = data .* pi/180
    sin_data = mean( sin.((data_to_rad .- low)*2*pi ./ (high - low)) )
    cos_data = mean( cos.((data_to_rad .- low)*2*pi ./ (high - low)) )
    R = minimum([1, hypot(sin_data, cos_data)])
    res = sqrt(-2*log(R))

    if !normalize
        res *= (high-low)/(2*pi)
    end

    return res
end



function _mostprobval(data::Array{T}) where T<:Real
    data = reshape(data, (1, :))
    data = convert(Array{Float64}, data)
    data_min = findmin(data)[1]
    data_max = findmax(data)[1]
    x_space = LinRange(data_min, data_max, 100)

    # compute the probability density function
    kde = PointwiseKDE(data)
    y_space = reshape(rand(kde, 100),100)

    return x_space[findmax(y_space)[2]]

end