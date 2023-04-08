#!/usr/local/bin julia
# coding=utf-8

# Utility functions for SAP.jl

# Ivan Melchor


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
  fb2(*args)
    
    Filter buterworth second order desgin
"""
function _fb2(x::Array{T}, fc::T, fs::J, lowpass::Bool; amort=0.47) where {T<:Real, J<:Real}

  a = tan(pi*fc/fs)
  b = 2*a*a - 2
  c = 1 - 2*amort*a + a*a
  d = 1 + 2*amort*a + a*a

  if lowpass
    a0 = a*a/d
    a1 = 2*a0
  else
    a0 = 1/d
    a1 = -2*a0
  end
  
  a2 = a0
  b1 = -b/d
  b2 = -c/d   
  
  ndata = size(x, 1)
  y = Array{T}(undef, ndata)
  y[1] = x[1]
  y[2] = x[2]

  for j in 3:ndata
    y[j] = a0*x[j] + a1*x[j-1] + a2*x[j-2] + b1*y[j-1] + b2*y[j-2]
  end

  return y
end


function _filter!(data::Array{T,2}, fs::J, fq_band::Vector{T}) where {T<:Real, J<:Real}

  nsta = size(data,1)
  for i in 1:nsta
    _filter!(data[i,:], fs, fq_band)
  end

end

function _filter!(data::Array{T,1}, fs::J, fq_band::Vector{T}) where {T<:Real, J<:Real}
    
    fl, fh = fq_band
    temp = _fb2(data, fh, fs, true)
    data = _fb2(temp, fl, fs, false)
    temp = reverse(data)
    data = _fb2(temp, fh, fs, true)
    temp = _fb2(data, fl, fs, false)
    data = reverse(temp)
  
end


function _filter(data::Array{T}, fs::J, fq_band::Vector{T}) where {T<:Real, J<:Real}
    U = deepcopy(data)
    _filter!(U, fs, fq_band)

    return U
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
    x = convert(Array{Union{Missing,Float32}}, data)
    z = convert(Array{Union{Missing,Float32}}, _standarize(data))
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