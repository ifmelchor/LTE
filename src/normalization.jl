#!/usr/local/bin julia
# coding=utf-8

# Utility functions for normalization


# Ivan Melchor 2023


"""
    _spec_white(*args)

Spectral withening as Benson et al. 2007
"""
function _spec_white!(A::AbstractArray{T,1}, fs::J, win_freq::T, fq_band::Vector{T}) where {T<:Real, J<:Integer}

    npts = size(A,1)
    R = rfft(A)

    deltaf = fs/npts  # frequency step
    halfwindow = convert(Int32, round(win_freq / deltaf / 2.0))

    _moving_avg!(R, halfwindow)

    A = irfft(R, npts)
    _filter!(A, fs, fq_band)

end

function _spec_white!(A::AbstractArray{T,2}, fs::J, win_freq::T, fq_band::Vector{T}) where {T<:Real, J<:Integer}
    nsta = size(A,1)
    for n in 1:nsta
        _spec_white!(A[n,:], fs, win_freq, fq_band)
    end

end

function _spec_white(data::AbstractArray{T}, fs::J, win_freq::T, fq_band::Vector{T}) where {T<:Real, J<:Integer}
    U = deepcopy(data)
    _spec_white!(U, fs, win_freq, fq_band)

    return U
end


"""
    _norm(*args)

Temporal normalization as Benson et al. 2007
"""
function _norm(A::AbstractArray{T}, fs::J, win_freq::T, fq_band::Vector{T}) where {T<:Real, J<:Integer}

    npts = size(A,1)
    R = rfft(A)

    deltaf = fs/npts  # frequency step
    halfwindow = convert(Int32, round(win_freq / deltaf / 2.0))

    _moving_avg(R, halfwindow)

    A_w = irfft(R, npts)

    return _filter(A_w, fq_band, fs)
end


# add demean
# add detrend

    



