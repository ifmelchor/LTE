#!/usr/local/bin julia
# coding=utf-8

# Utility functions for normalization


# Ivan Melchor 2023

using FFTW:rfft, irfft

"""
    _spec_white(*args)

Spectral withening as Benson et al. 2007
"""
function _spec_white(A::AbstractArray{T}, fs::J, win_freq::T, fq_band::Vector{T}) where {T<:Real, J<:Integer}

    npts = size(A,1)
    R = rfft(A)

    deltaf = fs/npts  # frequency step
    halfwindow = convert(Int64, round(win_freq / deltaf / 2.0))

    _moving_avg(R, halfwindow)

    A_w = irfft(R, npts)

    return _filter(A_w, fq_band, fs)
end


"""
    _norm(*args)

Temporal normalization as Benson et al. 2007
"""
function _norm(A::AbstractArray{T}, fs::J, win_freq::T, fq_band::Vector{T}) where {T<:Real, J<:Integer}

    npts = size(A,1)
    R = rfft(A)

    deltaf = fs/npts  # frequency step
    halfwindow = convert(Int64, round(win_freq / deltaf / 2.0))

    _moving_avg(R, halfwindow)

    A_w = irfft(R, npts)

    return _filter(A_w, fq_band, fs)
end


# add demean
# add detrend

    



