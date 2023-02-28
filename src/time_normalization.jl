#!/usr/local/bin julia
# coding=utf-8

# Utility functions for SAP.jl

# These functions were taken from SeisNoise julia package
# see more in github.com/tclements/SeisNoise.jl

# Ivan Melchor 2023

using DSP, FFTW


"""
    onebit!(A)

One-bit amplitude modification of Data A.
"""
function onebit!(A::AbstractArray)
  A .= sign.(A)
  return nothing
end

onebit(A::AbstractArray) = (U = deepcopy(A); onebit!(U);return U)



"""
  mute(A,factor)

Set high amplitudes in array `A` to zero.
Uses median of envelope of `A` to find outliers.
"""
function mute!(A::AbstractArray,factor::Real=3)
    T = eltype(A)
    envelope = abs.(hilbert(A))
    levels = mean(envelope,dims=1)
    level = factor .* median(levels)
    A[envelope .> level] .= T(0)
    return nothing
end
mute(A::AbstractArray,factor::Real=3) = (U = deepcopy(A); mute!(U,factor);
     return U)


"""
  clip(A,factor)
Truncate array A at `factor` times the root mean square of each column.
#Arguments
- `A::AbstractArray`: N-d time series array
- `factor::Real`:
- `f::Function`: Input statistical function (e.g. rms, var, std, mad)
- `dims`: Dimension of `A` to apply clipping (defaults to 1)
"""
function clip!(A::AbstractArray{T,N}, factor::Real; f::Function=std, dims=1) where {T,N}
    if N == 1
        high = f(A) .* factor
        clamp!(@view(A[:]), -high, high)
    else
        high = f(@view(A[:,:]), dims=dims) .* factor
        for ii = 1:size(A,2)
            clamp!(@view(A[:,ii]), -high[ii], high[ii])
        end
    end
    return nothing
end
clip(A::AbstractArray, factor::Real; f::Function=std, dims=1) = (U = deepcopy(A);
     clip!(U,factor,f=f,dims=dims);return U)


"""
    smooth(A, half_win)
Smooth array `A` with half-window `half_win` (defaults to 3).
"""
function smooth!(A::AbstractArray, half_win::Int=3, dims::Int=1)
    T = eltype(A)
    window_len = 2 * half_win + 1
    csumsize = tuple(collect(size(A)) .+ [i==1 ? 2 * (window_len - 1) + 1 : 0 for i in 1:ndims(A)]...)
    csum = similar(A,T,csumsize)
    csum[1:window_len,:] .= zero(T)
    csum[end - window_len + 1:end,:] .= zero(T)
    csum[window_len+1:end-window_len + 1,:] .= A
    csum .= cumsum(csum,dims=dims)
    weight = similar(A,T,size(A,1))
    weight[1:half_win] = T.(window_len ÷ 2 + 1:window_len - 1)
    weight[half_win + 1: end - half_win] .= T(window_len)
    weight[end-half_win:end] = T.(window_len:-1:window_len ÷ 2 + 1)
    A[:,:] .= (csum[window_len+half_win+1:end-half_win,:] .- csum[half_win+1:end-window_len-half_win,:]) ./ weight
    return nothing
end
smooth(A::AbstractArray,half_win::Int=3, dims::Int=1) =
      (U = deepcopy(A);smooth!(U,half_win,dims);return U)

    



