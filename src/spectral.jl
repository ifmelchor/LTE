#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

"""
    _psd(args)

Compute the psd using multitaper approach
"""
function _psd(data::AbstractArray{T}, fs::J, lwin::Union{J,Nothing}, nwin::Union{J,Nothing}, nadv::Union{T,Nothing}, fqr::Tuple{J,J}, NW::T, pad::T) where {T<:Real, J<:Int}

    npts = size(data)
    frmin, frmax = fqr
    nfs  = frmax-frmin+1
    
    psd  = zeros(Float32, nfs)
    K    = convert(Int32, 2*NW - 1)

    if !isnothing(nwin) && !isnothing(lwin)

        if isnothing(nadv)
            nadv = 0
        end
        
        for n in 1:nwin
            n0  = 1 + floor(Int32, lwin * nadv*(n-1))
            nf  = floor(Int32, n0 + lwin)
            s_n = multispec(data[n0:nf], ctr=true, dt=1/fs, NW=NW, K=K, pad=pad)
            psd .+= s_n.S[frmin:frmax]
            
            if n == 1
                freq = s_n.f[frmin:frmax]
            end
        end
        
        psd ./= nwin
    else

        s    = multispec(data, ctr=true, dt=1/fs, NW=NW, K=K, pad=pad)
        psd  = s.S[frmin:frmax]
        freq = s.f[frmin:frmax]
    end 

    return psd, freq
end


"""
    _crosscorr(args)

Compute the cross spectral correlation of two signals
"""
function _crosscorr(data_i::AbstractArray{T}, data_j::AbstractArray{T}, fs::J, NW::T, pad::T) where {T<:Real, J<:Int}
    
    K   = convert(Int32, 2*NW - 1)
    sij = multispec(data_i, data_j, outp=:spec, dt=1/fs, NW=NW, K=K, ctr=true, pad=pad, guts=true)
    Si  = sum(sij.coef[1].coef, dims=2)
    Sj  = sum(sij.coef[2].coef, dims=2)

    return Si .* conj(Sj)
end


"""
    _csm(args)

Compute the SVD of the average moving Hermitian covariance spectral matrix (CSM) of data.
This is useful for spectral width measurement and polarization analysis
"""
function _csm(data::AbstractArray{T}, fs::J, lwin::Union{J,Nothing}, nwin::Union{J,Nothing}, nadv::Union{T,Nothing}, fqr::Tuple{J,J}, NW::T, pad::T) where {T<:Real, J<:Int}

    # define npts and nro of components 
    # (for polarization analysis ncomp is 3)
    ncomp, npts = size(data)
    frmin, frmax = fqr
    nfs = frmax - frmin + 1
    covm = zeros(ComplexF64, ncomp, ncomp, nfs)
    
    # build half part of covm matrix
    if !isnothing(nwin) && !isnothing(lwin)

        if isnothing(nadv)
            nadv = 0.0
        end
        
        for n in 1:nwin
            n0  = 1 + floor(Int, lwin * nadv*(n-1))
            nf  = floor(Int, n0 + lwin)
            for i in 1:ncomp
                for j in i:ncomp
                    covm[i,j,:] .+= _crosscorr(data[i,n0:nf], data[j,n0:nf], fs, NW, pad)[frmin:frmax, 1]
                end
            end
        end

        covm ./= nwin
    else
        for i in 1:ncomp
            for j in i:ncomp
                covm[i,j,:] = _crosscorr(data[i,:], data[j,:], fs, NW, pad)[frmin:frmax,1]
            end
        end
    end
    
    # build full covm matrix
    for i in 1:ncomp
        for j in i:ncomp
            if i != j
                covm[j,i,:] = conj(covm[i,j,:])
            end
        end
    end

    # do singular value decomposition
    covm_svd = map(svd,[covm[:,:,i] for i in 1:nfs])

    return covm_svd
end