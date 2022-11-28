#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

using Multitaper


function _fqbds(ndata::Integer, fs::Integer, fq_band::Tuple{Float64,Float64}, pad=1.0)
    lengt, fftleng, halffreq = Multitaper.output_len(range(ndata), pad)
    freq = fs*range(0,1,length=fftleng+1)[1:halffreq]
    frmin = findmin(abs.(freq.-fq_band[1]))[2]
    frmax = findmin(abs.(freq.-fq_band[2]))[2]
    freq = freq[frmin:frmax]

    return freq, (frmin, frmax)
end


function _spectral(data::Array{Float64,1}, fs::Integer, fqr::Union{Tuple{Integer,Integer}, Nothing}, fq_band::Union{Tuple{Float64,Float64}}; NW=3.5, pad=1.0)

    # multitaper psd of data_i
    si = multispec(data, ctr=true, dt=1/fs, NW=NW, K=K, pad=pad)

    # limit to fq bounds
    if !ismissing(fqr)
        frmin = findmin(abs.(si.freq.-fq_band[1]))[2]
        frmax = findmin(abs.(si.freq.-fq_band[2]))[2]
    else
        frmin = fqr[1]
        frmax = fqr[2]
    end
    
    # get true freq and psd
    freq = si.freq[frmin:frmax]
    S = si.S[frmin:frmax]
    
    # dominant freq
    domfreq = freq[findmax(S)[2]]

    # centroid freq
    cenfreq = dot(freq,S)/sum(S)

    # energy
    erg = sum(S)

    return SParams(freq, S, erg, domfreq, cenfreq)
end


function _crosscorr(data_i::Array{Float64,1}, data_j::Array{Float64,1}, fs::Integer; NW=3.5, pad=1.0)
    
    K = convert(Int64, 2*NW - 1)
    sij = multispec(data_i, data_j, outp=:spec, dt=1/fs, NW=NW, K=K, ctr=true, pad=pad, guts=true)
    Si = sum(sij.coef[1].coef, dims=2)
    Sj = sum(sij.coef[2].coef, dims=2)

    return Si .* conj(Sj)
end


function _csm(data::Array{Float64,2}, fs::Integer, fqr::Union{Tuple{Integer,Integer}, Nothing}, fq_band::Tuple{Float64,Float64}; NW=3.5, pad=1.0)

    if !ismissing(fqr)
        ndata = size(data[1,:])
        freq, (frmin, frmax) = _fqbds(ndata, fs, fq_band, pad)
        nfs = size(freq)[1]
    else
        frmin = fqr[1]
        frmax = fqr[2]
        nfs = frmax-frmin+1
        freq = nothing
    end

    csm = Array{ComplexF64}(undef, 3, 3, nfs)
    for i in 1:3
        for j in i:3
            csm[i,j,:] = _crosscorr(data[i,:], data[j,:], fs, NW=NW, pad=pad)[frmin:frmax, 1]
            if i != j
                csm[j,i,:] = conj(csm[i,j,:])
            end
        end
    end

    # do singular value decomposition
    csm_svd = map(svd,[csm[:,:,i] for i in 1:size(csm)[3]])

    return freq, csm_svd
end