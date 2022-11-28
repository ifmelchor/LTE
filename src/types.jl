#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor


struct LTEBase
    seis_data :: Array{Float64}
    disp_data :: Union{Array{Float64, 1},Nothing}
    freq :: Array{Float64}
    fqminmax :: Tuple{Integer, Integer}
    fq_band :: Tuple{Float64, Float64}
    fs :: Integer
    nfs :: Integer
    NW :: Float64
    pad :: Float64
    nwin :: Integer
    lwin :: Integer
    nadv :: Float64
    nini :: Integer
    dsar :: Bool
    opt :: Bool
    polar :: Bool
    pe_order :: Integer
    pe_delta :: Integer
    opt_twin :: Float64
    opt_th :: Float64
end


struct Angles
    azimuth :: Float64
    elev :: Float64
    phyHH :: Float64
    phyVH :: Float64
end


struct PParams
    freq :: Union{Array{Float64}, Nothing}
    degree :: Array{Float64}
    rect :: Array{Float64}
    azimuth :: Array{Float64}
    elev  :: Array{Float64}
    phyHH :: Array{Float64}
    phyVH :: Array{Float64}
end


struct SParams
    freq :: Array{Float64}
    S :: Array{Float64}
    erg :: Float64
    dominant :: Float64
    centroid  :: Float64
end
