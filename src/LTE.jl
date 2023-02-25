#!/usr/local/bin julia
# coding=utf-8

# LTE package for SEISVO software

# Ivan Melchor

__precompile__()


module LTE

    using LinearAlgebra
    using ComplexValues
    using Multitaper
    
    export sta_run, net_run, polar_run

    include("types.jl")

    include("utils.jl")

    include("spectral.jl")

    include("polar.jl")

    include("process.jl")

end
