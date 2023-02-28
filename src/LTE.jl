#!/usr/local/bin julia
# coding=utf-8

# LTE package for SEISVO software

# Ivan Melchor

__precompile__()


module LTE

    using LinearAlgebra
    using ComplexValues
    using Multitaper
    using Statistics
    
    export sta_run, net_run, polar_run

    include("types.jl")

    include("dicts.jl")

    include("utils.jl")

    include("more_params.jl")

    include("spectral.jl")

    include("polar.jl")

    include("process.jl")

end
