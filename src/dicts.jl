#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor


"""
   _empty_lte(*args)

Genera un dict vacio para llenar durante el procesado.
"""
function _empty_dict(channels::Tuple, add_param::Bool, polar::Bool, nfs::J, nwin::J) where J<:Integer

    dict = Dict()
    
    for chan in channels
        dict[chan] = Dict()
        dict[chan]["specgram"] = Array{Float32}(undef, nwin, nfs)
        for attr in ("energy", "fq_dominant", "fq_centroid", "perm_entr")
            dict[chan][attr] = Array{Float32}(undef, nwin)
        end
    end

    if add_param
        dict["opt"] = Dict()
        for attr in ("vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf", "dsar")
            dict["opt"][attr] = Array{Float32}(undef, nwin)
        end
    end

    if polar
        dict["polar"] = Dict()
        for attr in ("degree", "rect", "azimuth", "elev", "phyhh", "phyvh")
            dict["polar"][attr] = Array{Float32}(undef, nwin, nfs)
        end
    end

    return dict
end


function _empty_dict(channels::Tuple, nfs::J, nwin::J) where J<:Integer
    dict = Dict()
    
    nsta = length(channels)
    
    for chan in channels
        dict[chan] = Array{Float32}(undef, nwin, nfs)
    end

    dict["csw"]   = Array{Float32}(undef, nwin, nfs)
    dict["vt"]    = Array{ComplexF32}(undef, nwin, nfs, nsta)

    return dict
end
