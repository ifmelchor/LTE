#!/usr/local/bin julia
# coding=utf-8

# Ivan Melchor


"""
   sta_run(*args)

Funcion LTE-station para calculcar espectrograma, polargrama, etc...
"""
function sta_run(data::Array{T,2}, channels::Vector{String}, fs::J, nwin::J, lwin::J, nswin::Union{J,Nothing}, lswin::Union{J,Nothing}, nadv::Union{T,Nothing}, fq_band::Vector{T}, NW::T, pad::T, add_param::Bool, polar::Bool, pe_order::J, pe_delta::J, ap_twin::T, ap_th::T) where {T<:Real, J<:Int}

    # channel info for polarization analysis:
    #  1 --> Z
    #  2 --> N
    #  3 --> E

    if size(channels, 1) != size(data, 1)
        throw(ArgumentError("nro of components and seismic-data-matrix raw's must be equal"))
    end

    # compute the frequency domain
    if !isnothing(lswin)
        freq, fqminmax = _fqbds(lswin, fs, fq_band, pad=pad)
    else
        freq, fqminmax = _fqbds(lwin, fs, fq_band, pad=pad)
    end
    
    # define base
    base = STABase(fs, fq_band, fqminmax, size(freq, 1), nwin, lwin, nswin, lswin, nadv, NW, pad, pe_order, pe_delta, ap_twin, ap_th)
    
    # create empty dict
    lte_dict = _empty_dict(channels, add_param, polar, base)

    for n in 1:nwin
        n0 = 1 + floor(Int64, lwin*(n-1))
        nf = floor(Int64, n0 + lwin)
        data_n = @view data[:, n0:nf]
        
        for (c, chan) in enumerate(channels)
            if add_param && c==1
                cspe, perm_entr, copt = _stacore(data_n[c, :], base, true)
                dict["opt"]["vlf"][n]  = copt.vlf
                dict["opt"]["lf"][n]   = copt.lf
                dict["opt"]["vlar"][n] = copt.vlar
                dict["opt"]["rsam"][n] = copt.rsam
                dict["opt"]["lrar"][n] = copt.lrar
                dict["opt"]["mf"][n]   = copt.mf
                dict["opt"]["rmar"][n] = copt.rmar
                dict["opt"]["hf"][n]   = copt.hf
                dict["opt"]["dsar"][n] = copt.dsar
            else
                cspe, perm_entr, _ = _stacore(data_n[c, :], base, false)
            end

            dict[chan]["specgram"][n,:]  = cspe.S
            dict[chan]["energy"][n]      = cspe.erg
            dict[chan]["fq_dominant"][n] = cspe.dominant
            dict[chan]["fq_centroid"][n] = cspe.centroid
            dict[chan]["perm_entr"][n]   = perm_entr
        end
        
        # only when three components are available
        if  polar
            polarP = _polar(data_n[:,:], base)
            dict["polar"]["degree"][n,:]  = polarP.degree
            dict["polar"]["rect"][n,:]    = polarP.rect
            dict["polar"]["azimuth"][n,:] = polarP.azimuth
            dict["polar"]["elev"][n,:]    = polarP.elev
            dict["polar"]["phyhh"][n,:]   = polarP.phyhh
            dict["polar"]["phyvh"][n,:]   = polarP.phyvh
        end
    
    end
    
    return lte_dict
end


"""
   polar_run(*args)
Funcion para calcular el polargrama
"""
function polar_run(data::Array{T,2}, fq_band::Vector{T}, fs::J, NW::T, pad::T, nswin::Union{J,Nothing}, lswin::Union{J,Nothing}, nadv::Union{T,Nothing}, full_return::Bool) where {T<:Real, J<:Int}

    if size(data, 1) != 3
        throw(ArgumentError("nro of components must be 3 for polarization analysis"))
    end

    # compute the frequency domain
    if !isnothing(lswin)
        freq, fqminmax = _fqbds(lswin, fs, fq_band, pad=pad)
    else
        lwin = size(data, 2)
        freq, fqminmax = _fqbds(lwin, fs, fq_band, pad=pad)
    end

    nfs = size(freq, 1)
    base = PolarBase(fs, fq_band, fqminmax, nfs, nswin, lswin, nadv, NW, pad)

    return freq, _polar(data, base, full_return=full_return)
end


"""
    _stacore(args)

Compute core parameters for LTE
"""
function _stacore(data::AbstractArray{T}, base::STABase, opt_params::Bool) where T<:Real

    # compute psd
    sxx, freq = _psd(data, base.fs, base.lswin, base.nswin, base.nadv, base.fqminmax, base.nfs, base.NW, base.pad)

    # dominant freq
    domfreq = freq[findmax(sxx)[2]]
            
    # centroid freq
    cenfreq = dot(freq,sxx)/sum(sxx)
            
    # energy
    erg = sum(sxx)

    # crate spectral struct
    spec_p = SParams(sxx, erg, domfreq, cenfreq)

    # permutation entropy
    perm_entr = _PE(data, base.pe_order, base.pe_delta, base.fq_band, base.fs)

    if opt_params
        opt_p = _optparams(data, base.fs, base.ap_twin, base.ap_th)
    else
        opt_p = nothing
    end

    return spec_p, perm_entr, opt_p
end



"""
    _netcore(args)

Compute core parameters for LTE
"""
function _netcore(data::AbstractArray{T}, base::STABase, opt_params::Bool) where T<:Real

    # compute psd
    sxx, _ = _psd(data, base.fs, base.lswin, base.nswin, base.nadv, base.fqminmax, base.nfs, base.NW, base.pad)
            
    # energy
    erg = sum(sxx)
    

    return spec_p, perm_entr, opt_p
end