#!/usr/local/bin julia
# coding=utf-8

# Ivan Melchor (ifmelchor@unrn.edu.ar)


"""
   sta_run(*args)

Funcion LTE-station para calculcar espectrograma, polargrama, etc...
"""
function sta_run(data::Array{T,2}, channels::String, fs::J, nwin::J, lwin::J, nswin::Union{J,Nothing}, lswin::Union{J,Nothing}, nadv::Union{T,Nothing}, fq_band::Tuple{T,T}, NW::T, pad::T, add_param::Bool, polar::Bool, pe_order::J, pe_delta::J, ap_twin::T, ap_th::T) where {T<:Real, J<:Integer}

    # convert channels to a Vector of strings and fq_band to a vector of floats
    channels =[String(ch) for ch in split(channels, "/")]
    fq_band  = collect(fq_band)

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
    base = Base(fs, fq_band, fqminmax, nswin, lswin, nadv, NW, pad, pe_order, pe_delta, ap_twin, ap_th)
    
    # create empty dict
    nfs = size(freq, 1)
    lte_dict = _empty_dict(channels, add_param, polar, nfs, nwin)

    for n in 1:nwin
        n0 = 1 + floor(Int32, lwin*(n-1))
        nf = floor(Int32, n0 + lwin)
        data_n = @view data[:, n0:nf]

        for (c, chan) in enumerate(channels)
            if add_param && c==1
                cspe, perm_entr, copt = _stacore(data_n[c, :], base, true)
                lte_dict["opt"]["vlf"][n]  = copt.vlf
                lte_dict["opt"]["lf"][n]   = copt.lf
                lte_dict["opt"]["vlar"][n] = copt.vlar
                lte_dict["opt"]["rsam"][n] = copt.rsam
                lte_dict["opt"]["lrar"][n] = copt.lrar
                lte_dict["opt"]["mf"][n]   = copt.mf
                lte_dict["opt"]["rmar"][n] = copt.rmar
                lte_dict["opt"]["hf"][n]   = copt.hf
                lte_dict["opt"]["dsar"][n] = copt.dsar
            else
                cspe, perm_entr, _ = _stacore(data_n[c, :], base, false)
            end

            lte_dict[chan]["specgram"][n,:]  = cspe.S
            lte_dict[chan]["energy"][n]      = cspe.erg
            lte_dict[chan]["fq_dominant"][n] = cspe.dominant
            lte_dict[chan]["fq_centroid"][n] = cspe.centroid
            lte_dict[chan]["perm_entr"][n]   = perm_entr
        end
        
        # only when three components are available
        if  polar
            polarP = _polar(data_n[:,:], base)
            lte_dict["polar"]["degree"][n,:]  = polarP.degree
            lte_dict["polar"]["rect"][n,:]    = polarP.rect
            lte_dict["polar"]["azimuth"][n,:] = polarP.azimuth
            lte_dict["polar"]["elev"][n,:]    = polarP.elev
            lte_dict["polar"]["phyhh"][n,:]   = polarP.phyhh
            lte_dict["polar"]["phyvh"][n,:]   = polarP.phyvh
        end
    
    end
    
    return lte_dict
end


"""
   polar_run(*args)
Funcion para calcular el polargrama
"""
function polar_run(data::Array{T,2}, fq_band::Tuple{T,T}, fs::J, NW::T, pad::T, nswin::Union{J,Nothing}, lswin::Union{J,Nothing}, nadv::Union{T,Nothing}, full_return::Bool) where {T<:Real, J<:Integer}

    fq_band  = collect(fq_band)

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

    pe_order = pe_delta = ap_twin = ap_th = nothing
    base = Base(fs, fq_band, fqminmax, nswin, lswin, nadv, NW, pad, pe_order, pe_delta, ap_twin, ap_th)

    return freq, _polar(data, base, full_return=full_return)
end


"""
   net_run(*args)

Funcion LTE-network para calculcar espectrograma y CSW (covariance spectral width) sugerido por Seydeux
"""
function net_run(data::Array{T,2}, channels::String, fs::J, nswin::J, lswin::J, nadv::T, fq_band::Tuple{T,T}, NW::T, pad::T) where {T<:Real, J<:Integer}

    channels =[String(ch) for ch in split(channels, "/")]
    fq_band  = collect(fq_band)
    
    nsta = size(data, 1)

    if size(channels, 1) != nsta
        throw(ArgumentError("nro of components and seismic-data-matrix raw's must be equal"))
    end

    # compute the frequency domain
    freq, fqminmax = _fqbds(lswin, fs, fq_band, pad=pad)
    
    pe_order = pe_delta = ap_twin = ap_th = nothing
    base = Base(fs, fq_band, fqminmax, nswin, lswin, nadv, NW, pad, pe_order, pe_delta, ap_twin, ap_th)

    # create empty dict
    nfs = size(freq, 1)
    lte_dict = _empty_dict(channels, nfs)

    # compute parameters
    _netcore(data, channels, base, lte_dict)

    return lte_dict
end


"""
    _stacore(args)

Compute core parameters for LTE
"""
function _stacore(data::AbstractArray{T}, base::Base, opt_params::Bool) where T<:Real

    # compute psd
    sxx, freq = _psd(data, base.fs, base.lswin, base.nswin, base.nadv, base.fqminmax, base.NW, base.pad)

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
function _netcore(data::AbstractArray{T}, channels::Vector{String}, base::Base, lte_dict::Dict) where T<:Real

    # save specgram
    for (c, chan) in enumerate(channels)
        sxx, _ = _psd(data[c,:], base.fs, base.lswin, base.nswin, base.nadv, base.fqminmax, base.NW, base.pad)
        lte_dict[chan] = sxx
    end

    # preprocess
    _spec_white!(data, base.fs, 0.33, base.fq_band)

    # compute CSM
    csm_svd = _csm(data, base.fs, base.lswin, base.nswin, base.nadv, base.fqminmax, base.NW, base.pad)
    
    # compute spectral width and save Vt
    nsta = size(channels,1)
    idx = 0:nsta-1
    for i in eachindex(csm_svd)
        s   = csm_svd[i]
        # spectral width
        lte_dict["csw"][i] = dot(idx, s.S)/sum(s.S)
        # first eigenvector
        lte_dict["vt"][i,:] = s.Vt[1,:]
    end

end