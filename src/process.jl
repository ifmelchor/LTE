#!/usr/local/bin julia
# coding=utf-8

# Ivan Melchor (ifmelchor@unrn.edu.ar)


"""
   sta_run(*args)

Funcion LTE-station para calculcar espectrograma, polargrama, etc...
"""
function sta_run(data::Array{T,2}, channels::Tuple, fs::J, nwin::J, lwin::J, wadv::T, nswin::Union{J,Nothing}, lswin::Union{J,Nothing}, swadv::Union{T,Nothing}, fq_band::Vector{T}, NW::T, pad::T, add_param::Bool, polar::Bool, pe_order::J, pe_delta::J, ap_twin::T, ap_th::T) where {T<:Real, J<:Integer}

    if length(channels) != size(data, 1)
        throw(ArgumentError("nro of components and seismic-data-matrix raw's must be equal"))
    end

    if polar && length(channels) != 3
        throw(ArgumentError("for polar analysis, the nro of components must be three"))
    end

    # compute the frequency domain
    if !isnothing(lswin)
        freq, fqminmax = _fqbds(lswin, fs, fq_band, pad=pad)
    else
        freq, fqminmax = _fqbds(lwin, fs, fq_band, pad=pad)
    end
    
    # define base
    base = Base(fs, fq_band, fqminmax, nswin, lswin, swadv, NW, pad, pe_order, pe_delta, ap_twin, ap_th)
    
    # create empty dict
    nfs = size(freq, 1)
    lte_dict = _empty_dict(channels, add_param, polar, nfs, nwin)

    for n in 1:nwin
        n0 = 1 + floor(Int32, lwin * wadv*(n-1))
        nf = floor(Int32, n0 + lwin)
        data_n = @view data[:, n0:nf]

        stop = false
        for (c, chan) in enumerate(channels)
            if add_param && c==1
                compute_param = true
            else
                compute_param = false
            end

            try
                cspe, perm_entr, copt = _stacore(data_n[c, :], base, compute_param)
                lte_dict[chan]["specgram"][n,:]  = cspe.S
                lte_dict[chan]["energy"][n]      = cspe.erg
                lte_dict[chan]["fq_dominant"][n] = cspe.dominant
                lte_dict[chan]["fq_centroid"][n] = cspe.centroid
                lte_dict[chan]["perm_entr"][n]   = perm_entr
                if compute_param
                    lte_dict["opt"]["vlf"][n]  = copt.vlf
                    lte_dict["opt"]["lf"][n]   = copt.lf
                    lte_dict["opt"]["vlar"][n] = copt.vlar
                    lte_dict["opt"]["rsam"][n] = copt.rsam
                    lte_dict["opt"]["lrar"][n] = copt.lrar
                    lte_dict["opt"]["mf"][n]   = copt.mf
                    lte_dict["opt"]["rmar"][n] = copt.rmar
                    lte_dict["opt"]["hf"][n]   = copt.hf
                    lte_dict["opt"]["dsar"][n] = copt.dsar  
                end
            catch e
                stop = true
                lte_dict[chan]["specgram"][n,:]  .= NaN
                lte_dict[chan]["energy"][n]      = NaN
                lte_dict[chan]["fq_dominant"][n] = NaN
                lte_dict[chan]["fq_centroid"][n] = NaN
                lte_dict[chan]["perm_entr"][n]   = NaN
                if compute_param
                    lte_dict["opt"]["vlf"][n]  = NaN
                    lte_dict["opt"]["lf"][n]   = NaN
                    lte_dict["opt"]["vlar"][n] = NaN
                    lte_dict["opt"]["rsam"][n] = NaN
                    lte_dict["opt"]["lrar"][n] = NaN
                    lte_dict["opt"]["mf"][n]   = NaN
                    lte_dict["opt"]["rmar"][n] = NaN
                    lte_dict["opt"]["hf"][n]   = NaN
                    lte_dict["opt"]["dsar"][n] = NaN
                end
            end
        end
        
        if polar
            if stop
                lte_dict["polar"]["degree"][n,:]  .= NaN 
                lte_dict["polar"]["rect"][n,:]    .= NaN 
                lte_dict["polar"]["azimuth"][n,:] .= NaN 
                lte_dict["polar"]["elev"][n,:]    .= NaN 
                lte_dict["polar"]["phyhh"][n,:]   .= NaN 
                lte_dict["polar"]["phyvh"][n,:]   .= NaN
            else
                try
                    polarP = _polar(data_n[:,:], base)
                    lte_dict["polar"]["degree"][n,:]  = polarP.degree
                    lte_dict["polar"]["rect"][n,:]    = polarP.rect
                    lte_dict["polar"]["azimuth"][n,:] = polarP.azimuth
                    lte_dict["polar"]["elev"][n,:]    = polarP.elev
                    lte_dict["polar"]["phyhh"][n,:]   = polarP.phyhh
                    lte_dict["polar"]["phyvh"][n,:]   = polarP.phyvh
                catch e
                    lte_dict["polar"]["degree"][n,:]  .= NaN 
                    lte_dict["polar"]["rect"][n,:]    .= NaN 
                    lte_dict["polar"]["azimuth"][n,:] .= NaN 
                    lte_dict["polar"]["elev"][n,:]    .= NaN 
                    lte_dict["polar"]["phyhh"][n,:]   .= NaN 
                    lte_dict["polar"]["phyvh"][n,:]   .= NaN
                end
            end
        end
    
    end
    
    return lte_dict
end


"""
   polar_run(*args)
Funcion para calcular el polargrama
"""
function polar_run(data::Array{T,2}, fq_band::Vector{T}, fs::J, NW::T, pad::T, nswin::Union{J,Nothing}, lswin::Union{J,Nothing}, nadv::Union{T,Nothing}, return_all::Bool, full_return::Bool) where {T<:Real, J<:Integer}

    if size(data, 1) != 3
        throw(ArgumentError("the nro of components must be three"))
    end

    # compute the frequency domain
    if !isnothing(lswin)
        freq, fqr = _fqbds(lswin, fs, fq_band, pad=pad)
    else
        lwin = size(data, 2)
        freq, fqr = _fqbds(lwin, fs, fq_band, pad=pad)
    end

    base = Base(fs, fq_band, fqr, nswin, lswin, nadv, NW, pad, nothing, nothing, nothing, nothing)

    return freq, _polar(data, base, return_all=return_all, full_return=full_return)
end


"""
   csw_run(*args)
Funcion para calcular el spectral width
"""
function csw_run(data::Array{T,2}, fq_band::Vector{T}, fs::J, NW::T, pad::T, nswin::J, lswin::J, nadv::T, return_vt::Bool; win_freq::T=0.33) where {T<:Real, J<:Integer}

    nsta = size(data, 1)
    if nsta < 3
        throw(ArgumentError("the minimum nro of components must be three"))
    end

    # compute the frequency domain and define Base param
    freq, fqr = _fqbds(lswin, fs, fq_band, pad=pad)
    nfs = size(freq, 1)

    # preprocess
    _spec_white!(data, fs, win_freq, fq_band)

    # compute CSM
    csm_svd = _csm(data, fs, lswin, nswin, nadv, fqr, NW, pad)
    csw = Array{Float32}(undef, nfs)

    if return_vt
        vt = Array{ComplexF32}(undef, nfs, nsta)
    end

    for i in eachindex(csm_svd)
        s   = csm_svd[i]
        csw[i] = dot(0:nsta-1, s.S) / sum(s.S) # spectral width
        
        if return_vt
            vt[i,:] = s.Vt[1,:] # first eigenvector
        end
    end

    if return_vt
        return freq, csw, vt
    else
        return freq, csw
    end
end


"""
   net_run(*args)

Funcion LTE-network para calculcar espectrograma y CSW (covariance spectral width) sugerido por Seydeux
"""
function net_run(data::Array{T,2}, channels::Tuple, fs::J, nwin::J, lwin::J, wadv::T, nswin::J, lswin::J, swadv::T, fq_band::Vector{T}, NW::T, pad::T) where {T<:Real, J<:Integer}

    if length(channels) != size(data, 1)
        throw(ArgumentError("nro of components and seismic-data-matrix raw's must be equal"))
    end

    # compute the frequency domain
    freq, fqminmax = _fqbds(lswin, fs, fq_band, pad=pad)
    
    # define base srtuct
    base = Base(fs, fq_band, fqminmax, nswin, lswin, nadv, NW, pad, nothing, nothing, nothing, nothing)

    # create empty dict
    nfs = size(freq, 1)
    lte_dict = _empty_dict(channels, nfs, nwin)

    for n in 1:nwin
        n0 = 1 + floor(Int32, lwin * wadv*(n-1))
        nf = floor(Int32, n0 + lwin)
        data_n = @view data[:, n0:nf]

        # compute parameters
        sxx_dict, csw_n, vt_n = _netcore(data_n, channels, base)

        # save into lte_dict
        for chan in channels
            lte_dict[chan][n,:] = sxx_dict[chan]
        end

        lte_dict["csw"][n,:] = csw_n
        lte_dict["vt"][n,:,:] = vt_n
    end

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
function _netcore(data::AbstractArray{T}, channels::Tuple, base::Base) where T<:Real

    # save specgram
    sxx_dict = Dict()
    for (c, chan) in enumerate(channels)
        sxx, _ = _psd(data[c,:], base.fs, base.lswin, base.nswin, base.nadv, base.fqminmax, base.NW, base.pad)
        sxx_dict[chan] = sxx
    end

    # preprocess
    _spec_white!(data, base.fs, 0.33, base.fq_band)

    # compute CSM
    csm_svd = _csm(data, base.fs, base.lswin, base.nswin, base.nadv, base.fqminmax, base.NW, base.pad)
    
    # define empty arrays
    nfs  = length(csm_svd)
    nsta = length(channels)
    csw_n = Array{Float32}(undef, nfs)
    vt_n  = Array{ComplexF32}(undef, nfs, nsta) 
    
    # compute spectral width and Vt
    for i in eachindex(csm_svd)
        s   = csm_svd[i]
        csw_n[i] = dot(0:nsta-1, s.S) / sum(s.S) # spectral width
        vt_n[i,:] = s.Vt[1,:] # first eigenvector
    end

    return sxx_dict, csw_n, vt_n
end