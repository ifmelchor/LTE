#!/usr/local/bin julia
# coding=utf-8

# Ivan Melchor


"""
   sta_lte_run(*args)

Funcion LTE-station para calculcar espectrograma, polargrama, etc...
"""
function lte_run(s_data::Array{Float64,2}, d_data::Union{Array{Float64,1}, Nothing}, channels::Vector{String}, fs::Int64, nwin::Int64, lwin::Int64, nswin::Union{Int64,Nothing}, lswin::Union{Int64,Nothing}, nadv::Union{Float64,Nothing}, fq_band::Vector{Float64}, NW::Float64, pad::Float64, add_param::Bool, polar::Bool, pe_order::Int64, pe_delta::Int64, ap_twin::Float64, ap_th::Float64)

    # channel info for polarization analysis:
    #  1 --> Z
    #  2 --> N
    #  3 --> E

    if size(channels, 1) != size(s_data, 1)
        throw(ArgumentError("nro of components and seismic-data-matrix raw's must be equal"))
    end

    # compute the frequency domain
    if !isnothing(lswin)
        freq, fqr = _fqbds(lswin, fs, fq_band, pad=pad)
    else
        freq, fqr = _fqbds(lwin, fs, fq_band, pad=pad)
    end
    
    # define base
    nfs = size(freq, 1)
    base = LTEBase(s_data, d_data, channels, fs, freq, fq_band, fqr, nfs, nwin, lwin, nswin, lswin, nadv, NW, pad, add_param, polar, pe_order, pe_delta, ap_twin, ap_th)
    
    # run lte
    lte = _starun(base)

    return lte
end


"""
   _empty_lte(*args)

Genera un dict vacio para llenar durante el procesado.
"""
function _empty_dict(base::LTEBase)
    dict = Dict()
    
    for chan in base.channels
        dict[chan] = Dict()
        dict[chan]["specgram"] = Array{Float64}(undef, base.nwin, base.nfs)
        for attr in ("energy", "fq_dominant", "fq_centroid", "perm_entr")
            dict[chan][attr] = Array{Float64}(undef, base.nwin)
        end
    end

    if base.add_param
        dict["opt"] = Dict()
        for attr in ("vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf", "dsar")
            dict["opt"][attr] = Array{Float64}(undef, base.nwin)
        end
    end

    if length(base.channels)==3 && base.polar
        dict["polar"] = Dict()
        for attr in ("degree", "rect", "azimuth", "elev", "phyhh", "phyvh")
            dict["polar"][attr] = Array{Float64}(undef, base.nwin, base.nfs)
        end
    end

    return dict
end


"""
   _starun(*args)

Core del procesamiento LTE
"""
function _starun(base::LTEBase)
    dict = _empty_dict(base) # define an empty dictionary

    for n in 1:base.nwin
        n0 = 1 + floor(Int64, base.lwin*(n-1))
        nf = floor(Int64, n0 + base.lwin)
        sdata_n = @view base.seis_data[:, n0:nf]

        if !isnothing(base.disp_data)
            ddata_n = @view base.disp_data[n0:nf]
        else
            ddata_n = nothing
        end
        
        for (c, chan) in enumerate(base.channels)
            if base.add_param && c==1
                cspe, perm_entr, copt = _core(sdata_n[c, :], ddata_n, base, true)
                dict["opt"]["vlf"][n]  = copt.vlf
                dict["opt"]["lf"][n]   = copt.lf
                dict["opt"]["vlar"][n] = copt.vlar
                dict["opt"]["rsam"][n] = copt.rsam
                dict["opt"]["lrar"][n] = copt.lrar
                dict["opt"]["mf"][n]   = copt.mf
                dict["opt"]["rmar"][n] = copt.rmar
                dict["opt"]["hf"][n]   = copt.hf
                if !isnothing(base.disp_data)
                    dict["opt"]["dsar"][n] = copt.dsar
                end
                
            else
                cspe, perm_entr, copt = _core(sdata_n[c, :], nothing, base, false)
            end

            dict[chan]["specgram"][n,:]  = cspe.S
            dict[chan]["energy"][n]      = cspe.erg
            dict[chan]["fq_dominant"][n] = cspe.dominant
            dict[chan]["fq_centroid"][n] = cspe.centroid
            dict[chan]["perm_entr"][n]   = perm_entr
        end
        
        # only when three components are available
        if  base.polar
            polar = _polar(sdata_n[:,:], base)
            dict["polar"]["degree"][n,:]  = polar.degree
            dict["polar"]["rect"][n,:]    = polar.rect
            dict["polar"]["azimuth"][n,:] = polar.azimuth
            dict["polar"]["elev"][n,:]    = polar.elev
            dict["polar"]["phyhh"][n,:]   = polar.phyhh
            dict["polar"]["phyvh"][n,:]   = polar.phyvh
        end
    
    end
    
    return dict
end


"""
    _core(args)

Compute core parameters for LTE
"""
function _core(s_data::Array{Float64,1}, d_data::Union{Array{Float64,1}, Nothing}, base::LTEBase, opt_params::Bool)

    # compute psd
    sxx, freq = _psd(s_data, base.fs, base.lswin, base.nswin, base.nadv, base.fqminmax, base.nfs, base.NW, base.pad)

    # dominant freq
    domfreq = freq[findmax(sxx)[2]]
            
    # centroid freq
    cenfreq = dot(freq,sxx)/sum(sxx)
            
    # energy
    erg = sum(sxx)

    # crate spectral struct
    spec_p = SParams(freq, sxx, erg, domfreq, cenfreq)

    # permutation entropy
    perm_entr = _PE(s_data, base.pe_order, base.pe_delta, base.fq_band, base.fs)

    if opt_params
        opt_p = _optparams(s_data, d_data, base.fs, twindow=base.ap_twin, threshold=base.ap_th)
    else
        opt_p = nothing
    end

    return spec_p, perm_entr, opt_p
end