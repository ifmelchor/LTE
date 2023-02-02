#!/usr/local/bin julia
# coding=utf-8

# Ivan Melchor


"""
   sta_lte_run(*args)

Funcion LTE-station para calculcar espectrograma, polargrama, etc...
"""
function lte_run(s_data::Array{Float64,2}, d_data::Union{Array{Float64,1}, Nothing}, channels::Vector{String}, fs::Integer, nwin::Integer, lwin::Integer, nsubwin::Union{Integer,Nothing}, lsubwin::Union{Integer,Nothing}, nadv::Union{Float64,Nothing}, fq_band::Tuple{Float64,Float64}, NW::Float64, pad::Float64, add_param::Bool, polar::Bool, pe_order::Integer, pe_delta::Integer, ap_twin::Float64, ap_th::Float64)

    # channel info for polarization analysis:
    #  1 --> Z
    #  2 --> N
    #  3 --> E

    if size(channels, 1) != size(s_data, 1)
        throw(ArgumentError("nro of components and seismic-data-matrix raw's must be equal"))
    end

    # compute the frequency domain
    if lswin
        freq, fqr = _fqbds(lswin, fs, fq_band, pad=pad)
    else
        freq, fqr = _fqbds(lwin, fs, fq_band, pad=pad)
    end
    
    # define base
    nfs = size(freq, 1)
    base = LTEBase(s_data, d_data, channels, fs, freq, fq_band, fqr, nfs, nwin, lwin, nsubwin, lsubwin, nadv, NW, pad, add_param, polar, pe_order, pe_delta, ap_twin, ap_th)
    
    # run lte
    lte = _starun(base)

    return lte
end


"""
   _empty_lte(*args)

Genera un dict vacio para llenar durante el procesado.
"""
function _empty_lte(base::LTEBase)
    dict = Dict()
    
    for chan in base.channels
        dict[chan] = Dict()
        dict[chan]["specgram"]    = Array{Float64}(undef, base.nwin, base.nfs)
        dict[chan]["energy"]      = Array{Float64}(undef, base.nwin)
        dict[chan]["fq_dominant"] = Array{Float64}(undef, base.nwin)
        dict[chan]["fq_centroid"] = Array{Float64}(undef, base.nwin)
        dict[chan]["perm_entr"]   = Array{Float64}(undef, base.nwin)
    end


    if base.add_param
        dict["opt"] = Dict()
        for attr in ("vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf", "dsar")
            dict["opt"][attr] = Array{Float64}(undef, base.nwin)
        end
    end


    if base.polar
        dict["polar"] = Dict()
        dict["polar"]["degree"]   = Array{Float64}(undef, base.nwin, base.nfs)
        dict["polar"]["rect"]     = Array{Float64}(undef, base.nwin, base.nfs)
        dict["polar"]["azimuth"]  = Array{Float64}(undef, base.nwin, base.nfs)
        dict["polar"]["elev"]     = Array{Float64}(undef, base.nwin, base.nfs)
        dict["polar"]["phyhh"]    = Array{Float64}(undef, base.nwin, base.nfs)
        dict["polar"]["phyvh"]    = Array{Float64}(undef, base.nwin, base.nfs)
    end

    return dict
end


"""
   _starun(*args)

Core del procesamiento LTE
"""
function _starun(base::LTEBase)
    dlte = _empty_lte(base) # define an empty dictionary

    for n in 1:base.nwin
        n0 = floor(Int, base.lwin*(n-1))
        nf = floor(Int, n0 + base.lwin)
        sdata_n = @view base.seis_data[:, n0:nf]

        if !isnothing(base.disp_data)
            ddata_n = @view base.disp_data[n0:nf]
        else
            ddata_n = nothing
        end
        
        for (c, chan) in enumerate(base.channels)
            if base.add_param & c==1
                cspe, perm_entr, copt = _core(sdata_n[c, :], ddata_n, base, true)
                dict["opt"]["vlf"]  = copt.vlf
                dict["opt"]["lf"]   = copt.lf
                dict["opt"]["vlar"] = copt.vlar
                dict["opt"]["rsam"] = copt.rsam
                dict["opt"]["lrar"] = copt.lrar
                dict["opt"]["mf"]   = copt.mf
                dict["opt"]["rmar"] = copt.rmar
                dict["opt"]["hf"]   = copt.hf
                dict["opt"]["dsar"] = copt.dsar
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
        if c==3 & base.polar
            polar = _polar(sdata_n, base)
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
function _core(s_data::Array{Float64,1}, d_data::Union{Array{Float64,1}}, base::LTEBase, opt_params::Bool)

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
    else:
        opt_p = nothing
    end

    return spec_p, perm_entr, opt_p
end