#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor


function lte_run(s_data::Array{Float64,1}, d_data::Union{Array{Float64,1}, Nothing}, fs::Integer, nwin::Integer, lwin::Integer, nini::Integer, fq_band::Tuple{Float64,Float64}, NW::Float64, pad::Float64, dsar::Bool, opt::Bool, polar::Bool, pe_order::Integer, pe_delta::Integer, opt_twin::Float64, opt_th::Float64)

    # comp info:
    #  1 --> Z
    #  2 --> 
    #  3 --> 

    ncomp = size(base.s_data)[1]
    ndata = size(base.s_data)[2]
    freq, fqr = _fqbds(ndata, fs, fq_band, pad)
    nfs = size(freq)[1]

    # define base
    base = LTEBase(s_data, d_data, freq, fqr, fq_band, fs, nfs, NW, pad, nwin, lwin, nadv, nini, dsar, opt, polar, pe_order, pe_delta, opt_twin, opt_th)
    lte = _run(base)

    return lte
end


function polargram(data::Array{Float64,2}, fs::Integer, nwin::Integer, lwin::Integer, nini::Integer, fq_band::Tuple{Float64,Float64}; NW=3.5, pad=1.0)

    ndata = size(base.s_data)[2]
    freq, fqr = _fqbds(ndata, fs, fq_band, pad)

    degree = Array{Float64}(undef, nwin, base.nfs)
    rect = Array{Float64}(undef, nwin, base.nfs)
    azimuth = Array{Float64}(undef, nwin, base.nfs)
    elev = Array{Float64}(undef, nwin, base.nfs)
    phyHH = Array{Float64}(undef, nwin, base.nfs)
    phyVH = Array{Float64}(undef, nwin, base.nfs)

    for n in 1:base.nwin
        ninikk = base.nini + base.nadv*(n-1)
        data_n = @view base.data[:, ninikk:ninikk+base.lwin]
        polar = _polar(data_n, fs, fqr, fq_band, NW=base.NW, pd=base.pad)
        degree[n,:] = polar.degree
        azimuth[n,:] = polar.azimuth
        elev[n,:] = polar.elev
        phyHH[n,:] = polar.phyHH
        phyVH[n,:] = polar.phyVH
    end
    
    return PParams(freq, degree, rect, azimuth, elev, phyHH, phyVH)
end


function _empty_lte(base::LTEBase)
    dict = Dict("spec"=>Dict(), "npe"=>Dict())
    for c in 1:base.ncomp
        dict["spec"][c] = Array{Float64}(undef, base.nwin, base.nfs)
        dict["erg"][c]  = Array{Float64}(undef, base.nwin)
        dict["dfq"][c]  = Array{Float64}(undef, base.nwin)
        dict["cfq"][c]  = Array{Float64}(undef, base.nwin)
        dict["pe"][c]   = Array{Float64}(undef, base.nwin)
    end

    if base.dsar
        dict["dsar"] = Array{Float64}(undef, base.nwin)
    end

    if base.opt
        for attr in ("vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf")
            dict[attr] = Array{Float64}(undef, base.nwin)
        end
    end

    if base.polar
        dict["polar"] = Dict()
        dict["polar"]["degree"] = Array{Float64}(undef, base.nwin, base.nfs)
        dict["polar"]["rect"] = Array{Float64}(undef, base.nwin, base.nfs)
        dict["polar"]["azimuth"] = Array{Float64}(undef, base.nwin, base.nfs)
        dict["polar"]["elev"] = Array{Float64}(undef, base.nwin, base.nfs)
    end

    return dict
end


function _run(base::LTEBase)
    dlte = _empty_lte(base) # build empty dictionary

    for n in 1:base.nwin
        ninikk = base.nini + base.nadv*(n-1)
        
        for c in 1:base.ncomp
            cdata = @view base.sdata[c, ninikk:ninikk+base.lwin]
            cspe = _spectral(cdata, base.fs, base.fqminmax, base.fq_band, NW=base.NW, pd=base.pad)
            dlte["spec"][c][n,:] = cspe.S
            dlte["erg"][c][n] = cspe.erg
            dlte["dfq"][c][n] = cspe.dominant
            dlte["cfq"][c][n] = cspe.centroid
            dlte["pe"][c][n] = _PE(cdata, base.pe_order, base.pe_delta, base.fq_band, base.fs)
            
            if c==1 # comp 1 --> Z
                if base.dsar 
                    dict["dsar"][n] = _dsar(cdata, base.fs)
                end

                if base.opt
                    copt = _optparams(cdata, base.fs, twindow=base.opt_twin, threshold=base.opt_th)
                    dict["vlf"]  = copt.vlf
                    dict["lf"] = copt.lf
                    dict["vlar"] = copt.vlar
                    dict["rsam"] = copt.rsam
                    dict["lrar"] = copt.lrar
                    dict["mf"] = copt.mf
                    dict["rmar"] = copt.rmar
                    dict["hf"] = copt.hf
                end
            end
        end
        
        if base.polar
            data_n = @view base.data[:, ninikk:ninikk+base.lwin]
            polar = _polar(data_n, base.fs, base.fqminmax, base.fq_band, NW=base.NW, pd=base.pad)
            dict["polar"]["degree"][n,:]  = polar.degree
            dict["polar"]["rect"][n,:]    = polar.rect
            dict["polar"]["azimuth"][n,:] = polar.azimuth
            dict["polar"]["elev"][n,:]    = polar.elev
        end
    
    end
    
    return dlte
end