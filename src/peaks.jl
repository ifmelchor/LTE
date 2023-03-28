#!/usr/local/bin julia
# coding=utf-8

# Utility functions for LTE.jl

# Ivan Melchor

using Peaks


function extract(spec::Array{T,3}, degree::Array{T}, rect::Array{T}, azimuth::Array{T}, elev::Array{T}, freq::Array{T}, pksth::PeakThresholds) where T<:Real

    nro_comp, total_time, nro_freq  = size(spec)

    pks_out  = [Dict() for i in 1:nro_comp]
    nro_Wpks = [0 for i in 1:nro_comp]
    nro_Lpks = [0 for i in 1:nro_comp]

    for t in 1:total_time
        cond = any(isnan.(spec[:,t,:])) # avoid extracting peaks if nans in data

        if !cond
            deg_t  = degree[t,:]
            
            for c in 1:nro_comp
                sxx_ct   = spec[c,t,:]
                sxx_scal = _scaler(sxx_ct)
                pks_loc, pksxx = findmaxima(sxx_scal)
                peakheights!(pks_loc, pksxx, minheight=pksth.sxx_th)

                if length(pks_loc) > 1
                    pks_out[c][t] = Dict("fq"=>[], "specgram"=>[], "degree"=>[], "rect"=>[], "azimuth"=>[], "elev"=>[])

                    for pk in pks_loc
                        # get frequency and power of the peak
                        fq_pk = freq[pk]
                        fq_min, fq_max = (fq_pk-pksth.fq_dt, fq_pk+pksth.fq_dt)
                        fleft  = findfirst(x -> x >= fq_min, freq)
                        fright = findfirst(fq_max .<= freq)

                        if isnothing(fleft)
                            fleft = 1
                        end

                        if isnothing(fright)
                            fright = length(freq)
                        end

                        pk_sxx_avg = mean(sxx_ct[fleft:fright])
                        
                        # compute degree of the peak
                        deg_t_pk   = deg_t[fleft:fright]
                        pk_pd_avg  = mean(deg_t_pk)
                        pk_pd_std  = std(deg_t_pk, corrected=false, mean=pk_pd_avg)

                        if pk_pd_avg > pksth.pdg_th && pk_pd_std < pksth.pdg_std

                            # save dominant peak
                            nro_Wpks[c] += 1
                            push!(pks_out[c][t]["fq"], fq_pk)
                            push!(pks_out[c][t]["specgram"], pk_sxx_avg)
                            push!(pks_out[c][t]["degree"], pk_pd_avg)

                            # get rectilinearity measure
                            rect_t_pk = rect[t,fleft:fright]
                            deg_cond = deg_t_pk .> pksth.pdg_th
                            rect_t_pk = rect_t_pk[deg_cond]

                            if length(rect_t_pk) >= pksth.npts_min
                                pk_rect_avg  = mean(rect_t_pk)
                                pk_rect_std  = std(rect_t_pk, corrected=false, mean=pk_rect_avg)

                                if pk_rect_std < pksth.rect_std
                                    push!(pks_out[c][t]["rect"], pk_rect_avg)
                                else
                                    push!(pks_out[c][t]["rect"], NaN)
                                end

                                if pk_rect_avg > pksth.rect_th && pk_rect_std < pksth.rect_std
                                    # this peak is linearly polarized
                                    nro_Lpks[c] += 1

                                    # get azimuth measures
                                    azim_t_pk = azimuth[t,fleft:fright]
                                    azim_t_pk = azim_t_pk[deg_cond]
                                    rect_cond = rect_t_pk .> pksth.rect_th
                                    azim_t_pk = azim_t_pk[rect_cond]

                                    # get elevation measures
                                    elev_t_pk = elev[t,fleft:fright]
                                    elev_t_pk = elev_t_pk[deg_cond]
                                    elev_t_pk = elev_t_pk[rect_cond]

                                    if length(azim_t_pk) >= pksth.npts_min
                                        pk_azim_avg = _circmean(azim_t_pk)
                                        pk_azim_std = _circstd(azim_t_pk)
                                        pk_elev_avg = mean(elev_t_pk)
                                        pk_elev_std = std(elev_t_pk, corrected=false, mean=pk_elev_avg)
                                        
                                        if pk_azim_std < pksth.azim_std
                                            push!(pks_out[c][t]["azimuth"], pk_azim_avg)
                                        else
                                            push!(pks_out[c][t]["azimuth"], NaN)
                                        end
                                        
                                        if pk_elev_std < pksth.elev_std
                                            push!(pks_out[c][t]["elev"], pk_elev_avg)
                                        else
                                            push!(pks_out[c][t]["elev"], NaN)
                                        end
                                    
                                    else
                                        push!(pks_out[c][t]["azimuth"], NaN)
                                        push!(pks_out[c][t]["elev"], NaN)
                                    end

                                else
                                    push!(pks_out[c][t]["azimuth"], NaN)
                                    push!(pks_out[c][t]["elev"], NaN)
                                end
                            else
                                push!(pks_out[c][t]["rect"], NaN)
                                push!(pks_out[c][t]["azimuth"], NaN)
                                push!(pks_out[c][t]["elev"], NaN)
                            end
                        end
                    end
                end
            end 
        end
    end

    # clean empty dicts
    for pkso in pks_out
        itoclean = []
        for di in pkso
            if length(di[2]["fq"]) < 1
                push!(itoclean, di[1])
            end
        end
        
        for i in itoclean
            delete!(pkso, i)
        end
    end

    return pks_out, (nro_Wpks, nro_Lpks)
end

