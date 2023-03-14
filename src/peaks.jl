#!/usr/local/bin julia
# coding=utf-8

# Utility functions for LTE.jl

# Ivan Melchor

using Peaks


function extract(spec::Array{T,3}, degree::Array{T}, rect::Array{T}, azimuth::Array{T}, elev::Array{T}, freq::Array{T}, peak_thresholds::Dict) where T<:Real

    nro_comp, total_time, nro_freq  = size(spec)

    pks_out  = [Dict() for i in 1:nro_comp]
    nro_Wpks = 0
    nro_Lpks = 0

    for t in 1:total_time
        cond != any(isnan.(spec[:,t,:])) # if cont is true, proceed

        if cond
            deg_t  = degree[t,:]
            
            for c in 1:nro_comp
                sxx_ct   = spec[c,t,:]
                sxx_scal = _scaler(sxx_ct)
                pks_loc, pksxx = findmaxima(sxx_scal)
                peakheights!(pks_loc, pksxx, minheight=peak_thresholds["sxx_th"])

                if any(pks)
                    pks_out[c][t] = Dict("fq"=>[], "specgram"=>[], "degree"=>[], "rect"=>[], "azimuth"=>[], "elev"=>[])

                    for pk in pks_loc
                        fq_pk = freq[pk]
                        fq_min, fq_max = (fq_pk - peak_thresholds["fq_dt"], fq_pk + peak_thresholds["fq_dt"])
                        fleft  = findfirst(x -> x >= fq_min, freq)
                        fright = findfirst(fq_max .<= freq)

                        pk_sxx_avg = mean(sxx_ct[fleft:fright])
                        
                        deg_t_pk = deg_t[fleft:fright]
                        pk_pd_avg  = mean(deg_t_pk)
                        pk_pd_std  = std(deg_t_pk, corrected=false, mean=pk_pd_avg)

                        if pk_pd_avg > peak_thresholds["pdg_th"] && pd_peak_std < peak_thresholds["pdg_std"]

                            # save dominant peak
                            nro_Wpks += 1
                            push!(pks_out[c][t]["fq"], fq_pk)
                            push!(pks_out[c][t]["specgram"], pk_sxx_avg)
                            push!(pks_out[c][t]["degree"], pk_pd_avg)

                            # get rectilinearity measure
                            rect_t_pk = rect[t,fleft:fright]
                            deg_cond = deg_t_pk .> peak_thresholds["pdg_th"]
                            rect_t_pk = rect_t_pk[deg_cond]

                            if any(rect_t_pk)
                                pk_rect_avg  = mean(rect_t_pk)
                                pk_rect_std  = std(rect_t_pk, corrected=false, mean=pk_rect_avg)

                                if pk_rect_avg > peak_thresholds["rect_th"] && pk_rect_std <  peak_thresholds["rect_std"]

                                    # this peak is linearly polarized
                                    nro_Lpks += 1
                                    push!(pks_out[c][t]["rect"], pk_rect_avg)

                                    # get azimuth measures
                                    azim_t_pk = azimuth[t,fleft:fright]
                                    azim_t_pk = azim_t_pk[deg_cond]
                                    rect_cond = rect_t_pk .> peak_thresholds["rect_th"]
                                    azim_t_pk = azim_t_pk[rect_cond]

                                    # get elevation measures
                                    elev_t_pk = elev[t,fleft:fright]
                                    elev_t_pk = elev_t_pk[deg_cond]
                                    elev_t_pk = elev_t_pk[rect_cond]

                                    if any(azim_t_pk)
                                        pk_azim_avg = _circmean(pi/180 .* azim_t_pk) * 180/pi
                                        pk_azim_std = _circstd(pi/180 .* azim_t_pk)
                                        pk_elev_avg = mean(elev_t_pk)
                                        pk_elev_std = std(elev_t_pk, corrected=false, mean=pk_elev_avg)

                                        if pk_azim_std < peak_thresholds["azim_std"]
                                            push!(pks_out[c][t]["azimuth"], pk_azim_avg)
                                        else
                                            push!(pks_out[c][t]["azimuth"], nothing)
                                        end

                                        if pk_elev_std < peak_thresholds["elev_std"]
                                            push!(pks_out[c][t]["elev"], pk_elev_avg)
                                        else
                                            push!(pks_out[c][t]["elev"], nothing)
                                        end
                                    end

                                else
                                    push!(pks_out[c][t]["rect"], nothing)
                                    push!(pks_out[c][t]["azimuth"], nothing)
                                    push!(pks_out[c][t]["elev"], nothing)
                                
                                end
                            end
                        end
                    end
                end
            end 
        end
    end

    return pks_out, (nro_Wpks, nro_Lpks)
end

