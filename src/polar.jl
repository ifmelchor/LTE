#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

"""
    _polar(*args)

Functions for computing polarization analysis in the frequency domain

"""
function _polar(data::AbstractArray{T}, base::Base; return_all::Bool=false, full_return::Bool=true) where T<:Real

    # compute de the SVD of the spectral covariance matrix of ZNE data
    csm_svd = _csm(data, base.fs, base.lswin, base.nswin, base.nadv, base.fqminmax, base.NW, base.pad, return_all=return_all)

    if return_all && !isnothing(base.nswin)
        frmin, frmax = base.fqminmax
        nfs = frmax-frmin+1
        degree = Array{Float32}(undef, base.nswin, nfs)

        if full_return
            rect = Array{Float32}(undef, base.nswin, nfs)
            azi  = Array{Float32}(undef, base.nswin, nfs)
            ele  = Array{Float32}(undef, base.nswin, nfs)
            pHH  = Array{Float32}(undef, base.nswin, nfs)
            pVH  = Array{Float32}(undef, base.nswin, nfs)
        end

        for n in 1:base.nswin
            degree[n,:] = [(3*sum(s.S.^2)-sum(s.S)^2)/(2*sum(s.S)^2) for s in csm_svd[n]]

            if full_return
                z_rot = map(_rotate_vector, [s.Vt[1,:] for s in csm_svd[n]])
                rect[n,:] = map(_rectiliniarity, z_rot)
                angles = map(_angles, z_rot)
                azi[n,:] = [a.azimuth for a in angles]
                ele[n,:] = [a.elev for a in angles]
                pHH[n,:] = [a.phyhh for a in angles]
                pVH[n,:] = [a.phyvh for a in angles]
                to_return = PParams(degree, rect, azi, ele, pHH, pVH) 
            else
                to_return = degree
            end
        end

    else
        # compute polar degree
        degree = [(3*sum(s.S.^2)-sum(s.S)^2)/(2*sum(s.S)^2) for s in csm_svd]
        
        if full_return
            # compute rectilinearity of Vt
            z_rot = map(_rotate_vector, [s.Vt[1,:] for s in csm_svd])
            rect  = map(_rectiliniarity, z_rot)
            
            # get angles
            angles = map(_angles, z_rot)
            azi = [a.azimuth for a in angles]
            ele = [a.elev for a in angles]
            pHH = [a.phyhh for a in angles]
            pVH = [a.phyvh for a in angles]

            to_return = PParams(degree, rect, azi, ele, pHH, pVH)    
        else
            to_return = degree
        end
    end

    return to_return

end


function _rotate_vector(z::Array{C}) where C<:Complex
    
    zrot_(a) = [zk * (cos(a) + sin(a)im) for zk in z]
    alla = range(0, 2*pi, 100)
    allz = [zrot_(a) for a in alla]
    comp = [abs(dot([x.re for x in z_rotk], [x.im for x in z_rotk])) for z_rotk in allz]
    phi_min = alla[findmin(comp)[2]]
    
    return zrot_(phi_min)
end


function _rectiliniarity(z::Array{C}) where C<:Complex
    
    # z = _rotate_vector(z)
    a = norm([x.re for x in z])
    b = norm([x.im for x in z])
    sminor, smajor = sort([a,b])
    
    return 1 - (sminor/smajor)
end


function _elevation(z::Array{C}) where C<:Complex
    # compute elevation (0 for vertical, 90 for horizontal)
    zVpol = Polar(z[1])
    zH = sqrt(z[2]^2 + z[3]^2)
    zHpol = Polar(zH)
    zHVpol = Polar(z[1]^2 + zH^2)

    th(m) = -0.5*zHVpol.ang + m*pi/2
    func(th) = (zVpol.mod*cos(th + zVpol.ang))^2 + (zHpol.mod*cos(th + zHpol.ang))^2
    th_V = th(findmax(map(func, [th(l) for l in 0:5]))[2]-1)

    if zH.im < 0
        zH = Complex(-1*zH.re, -1*zH.im)
    end

    ztmp = Complex(cos(th_V), -1*sin(th_V))
    tV = pi/2 - atan(abs( (z[1]*ztmp).re / (zH*ztmp).re))

    return tV*180/pi
end


function _angles(z::Array{C}) where C<:Complex
    # transform complex vectors
    zVpol = Polar(z[1])
    zNpol = Polar(z[2])
    zEpol = Polar(z[3])

    # decompose horizontal comp. and minimize func
    zHpol = Polar(z[2]^2 + z[3]^2)
    th(l) = -0.5*zHpol.ang + l*pi/2
    func(th) = (zNpol.mod*cos(th+zNpol.ang))^2 + (zEpol.mod*cos(th+zEpol.ang))^2
    th_H = th(findmax(map(func, [th(l) for l in 0:5]))[2]-1)

    # rotate zH and compute tH
    zN_rot = z[2] * cis(-th_H) #exp(-th_H*1im)
    zE_rot = z[3] * cis(-th_H) #exp(-th_H*1im)
    tH = atan(zE_rot.re/zN_rot.re)

    arg = z[1]*conj(z[3])
    if arg.re < 0
        if tH < 0
            tH += pi
        end
    else
        if tH > 0
            tH -= pi
        end
    end

    # measure clockwise from N
    tH = pi/2 - tH
    if tH < 0
        tH += 2*pi
    end

    # compute phiHH, which is the phase difference between horizontals
    pHH = (zEpol.ang - zNpol.ang) * 180/pi

    if pHH > 180
        pHH -= 360
    end

    if pHH < -180
        pHH += 360
    end

    # compute pVH, which is the phase difference between vertical and horizontal
    pVH = (th_H - zVpol.ang)*180/pi

    if pVH > 90
        pVH -= 180
    end

    if pVH < -90
        pVH += 180
    end

    # we only seek for the direction (no orientation)
    # if tH > pi:
    #     tH -= pi
    
    return Angles(tH*180/pi, _elevation(z), pHH, pVH)
end


function time_polar(u::Array{T,2}) where T<:Real

   # compute covariance matrix
   C = Array{T}(undef, (3, 3))
   for i in 1:3, j in 1:3
    C[i,j] = cov(u[i,:],u[j,:])
   end

   # compute eigenanalysis
   # and sort into ascending order 
   f  = eigen(C, sortby=+)
   d  = f.values # eigenvalues
   v3 = f.vectors[:,3] #eigenvector for the first eigenvalue

   # v3[1] is Z component
   # v3[2] is N component
   # v3[3] is E component

   # calculate the rectilinearity
   rL = 1 - (d[1]+d[2])/(2*d[3])

   # calculate the azimuth clockwise from N
   baz = atan(v3[3], v3[2]) * (180/pi)
   if baz < 0
    baz += 360
   end

   # 180 ambiguity
   if baz > 180
      baz -= 180
   end

   # calculate incidence angle. 0 is vertical, 90 horizontal
   inc = atan(abs( sqrt(v3[2]^2 + v3[3]^2)/v3[1] )) * (180/pi)


   return rL, baz, inc
end


function time_polar(u::Array{T,2}, win::J, olap::T=0.0) where {T<:Real, J<:Integer}

   npts = size(u, 2)
   nolp   = round(Int, win*olap)
   nadv = win - nolp
   N    = floor(Int, (npts - 1 - nolp) / nadv)
   no = 1
   rL  = Array{T}(undef, N)
   baz = Array{T}(undef, N)
   inc = Array{T}(undef, N)
   
   for i in 1:N
      un_pol = time_polar(u[:, no:no+win])
      rL[i]  = un_pol[1]
      baz[i] = un_pol[2]
      inc[i] = un_pol[3]
      no += nadv
   end

   return rL, baz, inc
end
