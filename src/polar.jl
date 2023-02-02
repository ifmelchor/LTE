#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

using ComplexValues

"""
    _polar(*args)

Perform polarization analysis in the frequency domain.
For increasing the confidence in the measure use lwin, nadv, and nwin and twoside (optional) that applies a moving-window approach
"""
function _polar(data::Array{Float64,2}, base::LTEBase)

    # compute de the SVD of the spectral covariance matrix of ZNE data
    freq, csm_svd = csm(data, base.fs, base.lswin, base.nswin, base.nadv, base.fqminmax, base.nfs, base.NW, base.pad)

    # compute polar degree
    degree = [(3*sum(s.S.^2)-sum(s.S)^2)/(2*sum(s.S)^2) for s in csm_svd]
    
    # compute rectilinearity of Vt
    z_rot = map(_rotate_vector, [s.Vt[1,:] for s in csm_svd])
    rect  = map(_rectiliniarity, z_rot)
    
    # get angles
    angles = map(_angles, z_rot)
    azi = [a.azimuth for a in angles]
    ele = [a.elev for a in angles]
    pHH = [a.phyHH for a in angles]
    pVH = [a.phyVH for a in angles]

    return PParams(freq, degree, rect, azi, ele, pHH, pVH)
end


function _rotate_vector(z::Array{ComplexF64})
    
    zrot_(a) = [zk * (cos(a) + sin(a)im) for zk in z]
    alla = range(0, 2*pi, 100)
    allz = [zrot_(a) for a in alla]
    comp = [abs(dot([x.re for x in z_rotk], [x.im for x in z_rotk])) for z_rotk in allz]
    phi_min = alla[findmin(comp)[2]]
    
    return zrot_(phi_min)
end


function _rectiliniarity(z::Array{ComplexF64})
    
    # z = _rotate_vector(z)
    a = norm([x.re for x in z])
    b = norm([x.im for x in z])
    sminor, smajor = sort([a,b])
    
    return 1 - (sminor/smajor)
end


function _elevation(z::Array{ComplexF64})
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


function _angles(z::Array{ComplexF64})
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
    zN_rot = z[2] * exp(-th_H*1im)
    zE_rot = z[3] * exp(-th_H*1im)
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
    pHH = (zEpol.ang - zNpol.ang)*180/pi

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

