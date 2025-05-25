flag = @isdefined NOISEINTERPOLATORSHEADERGUARD

if !flag
    # const vel_interp = [0,3,6,9,12,15,18,21,24,27]
    const dbHz_interp = [30+3*i for i in 1:6]
    const vel0 = linear_interpolation(dbHz_interp, [0.26,0.2,0.14,0.1,0.08,0.05],extrapolation_bc=Interpolations.Line())
    const vel3 = linear_interpolation(dbHz_interp, [1.05,0.71,0.41,0.27,0.19,0.1],extrapolation_bc=Interpolations.Line())
    const vel6 = linear_interpolation(dbHz_interp, [1.09,0.77,0.54,0.35,0.22,0.11],extrapolation_bc=Interpolations.Line())
    const vel9 = linear_interpolation(dbHz_interp, [1.46,0.77,0.56,0.43,0.25,0.1],extrapolation_bc=Interpolations.Line())
    const vel12 = linear_interpolation(dbHz_interp, [1.62,1.06,0.82,0.49,0.25,0.11],extrapolation_bc=Interpolations.Line())
    const vel15 = linear_interpolation(dbHz_interp, [1.18,1.06,0.85,0.53,0.27,0.11],extrapolation_bc=Interpolations.Line())
    const vel18 = linear_interpolation(dbHz_interp, [1.34,0.77,0.56,0.42,0.24,0.09],extrapolation_bc=Interpolations.Line())
    const vel21 = linear_interpolation(dbHz_interp, [1.26,1.14,0.75,0.35,0.25,0.1],extrapolation_bc=Interpolations.Line())
    const vel24 = linear_interpolation(dbHz_interp, [0.81,0.58,0.44,0.31,0.25,0.08],extrapolation_bc=Interpolations.Line())
    const vel27 = linear_interpolation(dbHz_interp, [0.9,0.61,0.51,0.3,0.31,0.09],extrapolation_bc=Interpolations.Line())
    doppler_interpolators = [vel0,vel3,vel6,vel9,vel12,vel15,vel18,vel21,vel24,vel27]
end


const NOISEINTERPOLATORSHEADERGUARD = 1

function measnoise!(R, recstate, svposs, svvels, cnr, valid)
    numsvs = div(length(cnr), 2)

    recpos = ECEF(recstate[1:3])
    velocity = min(norm(recstate[4:6]), 27)
    
    las = lookangles.(svposs, Ref(recpos))
    el = [la[2] for la in las]

    # Params:
    T = 50
    F = 10
    A = 30
    a = 20
    L2BETTERRATIO = 1

    cnr = rnxssi2dbhz(cnr)

    s = cnr[1:numsvs]
    R[1:numsvs, 1:numsvs] = diagm(noiseformula(s, el, T, F, A, a))
    s = cnr[numsvs+1:2*numsvs]
    R[numsvs+1:2*numsvs, numsvs+1:2*numsvs] = diagm(noiseformula(s, el, T, F, A, a)) / L2BETTERRATIO

    s = cnr
    for i in 1:2*numsvs
        R[i + 2*numsvs,i + 2*numsvs] = doppler_interpolators[floor(Int,velocity/3)+1](s[i])
        if s[i] < 30
            R[i + 2*numsvs,i + 2*numsvs] *= 3
        elseif s[i] > 48
            R[i + 2*numsvs,i + 2*numsvs] = 0.09
        end
    end

    R[isnan.(R)] .= 1e6
    R[:,:] .*= 100
    return nothing
end


function rnxssi2dbhz(rnxssi)
    return 4 .+ 6 .* rnxssi
end

function noiseformula(cnr, el, T, F, A, a)
    r = 1 .+ 1 ./ ((sin.(el)).^2) .* (10 .^ (.-(cnr.-T)./a) .* (((A)./(10 .^(.-(F.-T)./a)) .-1) .* (cnr.-T)./(F.-T) .+ 1))
    r[cnr .> T] .= 1.
    return r
end


#CREATE YOUR OWN NOISE MODEL.... THE COPIED ONE IS NOT WORKING
# seems that using SSI is not a good idea...


function measnoise_sim!(R, recstate, svposs, svvels, cnr, valid)
    numsvs = length(cnr)

    recpos = ECEF(recstate[1:3])
    velocity = min(norm(recstate[4:6]), 27)
    
    az, el = enulookangles(svposs, recpos)

    # Params:
    T = 50
    F = 10
    A = 30
    a = 20
    L2BETTERRATIO = 1

    cnr = rnxssi2dbhz(cnr)

    s = cnr
    R[1:numsvs, 1:numsvs] = diagm(noiseformula(s, el, T, F, A, a))

    # doppler weights
    for i in 1:numsvs
        R[i + numsvs,i + numsvs] = doppler_interpolators[floor(Int,velocity/3)+1](s[i])
        if s[i] < 30
            R[i + numsvs,i + numsvs] *= 3
        elseif s[i] > 48
            R[i + numsvs,i + numsvs] = 0.09
        end
    end

    R[isnan.(R)] .= 1e6
    R[:,:] .*= 100
    return nothing
end