flag = @isdefined INTERPOLATORSHEADERGUARD

if !flag
    using Interpolations
    const lats_interp = [15,30,45,60,75]
    const P₀s = linear_interpolation(lats_interp, [1013.25,1017.25,1015.75,1011.75,1013.00],extrapolation_bc=Interpolations.Line())
    const T₀s = linear_interpolation(lats_interp, [299.65 ,294.15 ,283.15 ,272.15 ,263.65],extrapolation_bc=Interpolations.Line())
    const e₀s = linear_interpolation(lats_interp, [26.31  ,21.79  ,11.66  ,6.78   ,4.11],extrapolation_bc=Interpolations.Line())
    const β₀s = linear_interpolation(lats_interp, [6.3e-3 ,6.05e-3,5.58e-3,5.39e-3,4.53e-3],extrapolation_bc=Interpolations.Line())
    const λ₀s = linear_interpolation(lats_interp, [2.77   ,3.15   ,2.57   ,1.81   ,1.55],extrapolation_bc=Interpolations.Line())

    const ΔPs = linear_interpolation(lats_interp, [0.,-3.75,-2.25,-1.75,-0.50],extrapolation_bc=Interpolations.Line())
    const ΔTs = linear_interpolation(lats_interp, [0.,7.00 ,11.00,15.00,14.50],extrapolation_bc=Interpolations.Line())
    const Δes = linear_interpolation(lats_interp, [0.,8.85 ,7.24 ,5.36 ,3.39],extrapolation_bc=Interpolations.Line())
    const Δβs = linear_interpolation(lats_interp, [0.,0.25e-3,0.32e-3,0.81e-3,0.62e-3],extrapolation_bc=Interpolations.Line())
    const Δλs = linear_interpolation(lats_interp, [0.,0.33,0.46,0.74,0.30],extrapolation_bc=Interpolations.Line())
end


const INTERPOLATORSHEADERGUARD = 1


function tropomodel(lla, az, el, time)
    k₁ = 77.604
    k₂ = 382000
    Rₐ = 287.054
    gₘ = 9.784
    g = 9.80665

    lat = lla.lat
    h = lla.alt
    if lat < 0
        Dmin = 211
    else 
        Dmin = 28
    end
    D = day(time)
    Δt = cos(2*pi*(D-Dmin)/365.25)

    # Look it up bro. It will work, bro
    P = P₀s(lat) + ΔPs(lat)*Δt
    T = T₀s(lat) + ΔTs(lat)*Δt
    e = e₀s(lat) + Δes(lat)*Δt
    β = β₀s(lat) + Δβs(lat)*Δt
    λ = λ₀s(lat) + Δλs(lat)*Δt

    T_z0_dry = 1e-6 * P  * Rₐ * k₁ / gₘ
    T_z0_wet = 1e-6 * k₂ * Rₐ * e / (T*((λ+1)*gₘ - β * Rₐ ))

    T_dry = (1-β*h/T)^(g/Rₐ*β) * T_z0_dry
    T_wet = (1-β*h/T)^((λ + 1)/Rₐ*β - 1) * T_z0_wet

    M = 1.001/sqrt(0.002001 + sin(el)^2)
    return M*(T_dry + T_wet)
end