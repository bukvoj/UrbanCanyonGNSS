function rho_clk_correction(rho::Real, id::Int, time::TimeDate, navdata::DataFrame, svpos::ECEF,svvel::ECEF)
    time = timedate2unix(time)
    return rho_clk_correction(rho, id, time, navdata, svpos, svvel)
end



function rho_clk_correction(rho::Real, id::Int, time::DateTime, navdata::DataFrame, svpos::ECEF,svvel::ECEF)
    time = datetime2unix(time)
    return rho_clk_correction(rho, id, time, navdata, svpos, svvel)
end

function rho_clk_correction(rho::Real, id::Int, time::Real, navdata::DataFrame, svpos::ECEF,svvel::ECEF)
    c = 299792458
    nav = navdata[navdata.SatelliteID .== id, :]
    nav = nav[argmin(abs.(timedate2unix.(nav.Time) .- time)),:]
    navtime = timedate2unix(nav.Time)

    Δ_rel = -2*svpos'*svvel/(c^2) #+delta Trel2 ?
    bias = nav.SVClockBias + (time-navtime)*nav.SVClockDrift + 1/2*nav.SVClockDriftRate*(time-navtime).^2
    rho = rho + c*(bias + Δ_rel)
    return rho
end

function rho_tgd_correction(rho::Real, id::Int, time::TimeDate, navdata::DataFrame, constellation::Char; freq=false)
    # instrumentation delay correction
    c = 299792458
    nav = navdata[navdata.SatelliteID .== id, :]
    nav = nav[argmin(abs.(nav.Time .- time)),:]
    
    if constellation == 'G'
        bias = nav.TGD 
        if freq
            bias *=  0.512
        end
    elseif constellation == 'C' && !freq
        bias = nav.TGD1
    elseif constellation == 'C' && freq
        bias = nav.TGD2
    elseif constellation == 'E'
        bias = nav.BGDE5bE1
    end 

    rho = rho + c*bias
    return rho
end



function secondofweek(t::DateTime)
    return 3600*hour(t) + 60*minute(t) + second(t) + 24*3600*(Dates.dayofweek(t) % 7)
end

function secondofweek(t::TimeDate)
    return 3600*hour(t) + 60*minute(t) + second(t) + 24*3600*(Dates.dayofweek(t) % 7)
end