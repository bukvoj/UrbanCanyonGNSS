function init_batch_processing(obs, nav; prevresults = nothing, useklobuchar=true, usetropospheric=true)
    
    obs = deepcopy(obs.data)
    nav = deepcopy(nav)

    ionoparams = nav.header.ionocorrections.GPS
    ionoalpha = [ionoparams.a0,ionoparams.a1, ionoparams.a2, ionoparams.a3]
    ionobeta = [ionoparams.b0,ionoparams.b1, ionoparams.b2, ionoparams.b3]

    # Remove SVs without ephemeris data
    rm_missing!(obs.GPS, nav.data.GPS)
    rm_missing!(obs.GALILEO, nav.data.GALILEO)
    rm_missing!(obs.BEIDOU, nav.data.BEIDOU)

    # Renumber SVs in obs, nav data to have them numbered from 1 to numsvs
    numsvs = renumber!(nav, obs)

    # Get SV positions and velocities and add them to the obs data
    svpos!(obs, nav; results=prevresults) 

    # Corrections
    applyclock!(obs,nav) # Apply satellite clock corrections
    if useklobuchar
        applyklobuchar!(obs, nav, ionoalpha, ionobeta)
    end
    if usetropospheric
        applytropo!(obs, nav) # Apply tropospheric corrections
    end

    # normalize doppler measurements
    normalizedoppler!(obs)

    # Connect the data into single dataframe
    times = vcat(obs.GPS.Time, obs.BEIDOU.Time, obs.GALILEO.Time)
    l1 = vcat(obs.GPS.C1C, obs.BEIDOU.C1I, obs.GALILEO.C1X)
    l2 = vcat(obs.GPS.C2X, obs.BEIDOU.C7I, obs.GALILEO.C7X)

    d1 = vcat(obs.GPS.D1C, obs.BEIDOU.D1I, obs.GALILEO.D1X)
    d2 = vcat(obs.GPS.D2X, obs.BEIDOU.D7I, obs.GALILEO.D7X)
    ssi1 = vcat(obs.GPS.C1C_SSI, obs.BEIDOU.C1I_SSI, obs.GALILEO.C1X_SSI)
    ssi2 = vcat(obs.GPS.C2X_SSI, obs.BEIDOU.C7I_SSI, obs.GALILEO.C7X_SSI)

    svpos = vcat(obs.GPS.SvPos, obs.BEIDOU.SvPos, obs.GALILEO.SvPos)
    svvel = vcat(obs.GPS.SvVel, obs.BEIDOU.SvVel, obs.GALILEO.SvVel)

    # Create the dataframe
    df = DataFrame(
        Time = times,
        SatelliteID = vcat(obs.GPS.SatelliteID, obs.BEIDOU.SatelliteID, obs.GALILEO.SatelliteID),
        L1 = l1,
        L2 = l2,
        D1 = d1,
        D2 = d2,
        L1_SSI = ssi1,
        L2_SSI = ssi2,
        SvPos = svpos,
        SvVel = svvel
    )
    # Sort using time
    sort!(df, [:Time], rev=false)
    return df
end



function rm_missing!(obsdata,navdata)
    ids = unique(obsdata.SatelliteID)
    filter!(:SatelliteID => in(navdata.SatelliteID), obsdata)
end

function svpos!(obs, nav; results = nothing)
    c = 299792458

    gpsoffsets = zeros(length(obs.GPS.SatelliteID)) * Nanosecond(1)
    beidoffsets = zeros(length(obs.BEIDOU.SatelliteID)) * Nanosecond(1)
    galoffets = zeros(length(obs.GALILEO.SatelliteID)) * Nanosecond(1)

    if !isnothing(results)
        function findoffset(t, df, prevend)
            for r in prevend:length(df.time)
                if df.time[r] == t
                    return Nanosecond(floor(df.x[r][7] / c * 1e9)), r
                end
            end
            println("no match for $t")
            return Nanosecond(0)
        end
        prevend = 1
        for i in eachindex(obs.GPS.Time)
            t = obs.GPS.Time[i]
            gpsoffsets[i],prevend = findoffset(t, results,prevend)
        end
        prevend = 1
        for i in eachindex(obs.BEIDOU.Time)
            t = obs.BEIDOU.Time[i]
            beidoffsets[i],prevend = findoffset(t, results,prevend)
        end
        prevend = 1
        for i in eachindex(obs.GALILEO.Time)
            t = obs.GALILEO.Time[i]
            galoffets[i],prevend = findoffset(t, results, prevend)
        end
    end

    ## GPS
    rho = copy(obs.GPS.C1C)
    rho[isnan.(rho)] .= copy(obs.GPS.C2X[isnan.(rho)])
    svs = getsvpos.(obs.GPS.Time + gpsoffsets, obs.GPS.SatelliteID, 'G', Ref(nav), rho)
    obs.GPS.SvPos = [i[1] for i in svs]
    obs.GPS.SvVel = [i[2] for i in svs]
    ## BEIDOU
    rho = copy(obs.BEIDOU.C1I)
    rho[isnan.(rho)] .= copy(obs.BEIDOU.C7I[isnan.(rho)])
    svs = getsvpos.(obs.BEIDOU.Time + beidoffsets, obs.BEIDOU.SatelliteID, 'C', Ref(nav), rho)
    obs.BEIDOU.SvPos = [i[1] for i in svs]
    obs.BEIDOU.SvVel = [i[2] for i in svs]
    ## GALILEO
    rho = copy(obs.GALILEO.C1X)
    rho[isnan.(rho)] .= copy(obs.GALILEO.C7X[isnan.(rho)])
    svs = getsvpos.(obs.GALILEO.Time + galoffets, obs.GALILEO.SatelliteID, 'E', Ref(nav), rho)
    obs.GALILEO.SvPos = [i[1] for i in svs]
    obs.GALILEO.SvVel = [i[2] for i in svs]
end


function applyclock!(obs, nav)
    obs.GPS.C1C = rho_clk_correction.(obs.GPS.C1C, obs.GPS.SatelliteID, obs.GPS.Time, Ref(nav.data.GPS), obs.GPS.SvPos, obs.GPS.SvVel)
    obs.GPS.C2X = rho_clk_correction.(obs.GPS.C2X, obs.GPS.SatelliteID, obs.GPS.Time, Ref(nav.data.GPS), obs.GPS.SvPos, obs.GPS.SvVel)
    obs.BEIDOU.C1I = rho_clk_correction.(obs.BEIDOU.C1I, obs.BEIDOU.SatelliteID, obs.BEIDOU.Time, Ref(nav.data.BEIDOU), obs.BEIDOU.SvPos, obs.BEIDOU.SvVel)
    obs.BEIDOU.C7I = rho_clk_correction.(obs.BEIDOU.C7I, obs.BEIDOU.SatelliteID, obs.BEIDOU.Time, Ref(nav.data.BEIDOU), obs.BEIDOU.SvPos, obs.BEIDOU.SvVel)
    obs.GALILEO.C1X = rho_clk_correction.(obs.GALILEO.C1X, obs.GALILEO.SatelliteID, obs.GALILEO.Time, Ref(nav.data.GALILEO), obs.GALILEO.SvPos, obs.GALILEO.SvVel)
    obs.GALILEO.C7X = rho_clk_correction.(obs.GALILEO.C7X, obs.GALILEO.SatelliteID, obs.GALILEO.Time, Ref(nav.data.GALILEO), obs.GALILEO.SvPos, obs.GALILEO.SvVel)

    obs.GPS.C1C = rho_tgd_correction.(obs.GPS.C1C, obs.GPS.SatelliteID, obs.GPS.Time, Ref(nav.data.GPS), 'G')
    obs.GPS.C2X = rho_tgd_correction.(obs.GPS.C2X, obs.GPS.SatelliteID, obs.GPS.Time, Ref(nav.data.GPS), 'G'; freq = true)
    obs.BEIDOU.C1I = rho_tgd_correction.(obs.BEIDOU.C1I, obs.BEIDOU.SatelliteID, obs.BEIDOU.Time, Ref(nav.data.BEIDOU), 'C')
    obs.BEIDOU.C7I = rho_tgd_correction.(obs.BEIDOU.C7I, obs.BEIDOU.SatelliteID, obs.BEIDOU.Time, Ref(nav.data.BEIDOU), 'C'; freq=true)
    obs.GALILEO.C1X = rho_tgd_correction.(obs.GALILEO.C1X, obs.GALILEO.SatelliteID, obs.GALILEO.Time, Ref(nav.data.GALILEO), 'E')
    obs.GALILEO.C7X = rho_tgd_correction.(obs.GALILEO.C7X, obs.GALILEO.SatelliteID, obs.GALILEO.Time, Ref(nav.data.GALILEO), 'E'; freq=true)
end

function applyklobuchar!(obs, nav, alpha, beta)
    # Constants:
    c = 299792458
    gpsf1 = 1575.42e6
    gpsf2 = 1227.60e6
    galf1 = 1575.42e6
    galf2 = 1207.14e6
    beif1 = 1561.098e6
    beif2 = 1207.14e6

    lla = LLA(50.09, 14.42,200)
    ## GPS
    las = lookangles.(obs.GPS.SvPos, Ref(lla))
    az = [la[1] for la in las]
    el = [la[2] for la in las]
    time = secondofweek.(obs.GPS.Time)
    ionodelays = klobuchar.(Ref(lla),az,el,time, Ref(alpha), Ref(beta))
    obs.GPS.C1C -= c * ionodelays
    obs.GPS.C2X -= c * ionodelays * (gpsf1/gpsf2)^2
    obs.GPS.ionodelays = ionodelays
    ## BEIDOU
    las = lookangles.(obs.BEIDOU.SvPos, Ref(lla))
    az = [la[1] for la in las]
    el = [la[2] for la in las]
    time = secondofweek.(obs.BEIDOU.Time)
    ionodelays = klobuchar.(Ref(lla),az,el,time, Ref(alpha), Ref(beta))
    obs.BEIDOU.C1I -= c * ionodelays
    obs.BEIDOU.C7I -= c * ionodelays * (beif1/beif2)^2
    obs.BEIDOU.ionodelays = ionodelays
    ## GALILEO
    las = lookangles.(obs.GALILEO.SvPos, Ref(lla))
    az = [la[1] for la in las]
    el = [la[2] for la in las]
    time = secondofweek.(obs.GALILEO.Time)
    ionodelays = klobuchar.(Ref(lla),az,el,time, Ref(alpha), Ref(beta))
    obs.GALILEO.C1X -= c * ionodelays
    obs.GALILEO.C7X -= c * ionodelays * (galf1/galf2)^2
    obs.GALILEO.ionodelays = ionodelays
end

function applytropo!(obs,nav)
    lla = LLA(50.09, 14.42,200)
    ## GPS
    las = lookangles.(obs.GPS.SvPos, Ref(lla))
    az = [la[1] for la in las]
    el = [la[2] for la in las]
    tropodelays = tropomodel.(Ref(lla),az,el,obs.GPS.Time)
    obs.GPS.C1C -= tropodelays
    obs.GPS.C2X -= tropodelays
    obs.GPS.tropodelays = tropodelays
    ## BEIDOU
    las = lookangles.(obs.BEIDOU.SvPos, Ref(lla))
    az = [la[1] for la in las]
    el = [la[2] for la in las]
    tropodelays = tropomodel.(Ref(lla),az,el,obs.BEIDOU.Time)
    obs.BEIDOU.C1I -= tropodelays
    obs.BEIDOU.C7I -= tropodelays
    obs.BEIDOU.tropodelays = tropodelays
    ## GALILEO
    las = lookangles.(obs.GALILEO.SvPos, Ref(lla))
    az = [la[1] for la in las]
    el = [la[2] for la in las]
    tropodelays = tropomodel.(Ref(lla),az,el,obs.GALILEO.Time)
    obs.GALILEO.C1X -= tropodelays
    obs.GALILEO.C7X -= tropodelays
    obs.GALILEO.tropodelays = tropodelays
end

function normalizedoppler!(obs)
    # Constants:
    c = 299792458
    gpsf1 = 1575.42e6
    gpsf2 = 1227.60e6
    galf1 = 1575.42e6
    galf2 = 1207.14e6
    beif1 = 1561.098e6
    beif2 = 1207.14e6

    obs.GPS.D1C *= c / gpsf1
    obs.GPS.D2X *= c / gpsf2
    obs.BEIDOU.D1I *= c / beif1
    obs.BEIDOU.D7I *= c / beif2
    obs.GALILEO.D1X *= c / galf1
    obs.GALILEO.D7X *= c / galf2
end


function renumber!(nav, obs)
    # Renumber SVs in obs, nav data
    obs.GPS.SatelliteID .+= 100
    obs.BEIDOU.SatelliteID .+= 200
    obs.GALILEO.SatelliteID .+= 300
    nav.data.GPS.SatelliteID .+= 100
    nav.data.BEIDOU.SatelliteID .+= 200
    nav.data.GALILEO.SatelliteID .+= 300

    gpsids = unique(nav.data.GPS.SatelliteID)
    galids = unique(nav.data.GALILEO.SatelliteID)
    beids = unique(nav.data.BEIDOU.SatelliteID)
    navsatids = vcat(gpsids,beids,galids)
    numsvs = length(navsatids)
    for i in 1:numsvs
        obs.GPS.SatelliteID[obs.GPS.SatelliteID .== navsatids[i]] .= i
        obs.BEIDOU.SatelliteID[obs.BEIDOU.SatelliteID .== navsatids[i]] .= i
        obs.GALILEO.SatelliteID[obs.GALILEO.SatelliteID .== navsatids[i]] .= i
        nav.data.GPS.SatelliteID[nav.data.GPS.SatelliteID .== navsatids[i]] .= i
        nav.data.BEIDOU.SatelliteID[nav.data.BEIDOU.SatelliteID .== navsatids[i]] .= i
        nav.data.GALILEO.SatelliteID[nav.data.GALILEO.SatelliteID .== navsatids[i]] .= i
    end
    return numsvs
end