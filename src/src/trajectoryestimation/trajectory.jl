# DEPENDENCIES:
# TODO

function trajectory(nav::RinexContent,obs::RinexContent; 
    useklobuchar=true,          # Use Klobuchar model
    usetropospheric=true,       # Use tropospheric model
    runmpestimation = false,    # Run multipath estimation algorithm
    elanglelookup = nothing,    # Elevation angle lookup table
    envmap=nothing,             # 3D environment map (meshgrid + track map)
    noisemodel = nothing,       # Function that modifies the covariance matrix - has the form: measnoise!(R, x, svpos, svvel, ssi, valid)
    prevresults=nothing,        # Previous results to use for the filter (if trajectory is ran for the second time, use the previous results for better time estimation)
    kf = nothing,               # Filter to use (default is IteratedExtendedKalmanFilter with CVM model initialized at 0)
    maskangle = 15,             # Mask angle for the filter
    mpwinlen = 10,              # Length of the moving window for the multipath estimation
    ssimask = 0,                # Mask for the SSI values
    PFA = 0.01,                 # False alarm probability for the multipath estimation
    σ = 1,                      # Standard deviation of healthy pseudorange measurements (for the multipath estimation)
    )

    obs = deepcopy(obs.data)
    nav = deepcopy(nav)

    params = nav.header.ionocorrections.GPS
    alpha = [params.a0,params.a1, params.a2, params.a3]
    beta = [params.b0,params.b1, params.b2, params.b3]

    # Remove SVs without ephemeris data
    rm_missing!(obs.GPS, nav.data.GPS)
    rm_missing!(obs.GALILEO, nav.data.GALILEO)
    rm_missing!(obs.BEIDOU, nav.data.BEIDOU)

    # Renumber SVs in obs, nav data to have them numbered from 1 to numsvs
    numsvs = renumber!(nav, obs)

    # Get SV positions and velocities
    svpos!(obs, nav; results=prevresults)

    # Corrections
    applyclock!(obs,nav) # Apply satellite clock corrections
    if doklobuchar
        applyklobuchar!(obs, nav, alpha, beta)
    end
    if dotropo
        applytropo!(obs, nav) # Apply tropospheric corrections
    end

    # normalize doppler measurements
    normalizedoppler!(obs)

    if isnothing(kf)
        # Init filter
        x0 = zeros(8)
        d0 = MvNormal(x0, I(8)*1e10 + zeros(8,8))
        cvmQ = [1/3 1/2; 1/2 1]
        Q = [4*cvmQ zeros(2,6); 
                zeros(2,2) 4*cvmQ zeros(2,4);
                zeros(2,4) 4*cvmQ zeros(2,2);
                zeros(2,6) diagm([100,1])]    
        
        R = zeros(4*numsvs,4*numsvs) + I(4*numsvs) * σ^2
        kf = IteratedExtendedKalmanFilter(cvm,meas,Q,R,d0;nu=0,ny=numsvs, Cjac=meas_jac,Ajac=cvm_jac)
    end


    colnames =[
        "time",
        "x",
        "Σ",
        "ids",
        "freq",
        "res",
        "meastype",
        "az",
        "el",
        "r",
        "ssi",
        "svpos",
        "svvel",
        "valid",
        "mpstrat"
    ]
    out = DataFrame([name=>[] for name in colnames])
    i = 1

    # preallocate memory
    y = zeros(4 * numsvs)
    valid = zeros(4 * numsvs)
    svpos = [ECEF(0.,0.,0.) for i in 1:numsvs]
    svvel = [ECEF(0.,0.,0.) for i in 1:numsvs]
    ssi = zeros(numsvs*2)
    res = zeros(4*numsvs)
    u = nothing

    inovations = zeros(2*numsvs, mpwinlen)
    biases = zeros(2*numsvs)

    # Compute chi-square critical value for 1% false alarm
    if runmpestimation
        α =  quantile(Chisq(mpwinlen), 1-PFA)
        mp_filters = [MP_SV(0,0,σ^2,1,Circbuffer(mpwinlen)) for i in 1:2*numsvs]
    end
    mps = zeros(2*numsvs)

    println("Starting the main loop.")
    uniquetimes = unique(obs.GPS.Time)
    for t in uniquetimes
        fillmeas!(t, y, ssi, valid, svpos, svvel, obs, numsvs)

        az, el = lookangles(svpos, ECEF(iekf.x[1],iekf.x[2],iekf.x[3]); unit=:deg)

        mask = vcat(el .< maskangle,el .< maskangle,el .< maskangle,el .< maskangle)
        valid[mask] .= 0
        mask = vcat(ssi .< ssimask, ssi .< ssimask)
        valid[mask] .= 0

        p = (svpos, svvel, valid)
        
        # Run filter step
        time = timedate2unix(t)
        predict!(kf,u,p,time)
        
        if !isnothing(noisemodel)
            noisemodel!(kf.R2, kf.x, svpos, svvel, ssi, valid)
        end
        
        if (!isnothing(elanglelookup)) && (i > 1)
            usempmap!(kf.R2, elanglelookup, kf, az, el) # change
        end
        
        if (!isnothing(envmap)) && (i > 1)
            println(i)
            usechunkmap!(kf.R2, envmap, ECEF(kf.x[1],kf.x[2],kf.x[3]), trackmap, svpos, valid)
        end
        
        # Apply the IEKF correction step
        p = (svpos, svvel, valid, biases)
        correct!(kf,u,y,p,time)

        kf.R2 = zeros(4*numsvs,4*numsvs) + I(4*numsvs) * σ^2

        res = y - meas(iekf.x,u,p,time) 
        if runmpestimation & (i > 10)
            res[valid.==0] .= 0

            for j in 1:2*numsvs
                if valid[j] == 1
                    m,r,mstrat = mp_step!(mp_filters[j], res[j] - biases[j], α)
                    biases[j] = m
                    iekf.R2[j,j] = r
                    mps[j] = mstrat
                end
            end
        end

        # Store results
        res[valid.==0] .= 0

        for j in eachindex(res)
            push!(out,[
                t,
                kf.x,
                kf.R,
                (j - 1) % numsvs + 1,
                ((j > numsvs) & (j < 2*numsvs+1)) | (j > 3*numsvs),
                res[j],
                j > 2*numsvs,
                az[(j - 1) % numsvs + 1],
                el[(j - 1) % numsvs + 1],
                iekf.R2[j,j],
                ssi[(j - 1) % (2*numsvs) + 1],
                svpos[(j - 1) % numsvs + 1],
                svvel[(j - 1) % numsvs + 1],
                valid[j],
                mps[(j - 1) % (2*numsvs) + 1]
            ])
        end
        i += 1
    end

    println("Finished the main loop.")
    states = []
    covariances = []
    for t in unique(out.times)
        x = out.x[out.time .== t][1]
        Σ = out.Σ[out.time .== t][1]
        push!(states, x)
        push!(covariances, Σ)
    end
    return states,covariances,out, obs #out and obs are debug dataframes
end

function fillmeas!(t, y, ssi, valid, svpos, svvel, obs, numsvs)
    # Get current observation
    currgps = @view obs.GPS[obs.GPS.Time .== t, :]
    currgal = @view obs.GALILEO[obs.GALILEO.Time .== t, :]
    currbei = @view obs.BEIDOU[obs.BEIDOU.Time .== t, :]

    y .*= NaN
    valid .= 0.

    y[currgps.SatelliteID] = currgps.C1C
    y[currgps.SatelliteID .+ numsvs] = currgps.C2X
    y[currgps.SatelliteID .+ 2*numsvs] = currgps.D1C
    y[currgps.SatelliteID .+ 3*numsvs] = currgps.D2X
    y[currbei.SatelliteID] = currbei.C1I
    y[currbei.SatelliteID .+ numsvs] = currbei.C7I
    y[currbei.SatelliteID .+ 2*numsvs] = currbei.D1I
    y[currbei.SatelliteID .+ 3*numsvs] = currbei.D7I
    y[currgal.SatelliteID] = currgal.C1X
    y[currgal.SatelliteID .+ numsvs] = currgal.C7X
    y[currgal.SatelliteID .+ 2*numsvs] = currgal.D1X
    y[currgal.SatelliteID .+ 3*numsvs] = currgal.D7X
    
    ssi[currgps.SatelliteID] = currgps.C1C_SSI
    ssi[currgps.SatelliteID .+ numsvs] = currgps.C2X_SSI
    ssi[currbei.SatelliteID] = currbei.C1I_SSI
    ssi[currbei.SatelliteID .+ numsvs] = currbei.C7I_SSI
    ssi[currgal.SatelliteID] = currgal.C1X_SSI
    ssi[currgal.SatelliteID .+ numsvs] = currgal.C7X_SSI
    valid[.!isnan.(y)] .= 1
    y[isnan.(y)] .= 0    

    svpos[currgps.SatelliteID] = currgps.SvPos
    svvel[currgps.SatelliteID] = currgps.SvVel
    svpos[currbei.SatelliteID] = currbei.SvPos
    svvel[currbei.SatelliteID] = currbei.SvVel
    svpos[currgal.SatelliteID] = currgal.SvPos
    svvel[currgal.SatelliteID] = currgal.SvVel
    return nothing
end






function usechunkmap!(R, map3d, pos::ECEF, trackmap::TrackMap, svpos, valid)
    numsvs = floor(Int, size(R,1) / 4)

    chunksize = map3d.chunksize
    ychunks = map3d.ychunks

    pt = proj2map(trackmap, pos)

    x = floor(Int, (pt[1] - map3d.xmin) / chunksize)
    y = floor(Int, (pt[2] - map3d.ymin) / chunksize)
    id = x * ychunks + y
    c = map3d.map
    chunks = [c[id], c[id + 1], c[id - 1], c[id + ychunks], c[id - ychunks], c[id + ychunks + 1], c[id - ychunks - 1], c[id + ychunks - 1], c[id - ychunks + 1]]
    
    # Find index of closest point on map
    ptt = gs.Point{🌐}(Cartesian3D{WGS84Latest}(pt[1], pt[2], 0))

    for i in 1:numsvs
        if (valid[i] == 1) || (valid[i+numsvs] == 1)
            for chunk in chunks
                for wall in chunk
                    if Meshes.intersects(Ray(ptt, gs.Point(trackmap.ecef2enu(svpos[i])...)-ptt), wall)
                        println("Chunkmap: ", i, " ", wall)
                        R[i,i] *= 100
                        R[i+numsvs,i+numsvs] *= 100
                        break
                    end
                end
            end
        end
    end

    
end