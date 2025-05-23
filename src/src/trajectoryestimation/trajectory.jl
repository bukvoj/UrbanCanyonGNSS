# DEPENDENCIES:
# TODO
using Geodesy, LowLevelParticleFilters, LinearAlgebra, Dates, Statistics, Distributions, TimesDates


include("../maps/trackmap.jl")
include("../maps/map3d.jl")
include("../maps/dist2map.jl")

include("../mpmitigation/mpestimation.jl")

function trajectory(data::DataFrame;       # Created by init_batch_processing() or the simulation script
    runmpestimation = false,    # Run multipath estimation algorithm
    elanglelookup = nothing,    # Elevation angle lookup table
    envmap=nothing,             # 3D environment map (meshgrid + track map)
    noisemodel = nothing,       # Function that modifies the covariance matrix - has the form: measnoise!(R, x, svpos, svvel, ssi, valid)
    kf = nothing,               # Filter to use (default is IteratedExtendedKalmanFilter with CVM model initialized at 0)
    maskangle = 15,             # Mask angle for the filter
    mpwinlen = 10,              # Length of the moving window for the multipath estimation
    ssimask = 0,                # Mask for the SSI values
    PFA = 0.01,                 # False alarm probability for the multipath estimation
    σ = 1,                      # Standard deviation of healthy pseudorange measurements (for the multipath estimation)
    )

    usel2 = ("l2" in names(data)) || ("L2" in names(data))
    numsvs = length(unique(data.SatelliteID))

    # Initialize the filter
    if isnothing(kf)
        # Init filter
        x0 = zeros(8)
        d0 = MvNormal(x0, I(8)*1e10 + zeros(8,8))
        cvmQ = [1/3 1/2; 1/2 1]
        Q = [4*cvmQ zeros(2,6); 
                zeros(2,2) 4*cvmQ zeros(2,4);
                zeros(2,4) 4*cvmQ zeros(2,2);
                zeros(2,6) diagm([100,1])]    
        
        R = zeros((1 + usel2) * 2*numsvs,(1 + usel2) * 2*numsvs) + I((1 + usel2) * 2*numsvs) * σ^2
        kf = IteratedExtendedKalmanFilter(cvm,meas,Q,R,d0;nu=0,ny=(1 + usel2)*2*numsvs, Cjac=meas_jac,Ajac=cvm_jac)
    end

    # Create the DF to store the results
    colnames =[
        "time",
        "x",
        "Σ",
        "ids",
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
    y = zeros((1 + usel2) * 2 * numsvs)
    valid = zeros((1 + usel2) * 2 * numsvs)
    svpos = [ECEF(0.,0.,0.) for i in 1:numsvs]
    svvel = [ECEF(0.,0.,0.) for i in 1:numsvs]
    ssi = zeros((1+usel2)*numsvs)
    res = zeros((1+usel2)*2*numsvs)
    u = nothing

    # # Compute chi-square critical value for PFA false alarm
    if runmpestimation
        α =  quantile(Chisq(mpwinlen), 1-PFA)
        mp_filters = [MP_SV(0,0,σ^2,1,Circbuffer(mpwinlen)) for i in 1:((1+usel2)*numsvs)]
    end
    mps = zeros((1+usel2)* numsvs)                  # Multipath mitigation strategy
    residuals = zeros((1+usel2)* numsvs, mpwinlen)  # Window for the multipath estimation
    B = zeros((1+usel2)* numsvs)                    # Estimated biases

    println("Starting the main loop.")
    uniquetimes = unique(data.Time)
    for t in uniquetimes
        unixtime = timedate2unix(t)
        fillmeas!(t, y, ssi, valid, svpos, svvel, data, numsvs)
        az, el = lookangles(svpos, ECEF(kf.x[1],kf.x[2],kf.x[3]); unit=:deg)
        if usel2
            valid[vcat(el .< maskangle,el .< maskangle,el .< maskangle,el .< maskangle)] .= 0
            valid[vcat(ssi .< ssimask, ssi .< ssimask)] .= 0
        else
            valid[vcat(el .< maskangle,el .< maskangle)] .= 0
            valid[ssi .< ssimask] .= 0
        end

        # Parameters for the filter
        p = (svpos, svvel, valid, B)
        
        # Run filter step
        predict!(kf,u,p,unixtime)

        # Select the IEKF weights for the measurements 
        if !isnothing(noisemodel)
            noisemodel!(kf.R2, kf.x, svpos, svvel, ssi, valid)
        else
            kf.R2 = zeros((1 + usel2) * 2 * numsvs, (1 + usel2) * 2 * numsvs) + I((1 + usel2) * 2*numsvs) * σ^2
        end
        
        # NLOS removal using the elevation angle lookup table
        if (!isnothing(elanglelookup)) && (i > 1)
            usempmap!(kf.R2, elanglelookup, kf, az, el) # change
        end
        
        # NLOS removal using the environment map (map partitioning in the thesis...)
        if (!isnothing(envmap)) && (i > 1)
            println(i)
            usechunkmap!(kf.R2, envmap.map3d, ECEF(kf.x[1],kf.x[2],kf.x[3]), envmap.trackmap, svpos, valid)
        end
        
        # Apply the IEKF correction step
        correct!(kf,u,y,p,unixtime)

        # Apply the multipath estimation algorithm
        res = y - meas(kf.x,u,p,unixtime) 
        if runmpestimation & (i > 10) # skip first 10 iterations ... we need to have good enough a priori estimates
            for j in 1:(1+usel2)*numsvs
                if valid[j] == 1
                    m,r,mstrat = mp_step!(mp_filters[j], res[j] - B[j], α)
                    B[j] = m
                    kf.R2[j,j] = r
                    mps[j] = mstrat
                end
            end
        end

        # Store results
        res[valid.==0] .= 0

        if usel2
            push!(out,[
                t,
                kf.x,
                kf.R,
                [1:numsvs; 1:numsvs; 1:numsvs; 1:numsvs],
                res,
                [[:L1 for i in 1:numsvs]; [:L2 for i in 1:numsvs]; [:D1 for i in 1:numsvs]; [:D2 for i in 1:numsvs]],
                vcat(az, az, az, az),
                vcat(el, el, el, el),
                kf.R2,
                vcat(ssi, ssi, ssi, ssi),
                vcat(svpos, svpos, svpos, svpos),
                vcat(svvel, svvel, svvel, svvel),
                valid,
                mps,
            ])
        else
            push!(out,[
                t,
                kf.x,
                kf.R,
                [1:numsvs; 1:numsvs],
                res,
                [[:L1 for i in 1:numsvs]; [:D1 for i in 1:numsvs]],
                vcat(az, az),
                vcat(el, el),
                kf.R2,
                vcat(ssi, ssi),
                vcat(svpos, svpos),
                vcat(svvel, svvel),
                valid,
                mps,
            ])
        end
        i += 1
    end

    println("Finished the main loop.")
    return out
end

function fillmeas!(t, y, ssi, valid, svpos, svvel, obs, numsvs)
    usel2 = ("l2" in names(obs)) || ("L2" in names(obs))
    # Get current observation
    curr = @view obs[obs.Time .== t, :]

    y .*= NaN
    valid .= 0.

    if usel2
        y[curr.SatelliteID] = curr.L1
        y[curr.SatelliteID .+ numsvs] = curr.L2
        y[curr.SatelliteID .+ 2*numsvs] = curr.D1
        y[curr.SatelliteID .+ 3*numsvs] = curr.D2

        ssi[curr.SatelliteID] = curr.L1_SSI
        ssi[curr.SatelliteID .+ numsvs] = curr.L2_SSI

    else
        y[curr.SatelliteID] = curr.L1
        y[curr.SatelliteID .+ numsvs] = curr.D1

        ssi[curr.SatelliteID] = curr.L1_SSI
    end
        
    valid[.!isnan.(y)] .= 1
    y[isnan.(y)] .= 0    

    svpos[curr.SatelliteID] = curr.SvPos
    svvel[curr.SatelliteID] = curr.SvVel
    return nothing
end






function usechunkmap!(R, map3d, pos::ECEF, trackmap::TrackMap, svpos, valid)
    numsvs = length(svpos)
    usel2 = (length(valid) == 4 * numsvs)

    chunksize = map3d.chunksize
    ychunks = map3d.ychunks

    # Convert to ENU and project to the map
    pt = proj2map(trackmap, pos)

    # Load the chunks around the receiver
    x = floor(Int, (pt[1] - map3d.xmin) / chunksize)
    y = floor(Int, (pt[2] - map3d.ymin) / chunksize)
    id = x * ychunks + y
    c = map3d.map
    chunks = [c[id], c[id + 1], c[id - 1], c[id + ychunks], c[id - ychunks], c[id + ychunks + 1], c[id - ychunks - 1], c[id + ychunks - 1], c[id - ychunks + 1]]
    
    # Convert to GeoStats point
    ptt = GeoStats.Point{🌐}(Cartesian3D{WGS84Latest}(pt[1], pt[2], 0))

    for i in 1:numsvs
        if (valid[i] == 1) || (valid[i+usel2*numsvs] == 1)
            chunkid = 1
            while chunkid <= 9
                chunk = chunks[chunkid]
                for wall in chunk
                    if Meshes.intersects(Ray(ptt, gs.Point(trackmap.ecef2enu(svpos[i])...)-ptt), wall)
                        println("Chunkmap: ", i, " ", wall)
                        R[i,i] *= 100
                        R[i+usel2*numsvs,i+usel2*numsvs] *= 100
                        chunkid = 10 # Ray blocked, break the loop
                        break
                    end
                end
            end
        end
    end    
end


