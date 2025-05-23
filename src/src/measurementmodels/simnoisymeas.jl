"""
genmeas(mpmeas, numsvs = 28)
mpmeas is a DataFrame obtained from the simulation of the multipath
this function adds noise and biases to the measurements
and returns a new DataFrame with the modified measurements
and the biases

This is to test the algorithm with different noise and biases
"""
function genmeas(mpmeas)
    numsvs = length(mpmeas.SvPos[1])
    c = 299792458
    
    # Copy mpmeaservations into Geodesy objects...
    colnames = ["Time", "SatelliteID", "RecPos", "RecVel", "SvPos", "SvVel", "L1", "D1", "GTMeas", "MPMode", "L1_SSI"]
    out = DataFrame([name=>[] for name in colnames])
    
    numsvs = length(mpmeas.SvPos[1])
    dnormal = Distributions.Normal(0, 1)
    dmultipath = Distributions.Normal(0, 3)
    dbiases = Uniform(-10,10)
    biases = rand(dbiases, numsvs)
    t_offset0 = c * 1e-5
    t_offsetrate = c*1e-6


    svpos = mpmeas.SvPos
    svvel = mpmeas.SvVel

    recpos = [ECEF(x.coords.x.val, x.coords.y.val, x.coords.z.val) for x in mpmeas.RecPos]
    recvel = vcat([recpos[2] - recpos[1]], diff(recpos))    
    
    for i in eachindex(mpmeas.Time)
        p = (svpos[i], svvel[i], ones(numsvs*2), zeros(numsvs))
        x = [recpos[i]..., recvel[i]..., t_offset0 + i*t_offsetrate, t_offsetrate, zeros(numsvs)...]
        L1 = [l.val for l in mpmeas.L1[i][1]] # REMOVE THE [1] LATER
        L1 = L1 + rand(dnormal, numsvs) + (mpmeas.MPMode[i] .!= :los).*rand(dmultipath, numsvs) + biases .+ t_offset0 .+ i*t_offsetrate
        D1 =  meas(x,nothing,p,i)[numsvs+1:end] + 0.1 .* rand(dnormal, numsvs) + 0.05 * (mpmeas.MPMode[i] .!= :los).*rand(dmultipath, numsvs)
        for j in 1:numsvs
            if mpmeas.MPMode[i][j] != :blocked
                push!(out, [TimeDate(mpmeas.Time[i]),j,recpos[i],recvel[i],svpos[i][j],svvel[i][j],L1[j],D1[j],mpmeas.GTMeas[i][j],mpmeas.MPMode[i][j],3])
            end 
        end
    end

    return out, biases, t_offset0, t_offsetrate
end