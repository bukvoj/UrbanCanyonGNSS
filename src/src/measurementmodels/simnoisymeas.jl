"""
genmeas(mpmeas, numsvs = 28)
mpmeas is a DataFrame obtained from the simulation of the multipath
this function adds noise and biases to the measurements
and returns a new DataFrame with the modified measurements
and the biases

This is to test the algorithm with different noise and biases
"""
function genmeas(mpmeas)
    numsvs = length(mpmeas.svpos[1])
    c = 299792458

    # Copy mpmeaservations into Geodesy objects...
    noisy_mpmeas = DataFrame()
    noisy_mpmeas.id = mpmeas.id
    noisy_mpmeas.time = mpmeas.time
    noisy_mpmeas.svpos .= Ref([ECEF(0.,0.,0.) for i in 1:numsvs])
    noisy_mpmeas.svvel .= Ref([ECEF(0.,0.,0.) for i in 1:numsvs])
    noisy_mpmeas.recpos .= [ECEF(x.coords.x.val, x.coords.y.val, x.coords.z.val) for x in mpmeas.recpos]
    noisy_mpmeas.rho .= Ref(zeros(numsvs))
    # convert to Geodesy objects
    for i in eachindex(noisy_mpmeas.svpos)
        noisy_mpmeas.svpos[i] = [ECEF(x.coords.x.val, x.coords.y.val, x.coords.z.val) for x in mpmeas.svpos[i]]
        noisy_mpmeas.rho[i] = [x.val for x in mpmeas.mp[i]]
        # noisy_mpmeas.rho[i] = [x.val for x in mpmeas.pure[i]]
    end
    noisy_mpmeas.recvel = vcat([noisy_mpmeas.recpos[2] - noisy_mpmeas.recpos[1]], diff(noisy_mpmeas.recpos))

    # Generate measurements with noise and biases
    dnormal = Distributions.Normal(0, 1)
    dmultipath = Distributions.Normal(0, 3)
    dbiases = Uniform(-10,10)
    biases = rand(dbiases, numsvs)
    t_offset0 = c * 1e-5
    t_offsetrate = c*1e-6
    noisy_mpmeas.doppler .= Ref(zeros(numsvs))
    noisy_mpmeas.pure = mpmeas.pure
    for i in eachindex(noisy_mpmeas.rho)
        p = (noisy_mpmeas.svpos[i], noisy_mpmeas.svvel[i], ones(numsvs*2))
        x = [noisy_mpmeas.recpos[i]..., noisy_mpmeas.recvel[i]..., t_offset0 + i*t_offsetrate, t_offsetrate, zeros(numsvs)...]
        noisy_mpmeas.rho[i] .+= rand(dnormal, numsvs) + (mpmeas.mpmode[i] .> 1).*rand(dmultipath, numsvs) + biases .+ t_offset0 .+ i*t_offsetrate
        noisy_mpmeas.doppler[i] =  meas_sim(x,nothing,p,i)[numsvs+1:end] + 0.1 .* rand(dnormal, numsvs) + 0.05 * (mpmeas.mpmode[i] .> 1).*rand(dmultipath, numsvs)
    end
    return noisy_mpmeas, biases, t_offset0, t_offsetrate
end