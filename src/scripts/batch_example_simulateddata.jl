
using GLMakie #  # for plotting
using DataFrames, Dates, TimesDates, Geodesy, Serialization

using RinexRead

# LOAD THE PACKAGE!
if !@isdefined(UrbanCanyonGNSS)
    include("../src/UrbanCanyonGNSS.jl")
    using .UrbanCanyonGNSS
else
    println("UrbanCanyonGNSS already loaded.")
end


# Input data files
simdatafile = "simulateddata/multipathdata.dat"


# Load the RINEX files
mpdata = deserialize(simdatafile)
println("Simulated multipath measurements loaded.")

# Prepare the observations for processing
noisymeasurements, biases, t_offset0, t_offsetrate = genmeas(mpdata)
println("Observations parsed.")

# Compute the trajectory (plain IEKF without any augmentations)
println("Computing results without multipath estimation...")
results = trajectory(noisymeasurements;       # Created by init_batch_processing() or the simulation script
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

# Compute results with the multipath estimation
println("Computing results with multipath estimation...")
results_withmp = trajectory(noisymeasurements;       # Created by init_batch_processing() or the simulation script
           runmpestimation = true,    # Run multipath estimation algorithm
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



# extract positions
pos = [ENU(state[1:3]...) for state in results.x]
posmp = [ENU(state[1:3]...) for state in results_withmp.x]

x = [x[1] for x in pos]
y = [x[2] for x in pos]
xmp = [x[1] for x in posmp]
ymp = [x[2] for x in posmp]


println("Plotting results...")
f,ax = plot(x, y, color = :blue, label = "Noisy measurements")
plot!(ax,xmp, ymp, color = :red, label = "Multipath estimation")

