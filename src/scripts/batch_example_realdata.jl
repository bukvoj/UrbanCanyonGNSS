
using GLMakie #  # for plotting


using DataFrames, Dates, TimesDates, Geodesy

using LowLevelParticleFilters

using RinexRead

# Input data files
navfile = "data/2024_10_7/rinex/line22_fromhostivartopohorelec.24N"
obsfile = "data/2024_10_7/rinex/line22_fromhostivartopohorelec.24O"

include("../src/utils.jl")
include("../src/visualisation/geoplot.jl")

include("../src/trajectoryestimation/initialization.jl")
include("../src/trajectoryestimation/trajectory.jl")

include("../src/motionmodels/cvm.jl")
include("../src/measurementmodels/measmodels.jl")
include("../src/measurementmodels/noise_model.jl") # Realini model



# Load the RINEX files
nav = rinexread(navfile)
obs = rinexread(obsfile)
println("RINEX files loaded.")

# Prepare the observations for processing
observations_parsed = init_batch_processing(obs, nav; prevresults = nothing, useklobuchar=true, usetropospheric=true)
println("Observations parsed.")

# Compute the trajectory (plain IEKF without any augmentations)
println("Computing results without multipath estimation...")
results = trajectory(observations_parsed;       # Created by init_batch_processing() or the simulation script
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
observations_parsed = init_batch_processing(obs, nav; prevresults = results, useklobuchar=true, usetropospheric=true)
results_withmp = trajectory(observations_parsed;       # Created by init_batch_processing() or the simulation script
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
pos = [ECEF(state[1:3]...) for state in results.x]
posmp = [ECEF(state[1:3]...) for state in results_withmp.x]

# Convert ECEF to LLA
llapos = LLA.(pos, Ref(wgs84))
lat = [p.lat for p in llapos]
lon = [p.lon for p in llapos]
# Convert ECEF to LLA with multipath estimation
llaposmp = LLA.(posmp, Ref(wgs84))
latmp = [p.lat for p in llaposmp]
lonmp = [p.lon for p in llaposmp]

# plot the results
println("Plotting results...")
m, route, fig, ax = geoplot(lat, lon; color = :blue, style = :scatter, size = (900,900), markersize = 5, label = "Trajectory")
m2, route2, fig2, ax2 = geoplot!(m, latmp, lonmp; color = :red, style = :scatter, markersize = 5, label = "Trajectory with multipath estimation")
# if it disappears, call translate!(route, 0, 0, 10) to move it back into the window... Tyler is weird like that

