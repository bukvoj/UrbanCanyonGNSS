using Geodesy, GeoStats
using Serialization, DataFrames
using NearestNeighbors
using Statistics, LinearAlgebra, Distributions
using Base.Threads, Distributed
using RinexRead, GNSSEphemeris, GNSSMultipathSim
using Dates, TimesDates

# RINEX files
NAV_FILE = "data/2024_10_7/rinex/line22_fromhostivartopohorelec.24N"
OBS_FILE = "data/2024_10_7/rinex/line22_fromhostivartopohorelec.24O"

# map files
OSM_TRACKMAP_FILE = "data/2024_10_7/osm/geojson/prague_tramtracks.geojson"
OSM_ROUTE_XML = "data/2024_10_7/osm/route_xml/lin22_hostivar2bilahora.xml"
OSM_BUILDINGS_DIR = "data/2024_10_7/osm/geojson"
OSM_BUILDINGS_PREFIX = "buildings_"

# Output file
OUTPUTFILE = "simulateddata/multipathdata.dat"


# Load the navigation and observation data from the RINEX files
nav = rinexread(NAV_FILE)
obs = rinexread(OBS_FILE)
println("Observations loaded.")


if !@isdefined(UrbanCanyonGNSS)
    include("../src/UrbanCanyonGNSS.jl")
    using .UrbanCanyonGNSS
else
    println("UrbanCanyonGNSS already loaded.")
end


# Load the map data
trackmap = loadmap(OSM_TRACKMAP_FILE, OSM_ROUTE_XML)
trackmap = resamplemap(trackmap; step = 5)
walls = wallsfromdirectory(OSM_BUILDINGS_DIR, trackmap.center; filestartswith = OSM_BUILDINGS_PREFIX)
println("3D map loaded.")


# Generate ground truth positions
flag = @isdefined results_already_generated
if !flag 
    results_already_generated = true
    
    observations_parsed = init_batch_processing(obs, nav; prevresults = nothing, useklobuchar=true, usetropospheric=true)
    results = trajectory(observations_parsed;       # Created by init_batch_processing() or the simulation script
               runmpestimation = false,    # Run multipath estimation algorithm
               elanglelookup = nothing,    # Elevation angle lookup table
               envmap=nothing,           # 3D environment map (meshgrid + track map)
               noisemodel = nothing,       # Function that modifies the covariance matrix - has the form: measnoise!(R, x, svpos, svvel, ssi, valid)
               kf = nothing,               # Filter to use (default is IteratedExtendedKalmanFilter with CVM model initialized at 0)
               maskangle = 15,             # Mask angle for the filter
               mpwinlen = 10,              # Length of the moving window for the multipath estimation
               ssimask = 0,                # Mask for the SSI values
               PFA = 0.01,                 # False alarm probability for the multipath estimation
               σ = 1                       # Standard deviation of healthy pseudorange measurements (for the multipath estimation)
        )

end


# Project positions to the map.... will be used as ground truth positions
df = DataFrame()
ecef = ECEF.([x[1:3] for x in results.x])
ecef = unique(ecef)
df.projected = proj2map.(Ref(trackmap), ecef)                                               # USES results from trajectory.jl
x = [x[1] for x in df.projected]
y = [x[2] for x in df.projected]
df.RecPos = [GeoStats.Point{🌐}(Cartesian3D{WGS84Latest}(x[i], y[i], 0.0)) for i in eachindex(x)] # convert to Point objects
df.Time = DateTime.(unique(results.time))
println("Ground truth positions calculated.")

# Calculate SV positions and velocities
obs = obs.data
numsvs = renumber!(nav, obs)
numgps = length(unique(nav.data.GPS.SatelliteID))
numbeidou = length(unique(nav.data.BEIDOU.SatelliteID))
numgalileo = length(unique(nav.data.GALILEO.SatelliteID))
svposses = []
svvels = []
println("Number of SVs: $numsvs")
for t in df.Time
    tmp = []
    for g in 1:numgps
        svpos = getsvpos(t, g, 'G', nav)
        push!(tmp, svpos)
    end
    for b in 1+numgps:numbeidou+numgps
        svpos = getsvpos(t, b, 'C', nav)
        push!(tmp, svpos)
    end
    for e in 1+numgps+numbeidou:numgps+numbeidou+numgalileo
        svpos = getsvpos(t, e, 'E', nav)
        push!(tmp, svpos)
    end
    tmppos = [trackmap.ecef2enu(x[1]) for x in tmp]
    tmpvel = [trackmap.ecef2enu(x[2]) for x in tmp]
    push!(svposses, tmppos)
    push!(svvels, tmpvel)
end
df.SvPos = svposses
df.SvVel = svvels
println("SV positions and velocities calculated.")
df.Id = 1:length(df.Time)


# START OFF PARALLEL SECTION
println("SPAWNING WORKERS")
addprocs(8)  # or julia -p n

@everywhere using DataFrames, Geodesy, GeoStats, GeoIO, LinearAlgebra, NearestNeighbors, Statistics
@everywhere using LinearAlgebra
@everywhere using GNSSMultipathSim

println("SPAWNING WORKERS DONE")
@everywhere walls = @fetchfrom 1 walls
println("SPAWNING WALLS DONE")

@everywhere function modifytherow(df)
    t = df.Time
    i = df.Id
    println("Processing measurement $i")
    sv = df.SvPos
    global walls
    sv = [GeoStats.Point{🌐}(Cartesian3D{WGS84Latest}(s[1],s[2],s[3])) for s in sv]
    mp, mpmode = getmpmeasurements(df.RecPos, sv, walls)
    pure = norm.(sv.-df.RecPos)
    return (Id = i, Time = t, SvPos = df.SvPos, RecPos = df.RecPos, L1 = mp, MPMode = mpmode, GTMeas = pure, SvVel = df.SvVel)
end



rows = eachrow(df) |> collect  # convert to vector of rows
row_data = [NamedTuple(r) for r in rows]
modified_rows = pmap(modifytherow, row_data) # Run the computation in parallel
new_df = DataFrame(modified_rows)
df = copy(new_df)
# END OF PARALLEL SECTION

open(OUTPUTFILE, "w") do io
    serialize(io, df)
end

for w in workers()
    println("KILLING WORKER $w")
    rmprocs(w)
end

