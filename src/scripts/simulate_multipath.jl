using Geodesy
using Serialization, DataFrames, XML
using NearestNeighbors
using Statistics, LinearAlgebra, Distributions
using Base.Threads, Distributed
using RinexRead, GNSSEphemeris, GNSSMultipathSim
using Dates, TimesDates
using LowLevelParticleFilters

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

include("../src/trajectoryestimation/initialization.jl")
include("../src/dataloaders/loadtrackmap.jl")
include("../src/dataloaders/loadbuildings.jl")

function ecef2enurotmat(lat, lon)
    lat = lat * π / 180.0
    lon = lon * π / 180.0
    return [
        -sin(lon) cos(lon) 0.0;
        -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat);
        cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)
    ]
end


trackmap = loadmap(OSM_TRACKMAP_FILE, OSM_ROUTE_XML)
trackmap = resamplemap(trackmap; step = 5)
walls = wallsfromdirectory(OSM_BUILDINGS_DIR, trackmap.center; filestartswith = OSM_BUILDINGS_PREFIX)
println("3D map loaded.")



flag = @isdefined results_already_generated
if !flag 
    results_already_generated = true
    
    # TODO
    TODO: get results from trajectory.jl

end


# Project positions to the map.... will be used as ground truth positions
df = DataFrame()
ecef = ECEF.([x[1:3] for x in results.x])
ecef = unique(ecef)
df.projected = proj2map.(Ref(trackmap), ecef)                                               # USES results from trajectory.jl
x = [x[1] for x in df.projected]
y = [x[2] for x in df.projected]
df.recpos = [gs.Point{🌐}(Cartesian3D{WGS84Latest}(x[i], y[i], 0.0)) for i in eachindex(x)] # convert to Point objects
println("Ground truth positions calculated.")



# Calculate SV positions and velocities in ENU coordinates
obs = obs.data
numsvs = renumber!(nav, obs)
numgps = length(unique(nav.data.GPS.SatelliteID))
numbeidou = length(unique(nav.data.BEIDOU.SatelliteID))
numgalileo = length(unique(nav.data.GALILEO.SatelliteID))
svposses = []
svvels = []
println("Number of SVs: $numsvs")
for t in df.time
    tmp = []
    for g in 1:numgps
        svpos = getsvpos(t + Date(2024,7,10), g, 'G', nav)
        push!(tmp, svpos)
    end
    for b in 1+numgps:numbeidou+numgps
        svpos = getsvpos(t + Date(2024,7,10), b, 'C', nav)
        push!(tmp, svpos)
    end
    for e in 1+numgps+numbeidou:numgps+numbeidou+numgalileo
        svpos = getsvpos(t + Date(2024,7,10), e, 'E', nav)
        push!(tmp, svpos)
    end
    tmppos = [trackmap.ecef2enu(x[1]) for x in tmp]
    tmpvel = [trackmap.ecef2enu(x[2]) for x in tmp]
    push!(svposses, tmppos)
    push!(svvels, tmpvel)
end
df.svpos = svposses
df.svvel = svvels
println("SV positions and velocities calculated.")

# Initialize the measurement data as empty
df.pure .= Ref([])
df.mp .= Ref([])
df.mpmode .= Ref([])
df.id = 1:length(df.time)

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
    t = df.time
    i = df.id
    println("Processing measurement $i")
    sv = df.svpos
    global walls
    sv = [GeoStats.Point{🌐}(Cartesian3D{WGS84Latest}(s[1],s[2],s[3])) for s in sv]
    mp, mpmode = getmpmeasurements(df.recpos, sv, walls)
    pure = norm.(sv.-df.recpos)
    return (id = i, time = t, svpos = sv, recpos = df.recpos, mpmeas = mp, mpmode = mpmode, gtmeas = pure, svvel = df.svvel)
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

