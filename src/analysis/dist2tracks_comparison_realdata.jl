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
navfile = "data/2024_10_7/rinex/line22_fromhostivartopohorelec.24N"
obsfile = "data/2024_10_7/rinex/line22_fromhostivartopohorelec.24O"

# elangle lookup table
elanglelookupfile = "simulateddata/mplookuptable.dat"

# trackmap file
OSM_TRACKMAP_FILE = "data/2024_10_7/osm/geojson/prague_tramtracks.geojson"
OSM_ROUTE_XML = "data/2024_10_7/osm/route_xml/lin22_hostivar2bilahora.xml"
trackmap = loadmap(OSM_TRACKMAP_FILE, OSM_ROUTE_XML)



# Load the precomputed environment maps...
try 
    global mplookup = deserialize(elanglelookupfile)
catch e
    println("Error loading the elevation angle lookup table: ", e)
    println("Creating a new lookup table: This might take hours!")
    include("../scripts/create_elangle_lookup_table.jl")
    global mplookup = deserialize(elanglelookupfile)
end
println("Elevation angle lookup table loaded.")



# Load the RINEX files
nav = rinexread(navfile)
obs = rinexread(obsfile)
println("RINEX files loaded.")

# Prepare the observations for processing
observations_parsed = init_batch_processing(obs, nav; prevresults = nothing, useklobuchar=true, usetropospheric=true)
println("Observations parsed.")

plainresults = trajectory(observations_parsed)
observations_parsed = init_batch_processing(obs, nav; prevresults = plainresults, useklobuchar=true, usetropospheric=true)
println("Plain results computed. Starting testing different algorithms...")

println("realini + no multipath estimation")
realini_iekf = trajectory(observations_parsed; runmpestimation = false, elanglelookup = nothing, envmap = nothing, noisemodel = measnoise!, kf = nothing, maskangle = 15, mpwinlen = 10, ssimask = 0, PFA = 0.01, σ = 1)
println("realini + multipath estimation")
realini_mp = trajectory(observations_parsed; runmpestimation = true, elanglelookup = nothing, envmap = nothing, noisemodel = measnoise!, kf = nothing, maskangle = 15, mpwinlen = 10, ssimask = 0, PFA = 0.01, σ = 1)
println("realini + multipath estimation + map partitioning")

println("realini + multipath estimation + elangle lookup")
realini_mpmap = trajectory(observations_parsed; runmpestimation = true, elanglelookup = mplookup, envmap = nothing, noisemodel = measnoise!, kf = nothing, maskangle = 15, mpwinlen = 10, ssimask = 0, PFA = 0.01, σ = 1)


println("uniform + no multipath estimation")
uniform_iekf = trajectory(observations_parsed; runmpestimation = false, elanglelookup = nothing, envmap = nothing, noisemodel = nothing, kf = nothing, maskangle = 15, mpwinlen = 10, ssimask = 0, PFA = 0.01, σ = 1)
println("uniform + multipath estimation")
uniform_mp = trajectory(observations_parsed; runmpestimation = true, elanglelookup = nothing, envmap = nothing, noisemodel = nothing, kf = nothing, maskangle = 15, mpwinlen = 10, ssimask = 0, PFA = 0.01, σ = 1)
println("uniform + multipath estimation + map partitioning")

println("uniform + multipath estimation + elangle lookup")
uniform_mpmap = trajectory(observations_parsed; runmpestimation = true, elanglelookup = mplookup, envmap = nothing, noisemodel = nothing, kf = nothing, maskangle = 15, mpwinlen = 10, ssimask = 0, PFA = 0.01, σ = 1)

println("DONE computing!")







function results2dist2map(results, trackmap)
    ecef = ECEF.([x[1:3] for x in results.x])
    # Calculate distances to map
    distances = dist2map.(Ref(trackmap), ecef)
    return distances
end


uniform_iekf_distances = results2dist2map(uniform_iekf, trackmap)
uniform_mp_distances = results2dist2map(uniform_mp, trackmap)
uniform_mpmap_distances = results2dist2map(uniform_mpmap, trackmap)
# uniform_chunk_distances = results2dist2map(uniform_rays, trackmap)

realini_iekf_distances = results2dist2map(realini_iekf, trackmap)
realini_mp_distances = results2dist2map(realini_mp, trackmap)
realini_mpmap_distances = results2dist2map(realini_mpmap, trackmap)
# realini_chunk_distances = results2dist2map(realini_rays, trackmap)

println("Distances calculated.")
println("Uniform iekf RMS dist: ", rms(uniform_iekf_distances))
println("Uniform mp RMS dist: ", rms(uniform_mp_distances))
println("Uniform mpmap RMS dist: ", rms(uniform_mpmap_distances))
# println("Uniform chunks RMS dist: ", rms(uniform_chunk_distances))

println("Realini iekf RMS dist: ", rms(realini_iekf_distances))
println("Realini mp RMS dist: ", rms(realini_mp_distances))
println("Realini mpmap RMS dist: ", rms(realini_mpmap_distances))
# println("Realini chunks RMS dist: ", rms(realini_chunk_distances))