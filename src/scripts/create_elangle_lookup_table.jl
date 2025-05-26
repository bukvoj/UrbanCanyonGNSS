# Create a lookup table for elevation angles based on the walls and waypoints in the map

OSM_TRACKMAP_FILE = "data/2024_10_7/osm/geojson/prague_tramtracks.geojson"
OSM_ROUTE_XML = "data/2024_10_7/osm/route_xml/lin22_hostivar2bilahora.xml"
OSM_BUILDINGS_DIR = "data/2024_10_7/osm/geojson"
OSM_BUILDINGS_PREFIX = "buildings_"


using Geodesy, GeoStats, GeoIO, LinearAlgebra, NearestNeighbors, Dates, TimesDates, Meshes
using LinearAlgebra

if !@isdefined(UrbanCanyonGNSS)
    include("../src/UrbanCanyonGNSS.jl")
    using .UrbanCanyonGNSS
else
    println("UrbanCanyonGNSS already loaded.")
end


include("../src/mpmitigation/elanglelookup.jl")

trackmap = loadmap(OSM_TRACKMAP_FILE, OSM_ROUTE_XML)
trackmap = resamplemap(trackmap; step = 5)
println("Track map loaded with ", length(trackmap.waypoints), " waypoints.")

walls = wallsfromdirectory(OSM_BUILDINGS_DIR, trackmap.center; filestartswith = OSM_BUILDINGS_PREFIX)
println("3D map loaded with ", length(walls), " walls.")

pointified_walls = [pointify(wall) for wall in walls]
heights = [max([wall[i].coords.z.val for i in 1:4]...) for wall in pointified_walls]
w2d = [[GeoStats.Point(wall[1].coords.x.val,wall[1].coords.y.val), GeoStats.Point(wall[2].coords.x.val,wall[2].coords.y.val)] for wall in pointified_walls]


using Distributed
println("SPAWNING WORKERS")
addprocs(6)  # or julia -p n
@everywhere using Geodesy, GeoStats, GeoIO, LinearAlgebra, NearestNeighbors, Dates, TimesDates
@everywhere gs = GeoStats
@everywhere using LinearAlgebra
@everywhere if !@isdefined(UrbanCanyonGNSS)
    include("../src/UrbanCanyonGNSS.jl")
    using .UrbanCanyonGNSS
else
    println("UrbanCanyonGNSS already loaded.")
end
println("WORKERS SPAWNED")

@everywhere function minelangles(wpid)
    println("WP: ",wpid," / ", length(trackmap.waypoints))
    wp = trackmap.waypoints[wpid]
    n = 32
    raylen = 100

    mpt = zeros(n)
    p = gs.Point(wp.enu[1:2]...)
    for bin in 1:n
        α = 2π * (bin - 1) / n
        α2 = 2π * bin / n
        p2 = gs.Point(p.coords.x.val + raylen * sin(α), p.coords.y.val + raylen * cos(α))
        p3 = gs.Point(p.coords.x.val + raylen * sin(α2), p.coords.y.val + raylen * cos(α2))
        ray = Segment(p, p2)
        ray2 = Segment(p, p3)

        for (wid, w) in enumerate(w2d)
            if (norm(w[1]-p).val >2* raylen) && (norm(w[2]-p).val > 2 * raylen)
                continue  # Skip if the wall is too far from the ray
            end
            intersections = [
                intersection(Segment(w...), ray).geom, 
                intersection(Segment(w...), ray2).geom,
                w[1], w[2],
                nothing
                ]
            if !intersects(w[1], Triangle(p,p2,p3))
                intersections[3] = nothing
            end
            if !intersects(w[2], Triangle(p,p2,p3))
                intersections[4] = nothing
            end
            filter!(!isnothing, intersections)
            if !isempty(intersections)
                dists = [norm(p - x).val for x in intersections]
                elangle = atan(heights[wid], min(dists...))
                if elangle > mpt[bin]
                    mpt[bin] = elangle
                end
            end
        end
    end
    return mpt
end


println("Fetching data for workers")
@everywhere w2d = @fetchfrom 1 w2d
@everywhere heights = @fetchfrom 1 heights
@everywhere trackmap = @fetchfrom 1 trackmap
println("Data fetched")

#rows = trackmap.waypoints |> collect  # convert to vector of rows
mplookuptable = pmap(minelangles, [1:length(trackmap.waypoints)...]) # Run the computation in parallel


println("Multipath lookup table created with ", length(mplookuptable), " waypoints.")
using Serialization
mplookup = ElangleLookup(trackmap.kdtree, trackmap.ecef2enu, hcat(mplookuptable...)')
serialize("simulateddata/mplookuptable.dat", mplookup)

for w in workers()
    println("KILLING WORKER $w")
    rmprocs(w)
end
