using Geodesy
using NearestNeighbors

struct Waypoint
    lla::LLA
    ecef::ECEF
    enu::ENU
end

struct TrackMap
    waypoints::Vector{Waypoint}
    center::LLA
    kdtree::KDTree
    lla2enu::Geodesy.CoordinateTransformations.ComposedTransformation{ENUfromECEF{Float64}, ECEFfromLLA}
    ecef2enu::ENUfromECEF{Float64}
end

function TrackMap(tmap)
    # map is a vector of LLA. We assume it is ordered.
    tmap = unique(tmap)
    center = mean([pt.lat, pt.lon, pt.alt] for pt in tmap)
    center = LLA(center...)

    enutrans = ENUfromLLA(center, wgs84)
    ecefenutrans = ENUfromECEF(ECEF(center, wgs84), wgs84)
    waypoints = [Waypoint(pt, ECEF(pt, wgs84), enutrans(pt)) for pt in tmap]

    kdtree = KDTree(hcat([[pt.enu[1], pt.enu[2]] for pt in waypoints]...))
    return TrackMap(waypoints, center, kdtree, enutrans, ecefenutrans)
end