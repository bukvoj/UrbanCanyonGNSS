function dist2map(map::TrackMap, point::ECEF)
    pt = map.ecef2enu(point)
    pt = [pt[1], pt[2]]


    idxs, dist = knn(map.kdtree, pt, 1, true)
    p1 = [map.waypoints[idxs[1]-1].enu...][1:2]
    p2 = [map.waypoints[idxs[1]].enu...][1:2]
    p3 = [map.waypoints[idxs[1]+1].enu...][1:2]

    proj1 = p1 + dot(pt - p1, p2 - p1) / dot(p2 - p1, p2 - p1) * (p2 - p1)
    proj2 = p2 + dot(pt - p2, p3 - p2) / dot(p3 - p2, p3 - p2) * (p3 - p2)

    d1 = norm(pt-proj1)
    d2 = norm(pt-proj2)

    if isnan(d1)
        println("d1 is nan", p1, p2, p3, pt)
    end
    if isnan(d2)
        println("d2 is nan", p1, p2, p3, pt)
    end

    return min(d1, d2)
end


function dist2map(map::TrackMap, point::LLA)
    pt = ECEF(point, wgs84)
    return dist2map(map, pt)
end

function dist2map(map::TrackMap, point::ENU)
    pt = [point[1], point[2]]

    idxs, dist = knn(map.kdtree, pt, 1, true)
    p1 = [map.waypoints[idxs[1]-1].enu...][1:2]
    p2 = [map.waypoints[idxs[1]].enu...][1:2]
    p3 = [map.waypoints[idxs[1]+1].enu...][1:2]

    proj1 = p1 + dot(pt - p1, p2 - p1) / dot(p2 - p1, p2 - p1) * (p2 - p1)
    proj2 = p2 + dot(pt - p2, p3 - p2) / dot(p3 - p2, p3 - p2) * (p3 - p2)

    d1 = norm(pt-proj1)
    d2 = norm(pt-proj2)

    if isnan(d1)
        println("d1 is nan", p1, p2, p3, pt)
    end
    if isnan(d2)
        println("d2 is nan", p1, p2, p3, pt)
    end

    return min(d1, d2)
end


function proj2map(trackmap::TrackMap, point::ECEF)
    """
    takes:
    - trackmap: TrackMap
    - point: ECEF point
    returns:
    - point on the trackmap closest to the point in ENU coordinates
    """

    pt = trackmap.ecef2enu(point)
    pt = [pt[1], pt[2]]

    idxs, dist = knn(trackmap.kdtree, pt, 2, true)
    p1 = [trackmap.waypoints[idxs[1]].enu...][1:2]
    p2 = [trackmap.waypoints[idxs[2]].enu...][1:2]

    proj = p1 + dot(pt - p1, p2 - p1) / dot(p2 - p1, p2 - p1) * (p2 - p1)

    # TODO transform from ENU to LLA or something that I can actually use....
end

function proj2map(trackmap::TrackMap, point::ENU)
    """
    takes:
    - trackmap: TrackMap
    - point: ENU point
    returns:
    - point on the trackmap closest to the point in ENU coordinates
    """
    pt = [point[1], point[2]]

    idxs, dist = knn(trackmap.kdtree, pt, 2, true)
    p1 = [trackmap.waypoints[idxs[1]].enu...][1:2]
    p2 = [trackmap.waypoints[idxs[2]].enu...][1:2]

    proj = p1 + dot(pt - p1, p2 - p1) / dot(p2 - p1, p2 - p1) * (p2 - p1)
    proj = ENU(proj[1], proj[2], 0.0)
    return proj
end