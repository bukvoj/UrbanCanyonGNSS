using GeoIO, GeoStats, XML, Geodesy

function loadmap(mapfile, routefile)
    """
    Load the map of the route.
    - mapfile: path to the geojson of the map downloaded from OSM
    - routefile: xml including the relation of the route
    """

    osmmap = GeoIO.load(mapfile)
    routexml = XML.read(routefile, XML.Node)

    function extracttrammap(trackmap, wayids)
        track = trackmap |> Filter(row -> parse(Int, split(row["@id"], "/")[2]) in wayids)
        return track
    end

    function extractids(route_xml)
        stopids = []
        wayids = []

        for item in children(route_xml[1][1])
            if item.tag == "member"
                if item.attributes["role"] == "stop"
                    push!(stopids, parse(Int, item.attributes["ref"]))
                elseif item.attributes["type"] == "way"
                    push!(wayids, parse(Int, item.attributes["ref"]))
                end
            end
        end
        return wayids
    end

    wayids = extractids(routexml)
    ways = extracttrammap(osmmap, wayids).geometry

    # order the ways
    ordered = connectropes(ways)

    # convert to LLA
    latlons = pointify(ordered)
    lon = [pt.coords.lon.val for pt in latlons]
    lat = [pt.coords.lat.val for pt in latlons]

    # construct the trackmap
    trackmap = [LLA{Float64}(lat[i], lon[i], 0) for i in eachindex(lat)]
    trackmap = unique(trackmap)
    trackmap = TrackMap(trackmap)
end



function connectropes(ways)
    tmp = [ways[i] for i in eachindex(ways)]
    for i in 1:(length(ways)-1)
        r = tmp[1]
        dist = Inf
        closest = 0
        connectiontype = "00"
        len = length(tmp)
        for j in 2:len
            r2 = tmp[j]
            r2start = pointify(r2)[1]
            r2end = pointify(r2)[end]
            
            dist11 = norm(pointify(r)[1] - r2start).val
            dist22 = norm(pointify(r)[end] - r2end).val
            dist12 = norm(pointify(r)[1] - r2end).val
            dist21 = norm(pointify(r)[end] - r2start).val
            if dist11 < dist
                dist = dist11
                closest = j
                connectiontype = "11"   
            end
            if dist22 < dist
                dist = dist22
                closest = j
                connectiontype = "22"
            end
            if dist12 < dist
                dist = dist12
                closest = j
                connectiontype = "12"
            end
            if dist21 < dist
                dist = dist21
                closest = j
                connectiontype = "21"
            end
        end
        secondrope = pointify(tmp[closest])
        if connectiontype == "11" || connectiontype == "22"
            secondrope = reverse(secondrope)
        end

        if connectiontype == "11" || connectiontype == "12"
            newrope = Rope(vcat(secondrope, pointify(r)))
        else
            newrope = Rope(vcat(pointify(r), secondrope))
        end
        tmp[1] = newrope
        deleteat!(tmp, closest)
    end
    return tmp
end


function resamplemap(trackmap::TrackMap; step = 5)
    """
    Resample the trackmap to a fixed step size.
    - trackmap: the trackmap to resample
    - step: the step size
    """
    pts = []
    for i in 1:length(trackmap.waypoints)-1
        
        p1 = trackmap.waypoints[i].enu
        p2 = trackmap.waypoints[i+1].enu

        push!(pts, LLA(p1, trackmap.center, wgs84))
        while norm(p1-p2) > step
            p1 = p1 + step * (p2-p1)/norm(p2-p1)
            push!(pts, LLA(p1, trackmap.center, wgs84))
        end
    end
    return TrackMap(pts)
end