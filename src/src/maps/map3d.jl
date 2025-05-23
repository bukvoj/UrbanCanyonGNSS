struct Map3D
    center::GeoStats.Point
    chunksize::Float64
    map::Vector{Any}
    xchunks::Int
    ychunks::Int
    xmin::Float64
    ymin::Float64
end

function chunkifywalls(walls, center; chunksize = 150)
    xmin = center.coords.x.val
    ymin = center.coords.y.val
    xmax = center.coords.x.val
    ymax = center.coords.y.val
    for wall in walls
        wall = pointify(wall)
        xmin = min(xmin, wall[1].coords.x.val, wall[2].coords.x.val)
        ymin = min(ymin, wall[1].coords.y.val, wall[2].coords.y.val)
        xmax = max(xmax, wall[1].coords.x.val, wall[2].coords.x.val)
        ymax = max(ymax, wall[1].coords.y.val, wall[2].coords.y.val)
    end

    xchunks = ceil(Int, (xmax - xmin) / chunksize)
    ychunks = ceil(Int, (ymax - ymin) / chunksize)
    chunked = [[] for _ in 1:xchunks * ychunks]

    for wall in walls
        if typeof(wall) != typeof(walls[1])
            continue
        end
        w = pointify(wall)
        x = floor(Int, (w[1].coords.x.val - xmin) / chunksize)
        y = floor(Int, (w[1].coords.y.val - ymin) / chunksize)
        push!(chunked[x * ychunks + y], wall)
    end
    return Map3D(center, chunksize, chunked, xchunks, ychunks, xmin, ymin)
end



struct EnvMap 
    trackmap::TrackMap
    map3d::Map3D
end