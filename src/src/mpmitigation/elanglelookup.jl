struct ElangleLookup2
    kdtree::KDTree
    ecef2enu::Any
    elangle::Matrix{Float64}
end

function usempmap!(R, mpmap, ecefpos, az::AbstractArray, el::AbstractArray)
    numsvs = length(az)
    measpersv = size(R, 1) ÷ numsvs

    enupos = mpmap.ecef2enu(ECEF(ecefpos...))  # Convert ECEF to ENU coordinates
    pt = enupos[1:2]  # Use only the x and y coordinates for the lookup
    # Find index of closest point on map
    idxs, dist = knn(mpmap.kdtree, pt, 1, true)

    # The lookup table was started at point 500
    id = idxs[1]
    if id < 1
        return
    end

    # number of bins 
    numbins = size(mpmap.elangle, 2)

    # find bin for each az-angle
    az = mod.(az .- 360/numbins/2, 360)
    bins = Int.(floor.(az/360 .* numbins)) .+ 1

    # increase noise for each sv with elevation lower than the mask given by mpmap
    for j in eachindex(el)
        if mpmap.elangle[id,bins[j]] > el[j]*pi/180
            for i in 0:(measpersv-1)
                R[i*numsvs + j, i*numsvs + j] *= 100 # increase noise for NLOS measurements
            end
        end
    end
end













function Meshes.center(r::Ray)
    println("Center of the ray: ", r.p)
    return r.p
end

function Meshes.intersects(r::Ray, q::Quadrangle)
    pts = pointify(q) # get the points of the quadrangle
    return Meshes.intersects(Meshes.Triangle(pts[1], pts[2], pts[3]), r) || Meshes.intersects(Meshes.Triangle(pts[1], pts[3], pts[4]), r) # check if the ray Meshes.intersects with the quadrangle
end

function Meshes.intersects(q::Quadrangle, r::Ray)
    return Meshes.intersects(r,q) # check if the ray Meshes.intersects with the triangle
end

function Meshes.intersection(r::Ray, q::Quadrangle)
    pts = pointify(q) # get the points of the quadrangle
    i1 = intersection(r, Meshes.Triangle(pts[1], pts[2], pts[3])) # check if the ray Meshes.intersects with the triangle
    if !isnothing(i1.geom) # check if the intersection point is valid
        return i1.geom # get the first intersection point
    else
        i2 = intersection(r, Meshes.Triangle(pts[1], pts[3], pts[4])) # check if the ray Meshes.intersects with the triangle
        if !isnothing(i2.geom) # check if the intersection point is valid
            return i2.geom # get the first intersection point
        else
            return nothing # no intersection point found
        end
    end
end

function Meshes.intersection(q::Quadrangle, r::Ray)
    return intersection(r,q) # check if the ray Meshes.intersects with the triangle
end