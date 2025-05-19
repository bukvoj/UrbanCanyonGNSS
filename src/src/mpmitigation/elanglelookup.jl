struct ElangleLookup
    knntree::KNNTree
    elangle::Matrix{Float64}
end

function usempmap!(R, mpmap, kf, az::AbstractArray, el::AbstractArray)
    numsvs = floor(Int, size(R,1) / 4)
    pos = kf.x[1:3]

    # Find index of closest point on map
    idxs, dist = knn(mpmap.kdtree, pt, 1, true)

    # The lookup table was started at point 500
    id = idxs[1]
    if id < 1
        return
    end

    # number of bins 
    numbins = size(mpmap, 2)

    # find bin for each az-angle
    az = mod.(az .- 360/numbins/2, 360)
    bins = Int.(floor.(az/360 .* numbins)) .+ 1

    # increase noise for each sv with elevation lower than the mask given by mpmap
    for j in eachindex(el)
        if mpmap.lookuptable[id][bins[j]] > el[j]*pi/180
            for i in 0:3
                R[i*numsvs + j, i*numsvs + j] *= 100 # increase noise for NLOS measurements
            end
        end
    end
end