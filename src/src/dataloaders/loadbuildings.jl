# DEPENDENCIES: GeoIO, GeoStats
using GeoIO, GeoStats
using Geodesy

function buildings2walls(buildings, refpoint_lla)
    refpoint_ecef = refpoint_lla |> Proj(Cartesian{WGS84Latest})

    # convert to ENU
    enu_buildings = buildings.geometry |> 
        Proj(Cartesian) |> 
        Translate(-to(refpoint_ecef)...) |> 
        Rotate(RotMatrix{3}(ecef2enurotmat(refpoint_lla.coords.lat.val, refpoint_lla.coords.lon.val))) |>
        Scale(1.0,1.0,1.0e-99) # SIMULATED EARTH IS FLAT

    # estimate the heights of buildings
    ls = buildings[:,"building:levels"]
    ls[ismissing.(ls)] .= "4"
    ls[ls .== "1"] .= "2.5"
    ls = parse.(Float64, ls)
    hs = 2. .+ ls*4

    walls = []
    for (id,b) in enumerate(enu_buildings)
        if !(typeof(b) <: PolyArea) # PolyAreas not supported - YET
            continue
        end
        ring = pointify(b.rings[1])
        l = length(ring)
        for i in 1:l
            wall = Ngon(ring[i], ring[(i+1)%l], ring[(i+1)%l] |> Translate(0,0,hs[id]), ring[i] |> Translate(0,0,hs[id]))
            if norm(ring[i] - ring[(i+1)%l]) < 1e-2u"m"
                continue
            end
            push!(walls, wall)
        end
    end
    println("Number of walls: ", length(unique(walls)))
    return unique(walls)
end

function buildingsfromdirectory(path::String; filestartswith::String = "", verbose::Bool = false)
    """
    Load buildings from a directory. The buildings are in GeoJSON format.
    The buildings are converted to walls in ENU coordinates.
    
    args:
    path::String: path to the directory
    filestartswith::String: prefix of the files to load
    verbose::Bool: if true, print the names of the files being loaded
    """
    bs = []
    for pth in readdir(path)
        if startswith(pth, filestartswith)
            if verbose
                println("Loading $pth")
            end
            push!(bs, GeoIO.load(path * "/" * pth))
        end
    end
    bs = vcat(bs...)
    return bs
end


function wallsfromdirectory(path::String, center_lla::LLA; filestartswith::String = "", verbose::Bool = false)
    wallsfromdirectory(path, [center_lla.lat,center_lla.lon,center_lla.alt]; filestartswith=filestartswith, verbose=verbose)
end

function wallsfromdirectory(path::String, center_lla::AbstractArray; filestartswith::String = "", verbose::Bool = false)
    """
    Load walls from a directory. The walls are in GeoJSON format.
    The walls are converted to ENU coordinates.
    
    args:
    path::String: path to the directory
    center_lla::vector: center of the walls in LLA coordinates
    filestartswith::String: prefix of the files to load
    verbose::Bool: if true, print the names of the files being loaded
    """
    bs = buildingsfromdirectory(path, filestartswith=filestartswith, verbose=verbose)
    refpoint = GeoStats.Point{🌐}(GeodeticLatLonAlt(center_lla[1], center_lla[2], center_lla[3]))
    walls = buildings2walls(bs, refpoint)
    return walls
end