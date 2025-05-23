function geoplot(lats,lons; color = :blue, style = :scatter, provider = Tyler.TileProviders.OpenStreetMap(:Mapnik), size = (900,900), markersize = 5, linewidth = 5, label = "")
    minlat = minimum(lats)
    maxlat = maximum(lats)
    minlon = minimum(lons)
    maxlon = maximum(lons)
    
    window = Rect(minlon, minlat, maxlon-minlon, maxlat-minlat)
    fig = Figure(; size = size, fonts = (; regular = "Times New Roman", bold = "Times New Roman", bold_italic = "Times New Roman", italic = "Times New Roman"))
    ax = Axis(fig[1,1])
    
    m = Tyler.Map(window; provider, figure=fig, axis=ax, crs=Tyler.wgs84)
    wait(m)
    if style == :scatter
        route = scatter!(m.axis, lons, lats; color = color, markersize = markersize, label = label)
    elseif style == :line
        route = lines!(m.axis, lons, lats; color = color, linewidth = linewidth, label = label)
    else
        error("Unknown style, plotting as scatter")
        route = scatter!(m.axis, lons, lats; color = color)
    end
    wait(m)
    return m, route, fig, ax # call translate!(route, 0, 0, 10) when it disappears... because... reasons...
end

function geoplot!(m,lats,lons; color = :red, style = :scatter,  markersize = 5, linewidth = 5, label = "")
    if style == :scatter
        route = scatter!(m.axis, lons, lats; color = color, markersize = markersize, label = label)
    elseif style == :line
        route = lines!(m.axis, lons, lats; color = color, linewidth = linewidth, label = label)
    else
        error("Unknown style, plotting as scatter")
        route = scatter!(m.axis, lons, lats; color = color)
    end
    
    wait(m)
    return m, route, nothing, nothing # call translate!(route, 0, 0, 10) when it disappears... because... reasons...
end


