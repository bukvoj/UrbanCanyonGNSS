using Pkg
Pkg.activate("..")

Pkg.add("GLMakie")
Pkg.add("CairoMakie")

Pkg.add("DataFrames")
Pkg.add("Geodesy")
Pkg.add("GeoStats")
Pkg.add("XML")
Pkg.add("NearestNeighbors")
Pkg.add("LowLevelParticleFilters")
Pkg.add("Statistics")
Pkg.add("Distributions")
Pkg.add("LinearAlgebra")
Pkg.add("Serialization")
Pkg.add("Dates")
Pkg.add("TimesDates")
Pkg.add("Tyler")
# My packages are not registered yet, so we use the git url
Pkg.add(url="https://github.com/bukvoj/RinexRead.jl")
Pkg.add(url="https://github.com/bukvoj/GNSSEphemeris.jl")
Pkg.add(url="https://github.com/bukvoj/Klobuchar.jl")
Pkg.add(url="https://github.com/bukvoj/GNSSMultipathSim.jl")
