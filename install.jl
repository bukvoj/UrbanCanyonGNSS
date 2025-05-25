using Pkg
Pkg.activate(".")

# Yeah, the packages are not registered, 
# calling resolve resolves into error
# so we remove them first, then resolve and instantiate
# and then add them again
try
    Pkg.rm("RinexRead")
    Pkg.rm("GNSSEphemeris")
    Pkg.rm("Klobuchar")
    Pkg.rm("GNSSMultipathSim")
catch e
    nothing
end
Pkg.resolve()
Pkg.instantiate()

# add the packages again
Pkg.add(url="https://github.com/bukvoj/RinexRead.jl")
Pkg.add(url="https://github.com/bukvoj/GNSSEphemeris.jl")
Pkg.add(url="https://github.com/bukvoj/Klobuchar.jl")
Pkg.add(url="https://github.com/bukvoj/GNSSMultipathSim.jl")

# Create directories for results and simulated data if they do not exist
if !isdir("results/")
    mkdir("results/")
end
if !isdir("simulateddata/")
    mkdir("simulateddata/")
end