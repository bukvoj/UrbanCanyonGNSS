using GLMakie, CairoMakie #  # for plotting
GLMakie.activate!() # Activate GLMakie for plotting
using DataFrames, Dates, TimesDates, Geodesy, Serialization

using LinearAlgebra

using RinexRead, Serialization

# LOAD THE PACKAGE!
if !@isdefined(UrbanCanyonGNSS)
    include("../src/UrbanCanyonGNSS.jl")
    using .UrbanCanyonGNSS
else
    println("UrbanCanyonGNSS already loaded.")
end


# NAV_FILE = "../data/2024_10_7/rinex/line22_fromhostivartopohorelec.24N"
# OBS_FILE = "../data/2024_10_7/rinex/line22_fromhostivartopohorelec.24O"
SIM_FILE = "simulateddata/multipathdata.dat"

println("Loading data files...")
# nav = rinexread(NAV_FILE)
# obs = rinexread(OBS_FILE)
sim = deserialize(SIM_FILE)
println("Observations loaded.")

start = now()
# Prepare the observations for processing
using Random # for reproducibility
Random.seed!(20001129)
noisymeasurements, biases, t_offset0, t_offsetrate = genmeas(sim)
println("Pure IEKF simulation")
simresults = trajectory(noisymeasurements)
println("Pure IEKF simulated results calculated.")

println("Calculating sim with mpestimation iekf...")
simresults_mp = trajectory(noisymeasurements; runmpestimation = true)
println("Simulated multipath estimation results calculated.")
println("Duration: ", now()-start)



enu_iekf = unique(ENU.([x[1:3] for x in simresults.x]))
enu_mp = unique(ENU.([x[1:3] for x in simresults_mp.x]))


purepos = unique(ENU.([[x[1:2]...,0.] for x in simresults.x]))
puresigma = unique(x[1:3,1:3] for x in simresults.Σ)
errs = purepos .- unique(noisymeasurements.RecPos)
errpuredir = errs ./ norm.(errs)
errpuresigmadir = [x' for x in errpuredir] .* puresigma .* errpuredir

mppos = unique(ENU.([[x[1:2]...,0.] for x in simresults_mp.x]))
mpsigma = unique(x[1:3,1:3] for x in simresults_mp.Σ)
errsmp = mppos .- unique(noisymeasurements.RecPos)
errmpdir = errsmp ./ norm.(errsmp)
errmpsigmadir = [x' for x in errmpdir] .* mpsigma .* errmpdir 

println("RMS error without MP estimation: ", rms(norm.(errs)))
println("RMS error with MP estimation:    ", rms(norm.(errsmp)))


f = Figure(size=(600, 600),fonts = (; regular = "Times New Roman", bold = "Times New Roman", bold_italic = "Times New Roman", italic = "Times New Roman"))
Label(f[0,:], text = "Position error and its covariance \n (IEKF with MP estimation)", fontsize = 15, tellwidth = false)
ax1 = Axis(f[1, 1]; xlabel="Time [s]", ylabel="Position \n error [m]", limits = (-5,155,-5,35), yticks=[0,10,20,30])
ax2 = Axis(f[2, 1]; xlabel="Time [s]", ylabel="σ^2 [m]", limits = (-5,155,-1,7), yticks=[0,3,6])

plot!(ax1, norm.(errsmp); color = :red)
plot!(ax2, errmpsigmadir; color = :red)
save("results/simcovmp.pdf", f , size = (300, 300), update = false, backend = CairoMakie)

f = Figure(size=(600, 600), fonts = (; regular = "Times New Roman", bold = "Times New Roman", bold_italic = "Times New Roman", italic = "Times New Roman"))
Label(f[0,:], text = "Position Error and Its Covariance \n (IEKF without MP estimation)", fontsize = 15, tellwidth = false)
ax1 = Axis(f[1, 1]; xlabel="Time [s]", ylabel="Position \n error [m]", limits = (-5,155,-5,35), yticks=[0,10,20,30])
ax2 = Axis(f[2, 1]; xlabel="Time [s]", ylabel="σ^2 [m]", limits = (-5,155,-1,7),yticks=[0,3,6])
plot!(ax1, norm.(errs); color = :blue)
plot!(ax2, errpuresigmadir; color = :blue)
save("results/simcovpure.pdf", f , size = (300, 300), update = false, backend = CairoMakie)

purevel = unique(ENU.([[x[4:5]...,0.] for x in simresults.x]))
noisymeasurements = combine(groupby(noisymeasurements, :Time)) do sdf
    first(sdf, 1)
end
velerr = purevel .- noisymeasurements.RecVel
velerrdir = velerr ./ norm.(velerr)
puresigma = unique(x[4:6,4:6] for x in simresults.Σ)
velerrsigmadir = [x' for x in velerrdir] .* puresigma .* velerrdir

mpvel = unique(ENU.([[x[4:5]...,0.] for x in simresults_mp.x]))
velerrmp = mpvel .- noisymeasurements.RecVel
velerrmpdir = velerrmp ./ norm.(velerrmp)
mpsigma = unique(x[4:6,4:6] for x in simresults_mp.Σ)
velerrmpsigmadir = [x' for x in velerrmpdir] .* mpsigma .* velerrmpdir


f = Figure(size=(600, 600),fonts = (; regular = "Times New Roman", bold = "Times New Roman", bold_italic = "Times New Roman", italic = "Times New Roman"))
Label(f[0,:], text = "Velocity error and its covariance \n (IEKF with MP estimation)", fontsize = 15, tellwidth = false)
ax1 = Axis(f[1, 1]; xlabel="Time [s]", ylabel="Velocity \n error [m]", limits = (-5,155,-0.1,0.8))
ax2 = Axis(f[2, 1]; xlabel="Time [s]", ylabel="σ^2 [m]", limits = (-5,155,-0.01,0.7))
plot!(ax1, norm.(velerrmp); color = :red)
plot!(ax2, velerrmpsigmadir; color = :red)
save("results/simcovmpvel.pdf", f , size = (300, 300), update = false, backend = CairoMakie)

f = Figure(size=(600, 600), fonts = (; regular = "Times New Roman", bold = "Times New Roman", bold_italic = "Times New Roman", italic = "Times New Roman"))
Label(f[0,:], text = "Velocity Error and Its Covariance \n (IEKF without MP estimation)", fontsize = 15, tellwidth = false)
ax1 = Axis(f[1, 1]; xlabel="Time [s]", ylabel="Velocity \n error [m]", limits = (-5,155,-0.1,0.8))
ax2 = Axis(f[2, 1]; xlabel="Time [s]", ylabel="σ^2 [m]", limits = (-5,155,-0.01,0.7))
plot!(ax1, norm.(velerr); color = :blue)
plot!(ax2, velerrsigmadir; color = :blue)
save("results/simcovpurevel.pdf", f , size = (300, 300), update = false, backend = CairoMakie)
