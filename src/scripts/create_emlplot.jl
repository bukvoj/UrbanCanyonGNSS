using CairoMakie
using GLMakie
GLMakie.activate!()



earlyx = [-3.5,-3,-1,1,3.5]
earlyy = [0,0,1,0,0]

latex = -earlyx
latey = earlyy


f = Figure(size=(600, 600),fonts = (; regular = "Times New Roman", bold = "Times New Roman", bold_italic = "Times New Roman", italic = "Times New Roman"))
ax = Axis(f[1,1], ylabel = "Correlation", xlabel = "Delay [chips]", xticks = [-3:1:3...],)
# lines!(ax,[-3,3],[0,0];color=:grey)
lines!(ax, earlyx,earlyy;color=:black)
lines!(ax, latex,latey;color=:black)
limits!(ax, -3.5, 3.5, -0.5, 1.5)
save("results/emlplot.pdf", f, size = (300, 200), backend = CairoMakie, update = false)

emlx = [-3.5,-3,-1,1,3,3.5]
emly = [0,0,1,-1,0,0]

f = Figure(size=(600, 600),fonts = (; regular = "Times New Roman", bold = "Times New Roman", bold_italic = "Times New Roman", italic = "Times New Roman"))
ax = Axis(f[1,1], ylabel = "Early-Late Correlation", xlabel = "Delay [chips]", xticks = [-3:1:3...],)
# lines!(ax,[-3,3],[0,0];color=:grey)
lines!(ax, emlx,emly;color=:black)
limits!(ax, -3.5, 3.5, -1.5, 1.5)
save("results/emlplot2.pdf", f, size = (300, 200), backend = CairoMakie, update = false)