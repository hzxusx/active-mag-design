using Plots
plotly() # Choose the Plotly.jl backend for web interactivity
# # plot(rand(5,5),linewidth=2,title="My Plot")


X, Y, Z = G.xyzs
sample_B_wire = min_state[3]
B_diff_norm, B_diff, I = evalB(G, B_object, sample_B_wire)

B_diff_mag = [norm(B_diff[:, i]) for i in 1:G.points]
B_diff_y = [B_diff[2, i] for i in 1:G.points]
B_diff_z = [B_diff[3, i] for i in 1:G.points]
B_diff_2d = reshape(B_diff_mag, size(Y, 1), size(Z, 1))
contour_p = contour(Y,Z,B_diff_2d,fill=true)
plot(contour_p)


# quiver_p = quiver(Y,Z,B_diff_y,B_diff_z)
# plot(quiver_p)
# plot(Z, B_diff_2d[:, 16])

# # ---
# using PyPlot
# x = linspace(0,2*pi,1000); y = sin.(3*x + 4*cos.(2*x))
# plot(x, y, color="red", linewidth=2.0, linestyle="--")