# [Harris Sheet](@id harris_sheet_example)

A current sheet model often used in space physics ($B_x(z) = B_0 \tanh(z/L)$).

```@example harris
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

B0 = 1.0       # Asymptotic field strength
L = 2.0        # Half-width
sheet = HarrisSheet(B0, L)

r = SVector(0.0, 0.0, 1.0)
B = sheet(r)
println("B at $r: $B [T]")
```

Visualizing the field reversal:

```@example harris
xs = range(-5, 5, length=51)
zs = range(-5, 5, length=51)

# B field is only in x direction, varying with z.
# But for streamplot let's show it in xz plane. By=0, Bz=0.
# Wait, Harris sheet usually has guide field or just Bx(z).
# Ideally B = (B0 tanh(z/L), 0, 0).
# Streamlines are just horizontal lines. Heatmap is more useful.

fig = Figure(size = (700, 600), fontsize=20)
ax = Axis(fig[1, 1], xlabel="x", ylabel="z", title="Harris Sheet Field (Bx)")

Bx_vals = [sheet(SVector(x, 0.0, z))[1] for x in xs, z in zs] # Bx is constant in x, varies in z

hm = heatmap!(ax, xs, zs, Bx_vals, colormap=:balance)
Colorbar(fig[1, 2], hm, label="Bx")

# Add arrows to show direction
arrows!(ax, xs[1:5:end], zs[1:5:end], [1.0 for _ in xs[1:5:end], _ in zs[1:5:end]], zeros(length(xs[1:5:end]), length(zs[1:5:end])), arrowsize=10, lengthscale=0.3, color=:black)

fig
```
