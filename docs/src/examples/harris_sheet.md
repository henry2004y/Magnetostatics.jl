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
zs = range(-2, 2, length=21)

Bx_vals = [sheet(SVector(x, 0.0, z))[1] for x in xs, z in zs]

fig = Figure(size = (800, 400), fontsize=20)
ax = Axis(fig[1, 1];
    xlabel="x", ylabel="z", aspect=DataAspect(), title="Harris Sheet Field (Bx)")

hm = heatmap!(ax, xs, zs, Bx_vals, colormap=:balance)
Colorbar(fig[1, 2], hm, label="Bx")

ps = [Point2f(x, z) for x in xs[5:5:end-5], z in zs[1:3:end]]
ns = [Vec2f(sheet(SVector(x, 0.0, z))[1], sheet(SVector(x, 0.0, z))[3])
    for x in xs[5:5:end-5], z in zs[1:3:end]]

arrows2d!(ax, vec(ps), vec(ns); lengthscale=0.6, color=vec(Bx_vals), colormap=:rain)

fig
```
