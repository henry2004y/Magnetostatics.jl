# [Magnetic Mirror](@id magnetic_mirror_example)

Two co-axial current loops separated by a distance.

```@example mirror
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

distance = 2.0
a = 1.0 # Radius
I = 100.0 # Current * windings

# Field at the center (0,0,0)
B = getB_mirror(0.0, 0.0, 0.0, distance, a, I)
println("B at center: $B [T]")
```

Visualizing the field in the xz-plane:

```@example mirror
xs = range(-2, 2, length=51)
zs = range(-2, 2, length=51)

function field_xz_mirror(x, z)
    B = getB_mirror(x, 0.0, z, distance, a, I)
    return Point2f(B[1], B[3])
end

fig = Figure(size = (700, 600), fontsize=20)
ax = Axis(fig[1, 1]; xlabel="x", ylabel="z", aspect=DataAspect(), title="Magnetic Mirror Field")

Bmag = [norm(getB_mirror(x, 0.0, z, distance, a, I)) for x in xs, z in zs]
hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1e-9), colormap=:plasma)
Colorbar(fig[1, 2], hm, label="log10(|B|)")

streamplot!(ax, field_xz_mirror, -2..2, -2..2; arrow_size = 8, linewidth = 1.5)

fig
```
