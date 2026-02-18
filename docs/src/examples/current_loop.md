# [Current Loop](@id current_loop_example)

Calculate the field of a circular current loop.

```@example loop
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

# Define loop parameters
radius = 1.0
current = 10.0
center = SVector(0.0, 0.0, 0.0)
normal = SVector(0.0, 0.0, 1.0)

# Create loop object
loop = CurrentLoop(radius, current, center, normal)

# Query field
r = SVector(0.0, 0.0, 0.5)
B = getB_loop(r, loop)
println("B at (0,0,0.5): $B [T]")
```

Visualizing the field in the xz-plane:

```@example loop
xs = range(-2, 2, length=51)
zs = range(-2, 2, length=51)

function field_xz_loop(x, z)
    B = getB_loop(SVector(x, 0.0, z), loop)
    # Return nan near wire? Usually streamplot handles huge values okay, but let's be safe
    return Point2f(B[1], B[3])
end

fig = Figure(size = (700, 600), fontsize=20)
ax = Axis(fig[1, 1]; xlabel="x", ylabel="z", aspect=DataAspect(), title="Current Loop Field")

# Heatmap of magnitude
Bmag = [norm(getB_loop(SVector(x, 0.0, z), loop)) for x in xs, z in zs]
hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1e-9), colormap=:plasma)
Colorbar(fig[1, 2], hm, label="log10(|B|)")

# Field lines
streamplot!(ax, field_xz_loop, -2..2, -2..2; arrow_size = 8, linewidth = 1.5, colormap = :turbo)

fig
```
