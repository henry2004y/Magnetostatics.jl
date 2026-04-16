# [Magnetic Mirror](@id magnetic_mirror_example)

Two co-axial current loops of radius $a$ carrying current $I$ separated by a distance $d$. The magnetic field is the superposition of the fields from the two individual loops:
```math
\mathbf{B}(\mathbf{r}) = \mathbf{B}_{\text{loop1}}(\mathbf{r}) + \mathbf{B}_{\text{loop2}}(\mathbf{r})
```

```@example mirror
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie, UniformStreamlines

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
ax = Axis(fig[1, 1];
    xlabel="x", ylabel="z", aspect=DataAspect(), title="Magnetic Mirror Field")

Bmag = [norm(getB_mirror(x, 0.0, z, distance, a, I)) for x in xs, z in zs]
hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1e-9), colormap=:plasma)
Colorbar(fig[1, 2], hm, label="log10(|B|)")

str = evenstream(xs, zs, (x, z) -> field_xz_mirror(x, z)[1], (x, z) -> field_xz_mirror(x, z)[2])
streamlines!(ax, str; linewidth = 1.5, with_arrows = true)

fig
```
