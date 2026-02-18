# [Magnetic Dipole](@id dipole_example)

Calculate the field of a magnetic dipole moment $\mathbf{M} = (0,0,1)$ at various points.

```@example dipole
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

M = SVector(0.0, 0.0, 1.0)
dipole = Dipole(M)

r = SVector(1.0, 0.0, 0.0) # Point on x-axis
B = dipole(r)
println("B at $r: $B [T]")
```

Visualizing the field in the xz-plane:

```@example dipole
xs = range(-2, 2, length=51)
zs = range(-2, 2, length=51)

function field_xz(x::T, z::T) where T
    B = dipole(SVector(x, zero(T), z))
    return Point(B[1], B[3])
end

fig = Figure(size = (700, 600), fontsize=20)
ax = Axis(fig[1, 1]; xlabel="x", ylabel="z", aspect=DataAspect(), title="Magnetic Dipole Field")

# Heatmap of magnitude
Bmag = [norm(dipole(SVector(x, 0.0, z))) for x in xs, z in zs]
hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1e-9), colormap=:plasma)
Colorbar(fig[1, 2], hm, label="log10(|B|)")

# Field lines
streamplot!(ax, field_xz, -2..2, -2..2; arrow_size = 8, linewidth = 1.5)

fig
```
