# [Tokamak Configurations](@id tokamak_example)

## Tokamak Coils

Magnetic field from a Tokamak topology consisting of 16 toroidal field coils and a plasma current.

```@example tokamak
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

a = 1.0  # Coil radius
b = 2.0  # Major radius offset
ICoils = 1000.0
IPlasma = 500.0

# Query point
x, y, z = 3.0, 0.0, 0.0
B = getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma)
println("B at ($x, $y, $z): $B [T]")
```

Visualizing the poloidal field (xz-plane):

```@example tokamak
xs = range(0.5, 4.5, length=51)
zs = range(-2, 2, length=51)

function field_xz_tokamak(x, z)
    B = getB_tokamak_coil(x, 0.0, z, a, b, ICoils, IPlasma)
    return Point2f(B[1], B[3])
end

Bmag = [norm(getB_tokamak_coil(x, 0.0, z, a, b, ICoils, IPlasma)) for x in xs, z in zs]

fig = Figure(size = (700, 600), fontsize=20)
ax = Axis(fig[1, 1], xlabel="x (R)", ylabel="z", aspect=DataAspect(), title="Tokamak Coil Field (Poloidal)")

hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1e-9), colormap=:plasma)
Colorbar(fig[1, 2], hm, label="log10(|B|)")

streamplot!(ax, field_xz_tokamak, xs[1]..xs[end], zs[1]..zs[end]; arrow_size = 8, linewidth = 1.5)

fig
```

## Tokamak with q-profile

Reconstruct the magnetic field distribution from a safety factor ($q$) profile.

```@example tokamak_q
# Define q-profile as a function of normalized radius r/a
q_profile(r_norm) = 1.1 + r_norm^2

a = 1.0   # Minor radius
R0 = 3.0  # Major radius
B0 = 2.0  # Toroidal field on axis

# Query point inside plasma
x, y, z = 3.5, 0.0, 0.0
B = getB_tokamak_profile(x, y, z, q_profile, a, R0, B0)
println("B at ($x, $y, $z): $B [T]")
```

Visualizing the q-profile field:

```@example tokamak_q
xs = range(R0 - a, R0 + a, length=51)
zs = range(-a, a, length=51)

function field_xz_tokamak_q(x, z)
    # Check if inside plasma
    r_local = sqrt((x - R0)^2 + z^2)
    if r_local > a
        return Point2f(NaN, NaN)
    end
    B = getB_tokamak_profile(x, 0.0, z, q_profile, a, R0, B0)
    return Point2f(B[1], B[3])
end

# Check if inside plasma
function get_Bmag_tokamak_q(x, z)
    if sqrt((x - R0)^2 + z^2) > a
        return NaN
    end
    return norm(getB_tokamak_profile(x, 0.0, z, q_profile, a, R0, B0))
end
Bmag = [get_Bmag_tokamak_q(x, z) for x in xs, z in zs]

fig = Figure(size = (700, 600), fontsize=20)
ax = Axis(fig[1, 1], xlabel="x (R)", ylabel="z", aspect=DataAspect(), title="Tokamak q-profile Field")

hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1e-9), colormap=:plasma)
Colorbar(fig[1, 2], hm, label="log10(|B|)")

streamplot!(ax, field_xz_tokamak_q, xs[1]..xs[end], zs[1]..zs[end]; arrow_size = 8, linewidth = 1.5)
fig
```
