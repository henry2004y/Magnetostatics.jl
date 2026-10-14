# [Guide-Field Current Sheet](@id guide_field_sheet_example)

In magnetic reconnection studies (such as at the Earth's magnetopause, in tokamak tearing modes, and in laboratory experiments), the reconnecting current sheet often carries an out-of-plane guide field ``B_g`` and a normal field ``B_n``:
```math
B_x(z) = B_0 \tanh\left(\frac{z}{L}\right), \quad B_y(z) = B_g, \quad B_z(z) = B_n
```
The corresponding current density is directed along the ``y``-axis:
```math
J_y(z) = \frac{B_0}{\mu_0 L} \operatorname{sech}^2\left(\frac{z}{L}\right)
```

## Basic Usage

```@example guide_field
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

B0 = 20.0e-9   # Asymptotic reconnecting field [T] (20 nT)
L = 500.0e3    # Sheet half-thickness [m] (500 km)
Bg = 10.0e-9   # Guide field [T] (10 nT)
Bn = 2.0e-9    # Normal field [T] (2 nT)

sheet = GuideFieldHarrisSheet(B0, L; Bg=Bg, Bn=Bn)

# Center of current sheet
r0 = SVector(0.0, 0.0, 0.0)
println("B at sheet center: ", sheet(r0), " [T]")
println("J at sheet center: ", current_density(sheet, r0), " [A/m²]")
```

## Magnetic Field and Shear Profile

```@example guide_field
zs = range(-3 * L, 3 * L, length=100)

Bx_vals = [sheet(SVector(0.0, 0.0, z))[1] * 1e9 for z in zs]
By_vals = [sheet(SVector(0.0, 0.0, z))[2] * 1e9 for z in zs]
Bz_vals = [sheet(SVector(0.0, 0.0, z))[3] * 1e9 for z in zs]
Bmag = [norm(sheet(SVector(0.0, 0.0, z))) * 1e9 for z in zs]

fig = Figure(size=(750, 500), fontsize=18)
ax = Axis(fig[1, 1], xlabel="z / L", ylabel="Field Components [nT]",
    title="Guide-Field Harris Current Sheet Profiles")

lines!(ax, zs ./ L, Bx_vals, label="B_x (reconnecting)", linewidth=3, color=:blue)
lines!(ax, zs ./ L, By_vals, label="B_y (guide field)", linewidth=3, color=:red, linestyle=:dash)
lines!(ax, zs ./ L, Bz_vals, label="B_z (normal)", linewidth=2, color=:green, linestyle=:dot)
lines!(ax, zs ./ L, Bmag, label="|B|", linewidth=3, color=:black)

axislegend(ax, position=:rt)
fig
```
