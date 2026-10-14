# [Parker Spiral (Interplanetary Magnetic Field)](@id parker_spiral_example)

The Parker spiral model (Parker 1958) describes the large-scale configuration of the interplanetary magnetic field (IMF) carried out from a rotating star by the stellar wind.

In spherical coordinates ``(r, \theta, \phi)`` centered on the Sun:
```math
B_r(r, \theta) = B_0 \left( \frac{r_0}{r} \right)^2
```
```math
B_\phi(r, \theta) = -B_0 \left( \frac{r_0}{r} \right)^2 \left( \frac{\Omega_\odot (r - r_0)}{v_{\text{sw}}} \right) \sin\theta
```
```math
B_\theta = 0
```
where ``r_0`` is a reference radius (e.g. ``1 \text{ AU} \approx 1.496 \times 10^{11}`` m), ``B_0`` is the radial field strength at ``r_0``, ``\Omega_\odot \approx 2.865 \times 10^{-6}`` rad/s is the solar angular rotation rate, and ``v_{\text{sw}} \approx 400`` km/s is the solar wind speed.

## Basic Usage

```@example parker
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

B0 = 5.0e-9       # 5 nT at 1 AU
r0 = 1.496e11     # 1 AU in meters
v_sw = 4.0e5      # 400 km/s
Omega = 2.865e-6  # Solar rotation rate [rad/s]

spiral = ParkerSpiral(B0, r0; v_sw=v_sw, Omega=Omega)

# Field at 1 AU on the equator
r_earth = SVector(r0, 0.0, 0.0)
B_earth = spiral(r_earth)
println("B at 1 AU: ", B_earth, " [T]")

# Spiral angle psi = atan(|B_phi| / B_r)
psi = rad2deg(atan(abs(B_earth[2]) / B_earth[1]))
println("Parker spiral angle at Earth: ", round(psi, digits=1), "°")
```

## Visualizing the Heliospheric Field in the Ecliptic Plane

```@example parker
# Grid in units of AU
AU = 1.496e11
xs = range(-3.0 * AU, 3.0 * AU, length=25)
ys = range(-3.0 * AU, 3.0 * AU, length=25)

pts = Point2f[]
dirs = Vec2f[]

for x in xs, y in ys
    r_val = hypot(x, y)
    if r_val > 0.3 * AU
        B = spiral(SVector(x, y, 0.0))
        push!(pts, Point2f(x / AU, y / AU))
        b_vec = normalize(Vec2f(B[1], B[2]))
        push!(dirs, b_vec)
    end
end

fig = Figure(size=(700, 650), fontsize=18)
ax = Axis(fig[1, 1], xlabel="x [AU]", ylabel="y [AU]", aspect=DataAspect(),
    title="Parker Spiral Magnetic Field in the Ecliptic Plane")

arrows!(ax, pts, dirs; color=:deepskyblue3, arrowsize=9, lengthscale=0.18)
scatter!(ax, [Point2f(0, 0)], color=:gold, markersize=18, label="Sun")
axislegend(ax, position=:rt)

fig
```
