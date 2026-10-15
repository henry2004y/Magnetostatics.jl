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
where ``r_0`` is the source-surface radius where the field is purely radial and starts to be wound up by the solar rotation (typically ``r_0 \approx 2.5 R_\odot``), ``B_0`` is the radial field strength at ``r_0``, ``\Omega_\odot \approx 2.865 \times 10^{-6}`` rad/s is the solar angular rotation rate, and ``v_{\text{sw}} \approx 400`` km/s is the solar wind speed.

!!! note
    Since ``B_\phi \propto (r - r_0)``, choosing ``r_0 = 1`` AU would make the field purely radial at Earth.
    To obtain the usual spiral angle of about ``45^\circ`` at 1 AU, place ``r_0`` at the source surface and scale ``B_0`` so that ``B_r`` matches the observed value at 1 AU.

## Basic Usage

```@example parker
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie, FieldTracer

AU = 1.496e11            # 1 AU in meters
R_sun = 6.96e8           # Solar radius in meters
r0 = 2.5 * R_sun         # Source-surface radius
B0 = 5.0e-9 * (AU / r0)^2 # Gives B_r = 5 nT at 1 AU
v_sw = 4.0e5             # 400 km/s
Omega = 2.865e-6         # Solar rotation rate [rad/s]

spiral = ParkerSpiral(B0, r0; v_sw, Omega)

# Field at 1 AU on the equator
r_earth = SVector(AU, 0.0, 0.0)
B_earth = spiral(r_earth)
println("B at 1 AU: ", B_earth, " [T]")

# Spiral angle psi = atan(|B_phi| / B_r)
psi = rad2deg(atan(abs(B_earth[2]) / B_earth[1]))
println("Parker spiral angle at Earth: ", round(psi, digits=1), "°")
```

## Visualizing the Heliospheric Field in the Ecliptic Plane

```@example parker
# Grid in units of AU
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

## Field Lines in 3D

Since ``B_\theta = 0``, every field line stays on a cone of constant colatitude and winds into an Archimedean-like spiral as it moves outward.
The spiral is tightest near the equator, where ``|B_\phi| / B_r \propto \sin\theta`` is largest, and opens up towards the poles.
We use FieldTracer.jl v0.2.4 or later to trace
``d\mathbf{x}/ds = \mathbf{B}/|\mathbf{B}|`` from the source surface at several
latitudes and longitudes. The callback accepts position and signed arc length;
the stationary field only needs position. Steps have a physical length in meters.
Tracing stops just inside 5 AU, without evaluating the field beyond that radius.
Pass an ordinary `Vector` seed (`collect(seed)` below), as v0.2.4's native
tracer does not support `SVector` seeds internally.

```@example parker
parker_field(p, s) = spiral(p)
outside_heliosphere(p, s) = norm(p) > 5.0 * AU

fig3d = Figure(size=(900, 800), fontsize=16)
ax3 = Axis3(fig3d[1, 1]; aspect=:data, xlabel="x [AU]", ylabel="y [AU]",
    zlabel="z [AU]", title="Parker Spiral Field Lines in 3D",
    azimuth=1.2π, elevation=0.35)

colormap = :coolwarm
latitudes = (-60, -30, 0, 30, 60)
for lat in latitudes, φ0 in range(0, 2π; length=7)[1:(end - 1)]
    θ = deg2rad(90 - lat)
    seed = r0 * SVector(sin(θ) * cos(φ0), sin(θ) * sin(φ0), cos(θ))
    result = trace(
        parker_field, collect(seed); ds = 0.002 * AU, max_length = 20.0 * AU,
        direction = "forward", isoutside = outside_heliosphere
    )
    @assert result.status == :outside
    line = result.points
    lines!(ax3, [Point3f(p / AU) for p in line]; color=fill(lat, length(line)),
        colormap, colorrange=(-60, 60), linewidth=2)
end

mesh!(ax3, Sphere(Point3f(0), 0.12f0); color=:gold) # Sun (not to scale)
φs = range(0, 2π; length=200)
lines!(ax3, cos.(φs), sin.(φs), zeros(200); color=:black, linestyle=:dash,
    linewidth=1.5) # Earth orbit
Colorbar(fig3d[1, 2]; limits=(-60, 60), colormap, label="Source latitude [°]")

fig3d
```

The dashed circle marks Earth's orbit at 1 AU.
Field lines rooted near the equator wrap around the Sun several times within 5 AU, while those from high latitudes remain nearly radial and form open, loosely wound helices along the rotation axis.
