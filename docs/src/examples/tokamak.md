# [Tokamak Configurations](@id tokamak_example)

## Tokamak Coils

Magnetic field from a Tokamak topology consisting of ``N`` toroidal field coils and a plasma current. The total field is the sum of contributions from all coils and the plasma ring:
```math
\mathbf{B}(\mathbf{r}) = \sum_{i=1}^N \mathbf{B}_{\text{coil}, i}(\mathbf{r}) + \mathbf{B}_{\text{plasma}}(\mathbf{r})
```

```@example tokamak
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie, UniformStreamlines, GeometryBasics

"""
Construct the topology of Tokamak.
"""
function get_tokamak_topology(R0, a)
    nθ = LinRange(0, 2π, 30)
    nζ = LinRange(0, 2π, 30)
    points = [let (sinθ, cosθ) = sincos(θ), (sinζ, cosζ) = sincos(ζ)
                  Point3f(
                    (R0 + a * cosθ) * cosζ, (R0 + a * cosθ) * sinζ, a * sinθ)
              end for θ in nθ, ζ in nζ]
    faces = decompose(QuadFace{GLIndex}, Tesselation(Rect(0, 0, 1, 1), size(points)))

    return GeometryBasics.Mesh(vec(points), faces)
end

const a = 1.0  # Coil radius
const b = 2.0  # Major radius offset
const ICoils = 1000.0
const IPlasma = 500.0
const R_major = 2a + b
const R_minor = 2a

# Query point
x, y, z = 3.0, 0.0, 0.0
B = getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma)
println("B at ($x, $y, $z): $B [T]")
```

Visualizing the poloidal field (xz-plane):

```@example tokamak
xs = range(0.5, 4.5, length=51)
zs = range(-2, 2, length=51)

function field_xz_tokamak(x)
    # Check if inside plasma
    r_local = sqrt((x[1] - R_major)^2 + x[2]^2)
    if r_local > R_minor
        return Point2f(NaN, NaN)
    end
    B = getB_tokamak_coil(x[1], 0.0, x[2], a, b, ICoils, IPlasma)
    return Point2f(B[1], B[3])
end

# Check if inside plasma
function get_Bmag_tokamak(x, z)
    if sqrt((x - R_major)^2 + z^2) > R_minor
        return NaN
    end
    return norm(getB_tokamak_coil(x, 0.0, z, a, b, ICoils, IPlasma))
end
Bmag = [get_Bmag_tokamak(x, z) for x in xs, z in zs]

fig = Figure(size = (700, 600), fontsize=20)
ax = Axis(fig[1, 1];
    xlabel="x (R)", ylabel="z", aspect=DataAspect(), title="Tokamak Coil Field (Poloidal)")

hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1e-9), colormap=:plasma)
Colorbar(fig[1, 2], hm, label="log10(|B|)")

str = evenstream(xs, zs, field_xz_tokamak)
streamlines!(ax, str; color = :white, linewidth = 1.0, with_arrows = true)

fig
```

Visualizing the 3D field lines:

```@example tokamak
# topology for coils: Major radius R_p = 2a + b, Minor radius R_coil = 2a
tor_mesh = get_tokamak_topology(R_major, R_minor)

fig = Figure(size = (800, 700), fontsize=20)
ax = Axis3(fig[1, 1]; title="3D Tokamak Coil Field Lines", aspect=:data)

wireframe!(ax, tor_mesh, color = (:blue, 0.1), linewidth = 0.5, transparency = true)

# Define the grid for field lines
xs_3d = range(-R_major - R_minor, R_major + R_minor, length=21)
ys_3d = range(-R_major - R_minor, R_major + R_minor, length=21)
zs_3d = range(-R_minor, R_minor, length=11)

function field_3d(p)
    x, y, z = p
    R = sqrt(x^2 + y^2)
    r_local = sqrt((R - R_major)^2 + z^2)
    if r_local > R_minor
        return SVector(0.0, 0.0, 0.0)
    end
    return getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma)
end

str = evenstream(xs_3d, ys_3d, zs_3d, field_3d; min_density=1.0, max_density=2.0)
streamlines!(ax, str; color = :gray, linewidth = 1.0)

fig
```

## Tokamak with q-profile

Reconstruct the magnetic field distribution from a safety factor (``q``) profile. The magnetic field components in local cylindrical coordinates ``(R, \zeta, z)`` relative to the major radius ``R_0`` and minor radius ``r`` are given by:
```math
\begin{aligned}
B_\zeta &= B_{\zeta 0} \frac{R_0}{R} \\
B_\theta &= \frac{r B_\zeta}{R_0 q(r/a)}
\end{aligned}
```
where ``B_{\zeta 0}`` is the toroidal field on axis, ``q(r/a)`` is the safety factor profile, and ``\theta`` is the poloidal angle.

```@example tokamak
# Define q-profile as a function of normalized radius r/a
q_profile(r_norm) = 1.1 + r_norm^2

# Minor radius a = 1.0 is already defined
const R0 = 3.0  # Major radius
const B0 = 2.0  # Toroidal field on axis

# Query point inside plasma
x, y, z = 3.5, 0.0, 0.0
B = getB_tokamak_profile(x, y, z, q_profile, a, R0, B0)
println("B at ($x, $y, $z): $B [T]")
```

Visualizing the q-profile field:

```@example tokamak
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
ax = Axis(fig[1, 1];
    xlabel="x (R)", ylabel="z", aspect=DataAspect(), title="Tokamak q-profile Field")

hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1e-9), colormap=:plasma)
Colorbar(fig[1, 2], hm, label="log10(|B|)")

str = evenstream(xs, zs, (x, z) -> field_xz_tokamak_q(x, z)[1], (x, z) -> field_xz_tokamak_q(x, z)[2])
streamlines!(ax, str; color = :white, linewidth = 1.5, with_arrows = true)
fig
```

Visualizing the 3D field lines with q-profile:

```@example tokamak
tor_mesh = get_tokamak_topology(R0, a)

fig = Figure(size = (800, 700), fontsize=20)
ax = Axis3(fig[1, 1]; title="3D Tokamak q-profile Field Lines", aspect=:data)

wireframe!(ax, tor_mesh, color = (:blue, 0.1), linewidth = 0.5, transparency = true)

# Define the grid for field lines
xs_3d = range(-R0 - a, R0 + a, length=21)
ys_3d = range(-R0 - a, R0 + a, length=21)
zs_3d = range(-a, a, length=11)

function field_3d(p)
    x, y, z = p
    R = sqrt(x^2 + y^2)
    r = sqrt((R - R0)^2 + z^2)
    if r > a
        return SVector(0.0, 0.0, 0.0)
    end
    return getB_tokamak_profile(x, y, z, q_profile, a, R0, B0)
end

str = evenstream(xs_3d, ys_3d, zs_3d, field_3d; min_density=1.0, max_density=2.0)
streamlines!(ax, str; color = :gray, linewidth = 1.0)

fig
```
