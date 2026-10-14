# [Magnetic Geometry and Plasma Drifts](@id magnetic_geometry_example)

In plasma physics and magnetic confinement, single-particle dynamics are governed by guiding center motion in inhomogeneous magnetic fields. `Magnetostatics.jl` provides differential geometry tools to calculate field unit vectors, gradient of field strength, field-line curvature, and drift velocities.

## Guiding-Center Theory

1. **Unit magnetic vector**:
```math
\mathbf{b} = \frac{\mathbf{B}}{|\mathbf{B}|}
```
2. **Magnetic field gradient**:
```math
\nabla |\mathbf{B}|
```
3. **Magnetic curvature vector**:
```math
\boldsymbol{\kappa} = (\mathbf{b} \cdot \nabla) \mathbf{b}
```
4. **Guiding-center drift velocities**:
```math
\mathbf{v}_{\nabla B} = \frac{m v_\perp^2}{2 q B^2} (\mathbf{b} \times \nabla |\mathbf{B}|), \quad \mathbf{v}_\kappa = \frac{m v_\parallel^2}{q B} (\mathbf{b} \times \boldsymbol{\kappa})
```

## Example: Magnetic Mirror Geometry

```@example geometry
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

# Define a magnetic mirror
distance = 2.0
radius = 1.0
current = 1000.0
mirror = MagneticMirror(distance, radius, current)

# Evaluate field and geometric invariants at off-axis position
r = SVector(0.5, 0.0, 0.5)
B = mirror(r)
b = b_unit(mirror, r)
gradB = grad_b(mirror, r)
kappa = magnetic_curvature(mirror, r)

println("B: ", B)
println("b: ", b)
println("∇|B|: ", gradB)
println("Curvature κ: ", kappa)

# Calculate guiding center drift velocities for a proton (q=e, m=mp)
drifts = drift_velocities(mirror, r; q=1.602e-19, m=1.673e-27, v_perp=1.0e5, v_par=1.0e5)
println("Grad-B drift: ", drifts.v_gradB)
println("Curvature drift: ", drifts.v_curv)
println("Total magnetic drift: ", drifts.v_drift)

# Mirror ratio and loss-cone angle
r_throat = SVector(0.0, 0.0, 1.0)
r_center = SVector(0.0, 0.0, 0.0)
Rm = mirror_ratio(mirror, r_center, r_throat)
theta_loss = rad2deg(loss_cone_angle(mirror, r_center, r_throat))
println("Mirror ratio Rm: ", round(Rm, digits=3))
println("Loss cone angle θ_loss: ", round(theta_loss, digits=2), "°")
```

Visualizing $|B|$ and the curvature vector field in the $xz$-plane:

```@example geometry
xs = range(-1.5, 1.5, length=31)
zs = range(-1.5, 1.5, length=31)

Bmag = [field_magnitude(mirror, SVector(x, 0.0, z)) for x in xs, z in zs]

fig = Figure(size=(750, 600), fontsize=18)
ax = Axis(fig[1, 1], xlabel="x [m]", ylabel="z [m]", aspect=DataAspect(),
    title="Magnetic Mirror: |B| & Field Curvature κ")

hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1e-9), colormap=:viridis)
Colorbar(fig[1, 2], hm, label="log10(|B|)")

# Overlay curvature arrows at select points
pts = Point2f[]
kappas = Vec2f[]
for x in range(-1.2, 1.2, length=11)
    for z in range(-1.2, 1.2, length=11)
        r_pt = SVector(x, 0.0, z)
        k = magnetic_curvature(mirror, r_pt)
        if norm(k) > 1e-4
            push!(pts, Point2f(x, z))
            push!(kappas, Vec2f(k[1], k[3]))
        end
    end
end
arrows!(ax, pts, kappas; color=:white, arrowsize=8, lengthscale=0.15)

fig
```
