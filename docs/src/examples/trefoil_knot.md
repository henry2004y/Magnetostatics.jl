# [Trefoil Knot](@id trefoil_knot_example)

A torus knot is a knot that lies on the surface of a torus. It is characterized by two coprime integers $p$ and $q$, where $p$ is the number of times the knot winds around the symmetry axis, and $q$ is the number of times it winds through the "hole" of the torus.

The parametric representation of a $(p, q)$ torus knot on a torus with major radius $R$ and minor radius $r$ is given by:
```math
\begin{aligned}
x(\phi) &= (R + r \cos(q\phi)) \cos(p\phi) \\
y(\phi) &= (R + r \cos(q\phi)) \sin(p\phi) \\
z(\phi) &= r \sin(q\phi)
\end{aligned}
```
where $\phi \in [0, 2\pi)$.

The trefoil knot is the simplest non-trivial knot, and it can be represented as a $(2, 3)$ torus knot.

In `Magnetostatics.jl`, the magnetic field of a trefoil knot is calculated by discretizing the knot into a series of straight wire segments and applying the Biot-Savart law to each segment.

```@example trefoil
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie, UniformStreamlines

# Define Trefoil Knot parameters
R, r = 1.0, 0.3
current = 10.0 # [A]
knot = TrefoilKnot(R, r, current)

# Discretize the knot into a wire with 300 segments
n_segments = 300
wire = Wire(knot, n_segments)

# Set up the grid for magnetic field calculation
xs = range(-2.0, 2.0, length=21)
ys = range(-2.0, 2.0, length=21)
zs = range(-1.0, 1.0, length=11)

solver = BiotSavart()

# Define the magnetic field function for streamlines
field(p) = solver(wire, p)

# Visualization
fig = Figure(size=(800, 700), fontsize=18)
ax = Axis3(fig[1,1], title="Magnetic Field of a Trefoil Knot",
    xlabel="x", ylabel="y", zlabel="z",
    aspect=:data)

# Plot the current-carrying wire
lines!(ax, wire.points, color=:red, linewidth=4, label="Wire")

# Plot 3D field lines
str = evenstream(xs, ys, zs, field; min_density=1, max_density=2)
streamlines!(ax, str; color = :gray, linewidth = 1.5)

fig
```
