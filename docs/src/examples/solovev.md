# [Solov'ev Grad-Shafranov Equilibrium](@id solovev_example)

The Solov'ev equilibrium is an exact analytical solution to the Grad-Shafranov equation:
```math
\Delta^* \psi = R \frac{\partial}{\partial R} \left( \frac{1}{R} \frac{\partial \psi}{\partial R} \right) + \frac{\partial^2 \psi}{\partial Z^2} = -\mu_0 R^2 p'(\psi) - F F'(\psi)
```
with linear pressure profile ``p'(\psi) = \text{const}`` and ``F F'(\psi) = 0``.

The poloidal flux function ``\psi(R, Z)`` is given by:
```math
\psi(R, Z) = \frac{\psi_0}{R_0^4} \left[ \frac{R^2 Z^2}{\kappa^2} + \frac{1}{4} (R^2 - R_0^2)^2 \right]
```
where ``R_0`` is the major radius of the magnetic axis, ``\kappa`` is the plasma elongation, and ``\psi_0 = \frac{R_0 B_0 \kappa}{2 q_0}`` sets the flux scale with on-axis safety factor ``q_0``.

The magnetic field components in cylindrical coordinates ``(R, \phi, Z)`` are:
```math
B_R = -\frac{1}{R} \frac{\partial \psi}{\partial Z}, \quad B_\phi = \frac{R_0 B_0}{R}, \quad B_Z = \frac{1}{R} \frac{\partial \psi}{\partial R}
```

## Basic Usage

```@example solovev
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

B0 = 2.0     # On-axis toroidal field [T]
R0 = 3.0     # Major radius [m]
kappa = 1.6  # Elongation
q0 = 1.05    # On-axis safety factor

eq = SolovevEquilibrium(B0, R0, kappa, q0)

# Field at the magnetic axis
r_axis = SVector(R0, 0.0, 0.0)
println("B on axis: ", eq(r_axis), " [T]")
println("Poloidal flux on axis: ", poloidal_flux(eq, r_axis), " [Wb]")
```

## Visualizing Flux Surfaces and Poloidal Field

```@example solovev
Rs = range(1.5, 4.5, length=61)
Zs = range(-2.0, 2.0, length=61)

psi_grid = [poloidal_flux(eq, SVector(R, 0.0, Z)) for R in Rs, Z in Zs]

fig = Figure(size=(750, 600), fontsize=18)
ax = Axis(fig[1, 1], xlabel="R [m]", ylabel="Z [m]", aspect=DataAspect(),
    title="Solov'ev Equilibrium: Poloidal Flux Surfaces ψ(R, Z)")

co = contour!(ax, Rs, Zs, psi_grid, levels=15, colormap=:plasma, linewidth=2)
Colorbar(fig[1, 2], colormap=:plasma, limits=extrema(psi_grid), label="ψ [Wb]")


# Mark magnetic axis
scatter!(ax, [Point2f(R0, 0.0)], color=:red, markersize=14, label="Magnetic Axis")
axislegend(ax, position=:rt)

fig
```
