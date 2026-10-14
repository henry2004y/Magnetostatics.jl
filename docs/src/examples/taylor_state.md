# [Taylor State and Force-Free Equilibria](@id taylor_state_example)

A Woltjer-Taylor state is a minimum-energy, force-free relaxed magnetic configuration satisfying:
```math
\nabla \times \mathbf{B} = \alpha \mathbf{B}
```
where ``\alpha`` is a constant. In cylindrical geometry ``(r, \theta, z)``, the solution (Lundquist field) is given by Bessel functions of the first kind:
```math
B_z(r) = B_0 J_0(\alpha r), \quad B_\theta(r) = B_0 J_1(\alpha r), \quad B_r = 0
```
This configuration is the foundational model for:
- **Reversed Field Pinches (RFP)**: where ``B_z`` reverses sign at the first zero of ``J_0``, ``\alpha r \approx 2.4048``.
- **Spheromaks** and relaxed plasma states.
- **Solar coronal flux ropes**.

## Basic Usage

```@example taylor
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

B0 = 1.0     # Central axial field [T]
alpha = 3.0  # Force-free parameter [1/m]
state = TaylorState(B0, alpha)

# Evaluate on axis
r_center = SVector(0.0, 0.0, 0.0)
println("B on axis: ", state(r_center), " [T]")

# Check reversal radius
r_rev = 2.4048255577 / alpha
println("Reversal radius: ", round(r_rev, digits=4), " m")
println("B at reversal: ", state(SVector(r_rev, 0.0, 0.0)))
```

## Radial Profiles of Magnetic Field and Current Density

```@example taylor
rs = range(0.0, 1.2, length=100)

Bz_vals = [state(SVector(r, 0.0, 0.0))[3] for r in rs]
Bt_vals = [state(SVector(r, 0.0, 0.0))[2] for r in rs]

fig = Figure(size=(750, 500), fontsize=18)
ax = Axis(fig[1, 1], xlabel="Radius r [m]", ylabel="Magnetic Field [T]",
    title="Taylor State (RFP) Radial Profiles: B_z and B_θ")

lines!(ax, rs, Bz_vals, label="B_z (J₀)", linewidth=3, color=:blue)
lines!(ax, rs, Bt_vals, label="B_θ (J₁)", linewidth=3, color=:red)
hlines!(ax, [0.0], linestyle=:dash, color=:black)
vlines!(ax, [r_rev], linestyle=:dot, color=:green, label="Reversal Surface")

axislegend(ax, position=:rt)
fig
```
