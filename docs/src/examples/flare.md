# [FLARE Experiment Reconnection Coils](@id flare_example)

The **Facility for Laboratory Reconnection Experiments (FLARE)** at the Princeton
Plasma Physics Laboratory (PPPL) is a user facility designed to study magnetic
reconnection in the high-Lundquist-number ($S \sim 10^5$), multi-island regime
bridging laboratory and space/astrophysical plasmas.

FLARE scales up the cylindrical flux-core reconnection concept established by the
Magnetic Reconnection Experiment (MRX):
- **Flux Cores**: Two toroidally symmetric flux cores positioned at major radius
  $R_{\text{pf}} = 0.75$ m and axial separation $z = \pm 0.40$ m. Each flux core
  contains internal **Poloidal Field (PF) coils**.
- **Equilibrium Field (EF) Coils**: External coils positioned at
  $R_{\text{ef}} \approx 1.30$ m and $z = \pm 0.70$ m that supply a vertical field
  $B_z$ to position the reconnection X-line and balance the plasma hoop force.
- **Plasma Current Sheet**: In driven reconnection, inductive flux transfer drives a
  midplane plasma current sheet ($I_p$) around $R \approx 0.75$ m, forming a macroscopic
  reconnection layer.
- **Guide Field**: An optional out-of-plane toroidal field $B_{\text{guide}}$ can be
  applied for guide-field / component reconnection studies.

The magnetic topology in the poloidal plane is characterized by the poloidal flux function:
```math
\psi(R, z) = R \, A_\phi(R, z)
```
where contours of $\psi$ represent poloidal magnetic field lines.

## Basic Setup

```@example flare
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

# Standard FLARE configuration
flare = FLARECoils(;
    R_pf = 0.75,
    z_pf = 0.40,
    I_pf = 45.0e3,  # 45 kA in flux core PF coils
    R_ef = 1.30,
    z_ef = 0.70,
    I_ef = 12.0e3,  # 12 kA in external EF coils
    R_p = 0.75,
    I_p = 0.0,
    B_guide = 0.0
)

# Field at the midplane X-line region
r_mid = SVector(0.75, 0.0, 0.0)
println("B at (0.75, 0, 0): ", flare(r_mid), " [T]")
println("Poloidal flux ψ: ", poloidal_flux(flare, r_mid), " [Wb]")
```

## Poloidal Flux Surfaces and X-Point Topology

In vacuum reconnection, the PF coils in the two flux cores produce an opposing radial
field component across the midplane $z = 0$, forming a magnetic null (X-point) where
$B_R = 0$ and $B_z = 0$.

```@example flare
Rs = range(0.3, 1.4, length=120)
Zs = range(-0.7, 0.7, length=140)

psi_vac = [poloidal_flux(flare, SVector(R, 0.0, Z)) for R in Rs, Z in Zs]

fig = Figure(size=(800, 650), fontsize=16)
ax = Axis(
    fig[1, 1],
    xlabel="R [m]",
    ylabel="Z [m]",
    aspect=DataAspect(),
    title="FLARE: Poloidal Flux Surfaces ψ(R, z) with X-Line Separators"
)

co = contour!(ax, Rs, Zs, psi_vac, levels=25, colormap=:viridis, linewidth=1.8)
Colorbar(fig[1, 2], co, label="ψ [Wb]")

# Mark flux cores (PF coils) and EF coils
scatter!(
    ax,
    [Point2f(0.75, 0.40), Point2f(0.75, -0.40)],
    color=:crimson,
    markersize=14,
    label="PF Flux Cores"
)
scatter!(
    ax,
    [Point2f(1.30, 0.70), Point2f(1.30, -0.70)],
    color=:dodgerblue,
    markersize=12,
    label="EF Coils"
)

axislegend(ax, position=:rt)
fig
```

## Radial X-Line Positioning with EF Coils

By varying the external equilibrium coil current $I_{\text{ef}}$, experimentalists at
FLARE can shift the radial position of the reconnection X-line along the midplane.

```@example flare
# Scan radial profile of Bz along the midplane for different I_ef
R_scan = range(0.4, 1.2, length=150)

fig2 = Figure(size=(800, 450), fontsize=16)
ax2 = Axis(
    fig2[1, 1],
    xlabel="R [m]",
    ylabel="B_z(R, 0) [T]",
    title="Midplane Vertical Magnetic Field vs. EF Coil Current"
)
hlines!(ax2, [0.0], color=:gray, linestyle=:dash)

ef_cases = [
    (0.0, "I_ef = 0 kA"),
    (12.0e3, "I_ef = 12 kA"),
    (24.0e3, "I_ef = 24 kA")
]
for (ief, label_str) in ef_cases
    fl = FLARECoils(;
        R_pf = 0.75, z_pf = 0.40, I_pf = 45.0e3,
        R_ef = 1.30, z_ef = 0.70, I_ef = ief
    )
    Bz_vals = [fl(SVector(R, 0.0, 0.0))[3] for R in R_scan]
    lines!(ax2, R_scan, Bz_vals, linewidth=2.2, label=label_str)
end

axislegend(ax2, position=:rb)
fig2
```

The zero-crossing $B_z(R, 0) = 0$ corresponds to the radial position of the X-line.
Increasing the EF coil current pushes the magnetic null inward or outward depending
on the polarity of $I_{\text{ef}}$.

## Reconnection Current Sheet and Guide Field

When plasma current $I_p$ is induced during a discharge, the X-point bifurcates into
a pair of Y-points bounding a finite-width current sheet. Furthermore, applying an
external toroidal field introduces a guide field $B_\phi$.

```@example flare
# FLARE with induced sheet current I_p = 10 kA and toroidal guide field B_g = 0.02 T
flare_reconn = FLARECoils(;
    R_pf = 0.75,
    z_pf = 0.40,
    I_pf = 45.0e3,
    R_ef = 1.30,
    z_ef = 0.70,
    I_ef = 12.0e3,
    R_p = 0.75,
    I_p = 10.0e3,
    B_guide = 0.02
)

# Total field magnitude along the midplane
B_mid_tot = [norm(flare_reconn(SVector(R, 0.0, 0.0))) for R in R_scan]

fig3 = Figure(size=(800, 450), fontsize=16)
ax3 = Axis(
    fig3[1, 1],
    xlabel="R [m]",
    ylabel="|B| [T]",
    title="Total Magnetic Field Magnitude on Midplane with Guide Field"
)
lines!(ax3, R_scan, B_mid_tot, color=:teal, linewidth=2.4, label="|B| (with guide field)")
axislegend(ax3, position=:rt)
fig3
```
