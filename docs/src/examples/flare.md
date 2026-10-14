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

## In-Plane Field Lines and Reconnecting Field B_R

In the poloidal $(R, Z)$ plane, the reconnecting magnetic field is primarily the radial
component $B_R(R, Z)$, which has opposite signs across the midplane $z = 0$.
The poloidal magnetic flux function $\psi(R, Z) = R A_\phi(R, Z)$ defines the in-plane
magnetic field lines as its contours, with the separatrix passing directly through the
reconnection X-point null.

```@example flare
# Helper to locate the midplane X-point null
function find_xpoint(fl; Rmin=0.5, Rmax=1.3, N=2000)
    Rs = range(Rmin, Rmax, length=N)
    Bzs = [fl(SVector(R, 0.0, 0.0))[3] for R in Rs]
    for i in 1:(N-1)
        if Bzs[i] * Bzs[i+1] <= 0
            return Rs[i] - Bzs[i] * (Rs[i+1] - Rs[i]) / (Bzs[i+1] - Bzs[i])
        end
    end
    return NaN
end

Rs = range(0.35, 1.35, length=140)
Zs = range(-0.65, 0.65, length=140)

BR_grid = zeros(length(Rs), length(Zs))
psi_grid = zeros(length(Rs), length(Zs))

for (j, Z) in enumerate(Zs), (i, R) in enumerate(Rs)
    r = SVector(R, 0.0, Z)
    B = flare(r)
    BR_grid[i, j] = B[1]
    psi_grid[i, j] = poloidal_flux(flare, r)
end

Rx = find_xpoint(flare)
psi_X = poloidal_flux(flare, SVector(Rx, 0.0, 0.0))
b_max = 0.045

fig1 = Figure(size=(850, 720), fontsize=15)
ax1 = Axis(
    fig1[1, 1],
    xlabel="Major Radius R [m]",
    ylabel="Height Z [m]",
    aspect=DataAspect(),
    title="FLARE: In-Plane Magnetic Field Lines with B_R Reconnecting Field"
)

# Colormesh / heatmap of reconnecting radial field B_R
hm1 = heatmap!(
    ax1, Rs, Zs, BR_grid,
    colormap=:balance,
    colorrange=(-b_max, b_max)
)

# Poloidal magnetic field lines
contour!(
    ax1, Rs, Zs, psi_grid,
    levels=32,
    color=(:black, 0.75),
    linewidth=1.2
)

# Separatrix passing through the X-point
contour!(
    ax1, Rs, Zs, psi_grid,
    levels=[psi_X],
    color=:crimson,
    linewidth=2.8,
    linestyle=:dash
)

# Coils and X-point
p_pf = scatter!(
    ax1,
    [Point2f(0.75, 0.40), Point2f(0.75, -0.40)],
    color=:crimson,
    strokecolor=:black,
    strokewidth=1.5,
    markersize=16
)
p_ef = scatter!(
    ax1,
    [Point2f(1.30, 0.70), Point2f(1.30, -0.70)],
    color=:deepskyblue,
    strokecolor=:black,
    strokewidth=1.5,
    markersize=14
)
p_xp = scatter!(
    ax1,
    [Point2f(Rx, 0.0)],
    marker=:xcross,
    color=:gold,
    strokecolor=:black,
    strokewidth=2.5,
    markersize=22
)

# Top banner legend completely avoiding any coil overlap
Legend(
    fig1[0, 1],
    [p_pf, p_ef, p_xp, LineElement(color=:crimson, linestyle=:dash, linewidth=2.8)],
    ["PF Flux Cores", "EF Coils", "X-Point (R=$(round(Rx, digits=2)) m)", "Separatrix"],
    orientation = :horizontal,
    framevisible = false,
    tellwidth = false
)

Colorbar(
    fig1[1, 2], hm1,
    label="Radial Reconnecting Field B_R [T]"
)

fig1
```

## Radial X-Line Steering with EF Coils

By varying the external equilibrium coil current $I_{\text{ef}}$, experimentalists at
FLARE can shift the radial position of the reconnection X-line and control the private
flux volume:
- $I_{\text{ef}} < 0$: Pulls the X-line inward toward the flux cores ($R_X \approx 0.88$ m).
- $I_{\text{ef}} = 0$: Natural null position ($R_X \approx 0.97$ m).
- $I_{\text{ef}} > 0$: Pushes the X-line outward toward the outer wall
  ($R_X \approx 1.11$ m).

```@example flare
cases = [
    (-10.0e3, "Inward Shift (I_ef = -10 kA)"),
    (0.0, "Natural Null (I_ef = 0 kA)"),
    (10.0e3, "Outward Shift (I_ef = +10 kA)")
]

fig2 = Figure(size=(1250, 520), fontsize=14)

p_pf2, p_ef2, p_xp2 = nothing, nothing, nothing

for (idx, (ief, title_str)) in enumerate(cases)
    fl_case = FLARECoils(;
        R_pf = 0.75, z_pf = 0.40, I_pf = 45.0e3,
        R_ef = 1.30, z_ef = 0.70, I_ef = ief
    )
    psi_case = [poloidal_flux(fl_case, SVector(R, 0.0, Z)) for R in Rs, Z in Zs]
    rx_case = find_xpoint(fl_case)
    psi_x_case = poloidal_flux(fl_case, SVector(rx_case, 0.0, 0.0))

    ax = Axis(
        fig2[1, idx],
        xlabel="R [m]",
        ylabel=idx == 1 ? "Z [m]" : "",
        aspect=DataAspect(),
        title="$(title_str)\nX-Point at R = $(round(rx_case, digits=2)) m"
    )
    contour!(ax, Rs, Zs, psi_case, levels=25, colormap=:viridis, linewidth=1.4)
    contour!(
        ax, Rs, Zs, psi_case,
        levels=[psi_x_case],
        color=:red,
        linewidth=2.4,
        linestyle=:dash
    )

    global p_pf2 = scatter!(
        ax, [Point2f(0.75, 0.40), Point2f(0.75, -0.40)],
        color=:crimson, markersize=14
    )
    global p_ef2 = scatter!(
        ax, [Point2f(1.30, 0.70), Point2f(1.30, -0.70)],
        color=:deepskyblue, markersize=12
    )
    global p_xp2 = scatter!(
        ax, [Point2f(rx_case, 0.0)],
        marker=:xcross, color=:red, markersize=18, strokewidth=2
    )
end

Legend(
    fig2[0, 1:3],
    [p_pf2, p_ef2, p_xp2, LineElement(color=:red, linestyle=:dash, linewidth=2.4)],
    ["PF Coils", "EF Coils", "X-Point Null", "Separatrix"],
    orientation = :horizontal,
    framevisible = false,
    tellwidth = false
)

fig2
```

## Current Sheet Bifurcation & Y-Point Formation

When the reconnection drive induces a plasma current sheet ($I_p$) along the midplane,
the isolated vacuum X-point splits into a macroscopic current layer bounded by two
Y-points, squeezing the flux surfaces and creating a central magnetic island.

```@example flare
fl_vac2 = FLARECoils(;
    R_pf = 0.75, z_pf = 0.40, I_pf = 45.0e3,
    R_ef = 1.30, z_ef = 0.70, I_ef = 10.0e3,
    I_p = 0.0
)
fl_sheet2 = FLARECoils(;
    R_pf = 0.75, z_pf = 0.40, I_pf = 45.0e3,
    R_ef = 1.30, z_ef = 0.70, I_ef = 10.0e3,
    R_p = 0.78, I_p = 15.0e3
)

psi_vac2 = [poloidal_flux(fl_vac2, SVector(R, 0.0, Z)) for R in Rs, Z in Zs]
psi_sheet2 = [poloidal_flux(fl_sheet2, SVector(R, 0.0, Z)) for R in Rs, Z in Zs]

fig3 = Figure(size=(1100, 540), fontsize=15)

ax3a = Axis(
    fig3[1, 1],
    xlabel="R [m]",
    ylabel="Z [m]",
    aspect=DataAspect(),
    title="Vacuum Magnetic Topology (I_p = 0)\nSingle Classical X-Point Null"
)
contour!(ax3a, Rs, Zs, psi_vac2, levels=25, colormap=:magma, linewidth=1.5)
p_pf3 = scatter!(
    ax3a, [Point2f(0.75, 0.40), Point2f(0.75, -0.40)],
    color=:crimson, markersize=14
)
p_ef3 = scatter!(
    ax3a, [Point2f(1.30, 0.70), Point2f(1.30, -0.70)],
    color=:deepskyblue, markersize=12
)
p_xp3 = scatter!(
    ax3a, [Point2f(Rx, 0.0)],
    marker=:xcross, color=:gold,
    strokecolor=:black, strokewidth=2, markersize=20
)

ax3b = Axis(
    fig3[1, 2],
    xlabel="R [m]",
    ylabel="Z [m]",
    aspect=DataAspect(),
    title="Driven Reconnection (I_p = 15 kA)\nInduced Current Layer & Y-Point Formation"
)
contour!(ax3b, Rs, Zs, psi_sheet2, levels=25, colormap=:magma, linewidth=1.5)
scatter!(
    ax3b, [Point2f(0.75, 0.40), Point2f(0.75, -0.40)],
    color=:crimson, markersize=14
)
scatter!(
    ax3b, [Point2f(1.30, 0.70), Point2f(1.30, -0.70)],
    color=:deepskyblue, markersize=12
)
p_curr = scatter!(
    ax3b, [Point2f(0.78, 0.0)],
    marker=:star5, color=:gold,
    strokecolor=:black, strokewidth=1.5, markersize=20
)

Legend(
    fig3[0, 1:2],
    [p_pf3, p_ef3, p_xp3, p_curr],
    ["PF Coils", "EF Coils", "X-Point Null", "Reconnection Current Sheet (I_p)"],
    orientation = :horizontal,
    framevisible = false,
    tellwidth = false
)

fig3
```

## Reconstructing from Sparse Probe Measurements

An axisymmetric field can be reconstructed from sparse measurements of the
cylindrical components $(B_R, B_\phi, B_Z)$ in the poloidal plane. The poloidal
field is represented by the flux function $\psi(R,Z)$, so its divergence-free
condition is satisfied by construction. Here, we sample the existing FLARE coil
model at a small set of locations to create synthetic probe data. In an experiment,
replace these samples with calibrated measurements and their uncertainties.

The reconstruction also uses flux values around the edge of its bounded domain.
These are generated from the known reference field here; experimental edge values
should come from boundary diagnostics or an appropriate wall model. Boundary flux
alone does not determine the interior field, so it must be combined with probe data
and regularization. Because this finite basis cannot exactly represent the coil
model's edge flux everywhere, the example treats these edge values as uncertain.

```@example flare
probe_rz = [(R, Z) for R in range(0.48, 1.22; length=6)
                    for Z in range(-0.52, 0.52; length=5)]
probe_positions = reduce(vcat, [reshape(collect(p), 1, 2) for p in probe_rz])
probe_B = reduce(vcat, [
    reshape(collect(flare(SVector(R, 0.0, Z))), 1, 3) for (R, Z) in probe_rz
])
for i in axes(probe_B, 1)
    probe_B[i, 1] += 2.0e-4 * sin(i)
    probe_B[i, 2] += 2.0e-4 * cos(i)
    probe_B[i, 3] += 2.0e-4 * sin(2i)
end
probe_uncertainties = fill(5.0e-4, size(probe_B))

Rmin, Rmax = 0.4, 1.3
Zmin, Zmax = -0.6, 0.6
boundary_rz = vcat(
    [(Rmin, Z) for Z in range(Zmin, Zmax; length=9)],
    [(Rmax, Z) for Z in range(Zmin, Zmax; length=9)],
    [(R, Zmin) for R in range(Rmin, Rmax; length=9)],
    [(R, Zmax) for R in range(Rmin, Rmax; length=9)]
)
boundary_positions = reduce(vcat, [reshape(collect(p), 1, 2) for p in boundary_rz])
boundary_flux = [
    poloidal_flux(flare, SVector(R, 0.0, Z)) for (R, Z) in boundary_rz
]

reconstructed = reconstruct_axisymmetric_field(
    probe_positions, probe_B;
    domain = (Rmin, Rmax, Zmin, Zmax),
    length_scale = 0.3,
    regularization = 1.0e-10,
    uncertainties = probe_uncertainties,
    boundary_positions,
    boundary_flux,
    boundary_uncertainties = fill(2.0e-5, length(boundary_flux))
)
```

Compare the reconstructed flux contours with the reference field. The dots mark
the sparse probe locations used in the fit.

```@example flare
R_view = range(Rmin, Rmax; length=100)
Z_view = range(Zmin, Zmax; length=100)
psi_reference = [
    poloidal_flux(flare, SVector(R, 0.0, Z)) for R in R_view, Z in Z_view
]
psi_reconstructed = [
    poloidal_flux(reconstructed, SVector(R, 0.0, Z)) for R in R_view, Z in Z_view
]

fig4 = Figure(size=(1000, 440), fontsize=14)
ax4a = Axis(fig4[1, 1], xlabel="R [m]", ylabel="Z [m]",
    aspect=DataAspect(), title="Reference flux")
ax4b = Axis(fig4[1, 2], xlabel="R [m]", ylabel="Z [m]",
    aspect=DataAspect(), title="Sparse-probe reconstruction")
contour!(ax4a, R_view, Z_view, psi_reference, levels=24, colormap=:viridis)
contour!(ax4b, R_view, Z_view, psi_reconstructed, levels=24, colormap=:viridis)
scatter!(ax4b, first.(probe_rz), last.(probe_rz), color=:red, markersize=7)
fig4
```

`B_\phi` is fit independently from its probe samples. Known toroidal current
measurements can also be supplied with `current_positions` and `current_density`;
the fit uses the Grad-Shafranov relation for the poloidal flux. Reconstruction is
restricted to the specified domain and does not assume zero field at its boundary.
