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

## Instrument Geometry

![FLARE coil geometry: toroidal PF flux cores and external EF coils surrounding
an illustrative midplane plasma current layer, with a matching R-Z section
and nominal coil coordinates.](../assets/flare_geometry.svg)

*Idealized geometry used in this demo, not an engineering drawing of the device.*
The left panel shows how the coils encircle the symmetry axis; the right panel
cuts through those rings at fixed toroidal angle to give the $(R, Z)$ view used
in the field plots below. The red rings represent the flux cores containing the
PF coils, and the larger blue rings represent the external EF coils. The orange
ring illustrates the optional midplane plasma current layer, not a hardware coil;
it is absent from the vacuum-field basic setup (`I_p = 0`). Coil thicknesses and
the plasma layer width are schematic; the vacuum vessel, supports, and diagnostics
are omitted.

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

- **Negative EF current** ($I_{\text{ef}} < 0$): Pulls the X-line inward toward
  the flux cores ($R_X \approx 0.88$ m).
- **Zero EF current** ($I_{\text{ef}} = 0$): Natural null position
  ($R_X \approx 0.97$ m).
- **Positive EF current** ($I_{\text{ef}} > 0$): Pushes the X-line outward toward the outer wall
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
condition is satisfied by construction, and the toroidal current density follows
from the Grad-Shafranov operator, $\mu_0 R J_\phi = -\Delta^*\psi$.

To test both the field and the current recovery, the reference field combines the
FLARE coils with a distributed 15 kA toroidal current sheet between the flux cores,

```math
J_\phi(R, Z) \propto \exp\left[-\left(\frac{R - R_c}{a}\right)^2\right]
\operatorname{sech}^2\left(\frac{Z}{\delta}\right),
```

with $R_c = 0.95$ m, $a = 0.15$ m and $\delta = 0.05$ m. The sheet is discretized
into coaxial loops with Gaussian cores of width $\sigma = 1$ cm. Each core adds the
exact two-dimensional correction to the loop field, so $\mathbf{B}$ is smooth and
$\nabla \times \mathbf{B} = \mu_0 \mathbf{J}$ holds to high accuracy.

```@example flare
struct GaussianCoreSheet
    loops::Vector{CurrentLoop{Float64}}
    σ::Float64
end

function GaussianCoreSheet(; Rc, a, δ, Ip, h=0.01, σ=0.01)
    # The offset keeps the cores away from the regular evaluation grids
    nodes = [(R + 0.00317, Z + 0.00213) for R in (Rc - 3a):h:(Rc + 3a)
                                         for Z in -4δ:h:4δ]
    w = [exp(-((R - Rc) / a)^2) * sech(Z / δ)^2 for (R, Z) in nodes]
    loops = [CurrentLoop(R, Ip * wk / sum(w), SVector(0.0, 0.0, Z),
                 SVector(0.0, 0.0, 1.0)) for ((R, Z), wk) in zip(nodes, w)]
    return GaussianCoreSheet(loops, σ)
end

function (s::GaussianCoreSheet)(p)
    R = hypot(p[1], p[2])
    B = zero(SVector{3, Float64})
    for loop in s.loops
        B += getB_loop(p, loop)
        dR, dZ = R - loop.radius, p[3] - loop.center[3]
        d2 = dR^2 + dZ^2
        c = -Magnetostatics.μ₀ * loop.current * exp(-d2 / (2s.σ^2)) / (2π * d2)
        B += SVector(c * dZ * p[1] / R, c * dZ * p[2] / R, -c * dR)
    end
    return B
end

function toroidal_current(s::GaussianCoreSheet, R, Z)
    return sum(s.loops) do loop
        d2 = (R - loop.radius)^2 + (Z - loop.center[3])^2
        loop.current * exp(-d2 / (2s.σ^2))
    end / (2π * s.σ^2)
end

sheet = GaussianCoreSheet(; Rc=0.95, a=0.15, δ=0.05, Ip=15.0e3)
reference(p) = flare(p) + sheet(p)
nothing # hide
```

We sample this field with magnetic probes and add small deterministic errors to
create synthetic measurements. In an experiment, replace these samples with
calibrated measurements and their uncertainties. A few modeling choices are
essential for an accurate reconstruction:
- **Exclude the sources.** The PF filaments at $(0.75, \pm 0.40)$ m are singular,
  so the reconstruction domain covers only the plasma region between the flux
  cores, where probes can actually be placed.
- **Represent the external coil field.** Virtual current filaments placed outside
  the domain (`filament_positions`) give exact current-free solutions inside it and
  capture the steep field near the flux cores. Their currents are fit to the data,
  so neither the coil currents nor their positions are assumed.
- **Represent the plasma current.** A bump in $\psi$ with compact support, such as
  the default `:wendland_c2` radial basis function, carries no net toroidal current,
  since $\oint \mathbf{B}\cdot d\mathbf{l} = 0$ around its support. We compare two
  remedies:
  1. *Gaussian kernel* (`kernel = :gaussian`) with a 6 cm length scale. Its wide
     tails can place the return current of each bump outside the domain, but every
     basis function touches every probe, so the linear system is dense.
  2. *Current elements* (`current_element_positions`) on a 5 cm grid, combined with
     the default sparse Wendland basis. Each element is a coaxial loop with a
     smooth core, so it carries net current and the sum of the fitted element
     currents estimates the plasma current.
- **Regularize the current.** Imposing $J_\phi = 0$ at collocation points with a
  loose uncertainty of $10^5$ A/m² (below the $\sim 10^6$ A/m² peak) selects the
  smallest current distribution consistent with the probes. A much tighter
  uncertainty would instead enforce a vacuum field.

```@example flare
using LinearSolve: QRFactorization

Rmin, Rmax = 0.45, 1.25
Zmin, Zmax = -0.28, 0.28
rz_matrix(points) = reduce(vcat, [[R Z] for (R, Z) in points])
grid_rz(h) = [(R, Z) for R in Rmin:h:Rmax for Z in Zmin:h:Zmax]

function synthetic_probes(probe_rz)
    B = reduce(vcat, [
        permutedims(collect(reference(SVector(R, 0.0, Z)))) for (R, Z) in probe_rz
    ])
    for i in axes(B, 1)
        B[i, 1] += 2.0e-4 * sin(i)
        B[i, 2] += 2.0e-4 * cos(i)
        B[i, 3] += 2.0e-4 * sin(2i)
    end
    return B
end

# Virtual filaments evenly spaced on a rectangle 0.1 m outside the domain
function rectangle_points(r0, r1, z0, z1, n)
    corners = [(r0, z0), (r1, z0), (r1, z1), (r0, z1), (r0, z0)]
    sides = [hypot(b[1] - a[1], b[2] - a[2])
             for (a, b) in zip(corners[1:end-1], corners[2:end])]
    ends = cumsum(sides)
    return map(range(0, ends[end]; length=n + 1)[1:end-1]) do s
        k = findfirst(>(s), ends)
        t = 1 - (ends[k] - s) / sides[k]
        a, b = corners[k], corners[k + 1]
        (a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]))
    end
end
filament_rz = rectangle_points(Rmin - 0.1, Rmax + 0.1, Zmin - 0.1, Zmax + 0.1, 48)
element_rz = grid_rz(0.05)
collocation_rz = grid_rz(0.015)

function reconstruct(probe_rz, method)
    basis = if method === :gaussian
        (; kernel=:gaussian, length_scale=0.06,
            basis_centers=rz_matrix(grid_rz(0.03)))
    else
        (; length_scale=0.12, basis_centers=rz_matrix(grid_rz(0.06)),
            current_element_positions=rz_matrix(element_rz))
    end
    probe_B = synthetic_probes(probe_rz)
    return reconstruct_axisymmetric_field(
        rz_matrix(probe_rz), probe_B;
        domain = (Rmin, Rmax, Zmin, Zmax),
        basis...,
        regularization = 1.0e2,
        uncertainties = fill(5.0e-4, size(probe_B)),
        filament_positions = rz_matrix(filament_rz),
        current_positions = rz_matrix(collocation_rz),
        current_density = zeros(length(collocation_rz)),
        current_uncertainties = fill(1.0e5, length(collocation_rz)),
        algorithm = QRFactorization()
    )
end
nothing # hide
```

The comparison helper evaluates the in-plane magnetic field
$\mathbf{B}_p = (B_R, B_Z)$ and the toroidal current density on a fine grid, prints
error metrics, and plots one column per field. The top row shows $|\mathbf{B}_p|$
with field lines, where the red dots mark the probe locations used in the fit. The
middle row shows the in-plane field error, or the probe layout for the reference.
The third row compares $J_\phi$, which is never measured directly; gray dots mark
the current elements. The bottom row shows radial and vertical cuts of $J_\phi$
through the sheet.

```@example flare
using UniformStreamlines

R_view = range(Rmin, Rmax; length=81)
Z_view = range(Zmin, Zmax; length=57)
dA = step(R_view) * step(Z_view)
inplane(field, R, Z) = (B = field(SVector(R, 0.0, Z)); SVector(B[1], B[3]))
toroidal_current(field::AxisymmetricFieldReconstruction, R, Z) =
    current_density(field, SVector(R, 0.0, Z))[2]
Bp_reference = [inplane(reference, R, Z) for R in R_view, Z in Z_view]
J_reference = [toroidal_current(sheet, R, Z) for R in R_view, Z in Z_view]
Bp_max = maximum(norm, Bp_reference)
J_max = maximum(J_reference) / 1e6

Bp_label = rich("|B", subscript("p"), "| [T]")
J_label = rich("J", subscript("ϕ"), " [MA/m²]")
error_label = rich("|ΔB", subscript("p"), "| / max|B", subscript("p"), "|")

function compare_reconstructions(probe_rz, fits)
    fig = Figure(size=(1600, 1900), fontsize=22)
    columns = [("Reference", "coils + current sheet", reference, sheet);
               [("Reconstruction", title, fit, fit) for (title, fit) in fits]]
    ax_R = Axis(fig[5, 1:2]; xlabel="R [m]", ylabel=J_label,
        title="Midplane profile (Z = 0)")
    ax_Z = Axis(fig[5, 3]; xlabel="Z [m]", ylabel=J_label,
        title="Vertical profile (R = 0.95 m)")
    probe_style = (color=:red, strokecolor=:black, strokewidth=0.8, markersize=9)
    profile_styles = [(:black, :solid), ("#0072B2", :dash), ("#D55E00", :dashdot)]
    header_colors = [(:gray60, 0.35); fill((:steelblue, 0.25), length(fits))]
    for (col, (kind, title, field, source)) in enumerate(columns)
        Bp = [inplane(field, R, Z) for R in R_view, Z in Z_view]
        J = [toroidal_current(source, R, Z) for R in R_view, Z in Z_view]
        Bp_error = norm.(Bp .- Bp_reference)
        if col > 1
            println(rpad(title, 34),
                "B_p RMS ", round(100 * norm(Bp_error) / norm(norm.(Bp_reference)),
                    digits=2), " %, ",
                "max ", round(100 * maximum(Bp_error) / Bp_max, digits=2), " %; ",
                "J_ϕ L2 ", round(100 * norm(J - J_reference) / norm(J_reference),
                    digits=1), " %; ",
                "I_p ", round(sum(J) * dA / 1e3, digits=2), " kA")
        end

        Box(fig[1, col]; color=header_colors[col], strokevisible=false)
        Label(fig[1, col], rich(rich(kind; font=:bold, fontsize=26), "\n", title);
            padding=(0, 0, 8, 8), justification=:center, tellwidth=false)

        ax = Axis(fig[2, col]; xlabel="R [m]", ylabel=col == 1 ? "Z [m]" : "",
            aspect=DataAspect(), title=rich("|B", subscript("p"), "| and field lines"))
        heatmap!(ax, R_view, Z_view, norm.(Bp);
            colormap=:cividis, colorrange=(0, Bp_max))
        str = evenstream(R_view, Z_view, (R, Z) -> inplane(field, R, Z);
            min_density=2.0, max_density=4.0)
        streamlines!(ax, str; color=:white, linewidth=1.2, with_arrows=true)
        scatter!(ax, first.(probe_rz), last.(probe_rz); probe_style...)

        ax = Axis(fig[3, col]; xlabel="R [m]", ylabel=col == 1 ? "Z [m]" : "",
            aspect=DataAspect(),
            title=col == 1 ? "Probe layout ($(length(probe_rz)) probes)" :
                rich("Field error ", error_label))
        heatmap!(ax, R_view, Z_view, Bp_error ./ Bp_max;
            colormap=:inferno, colorrange=(0, 0.05))
        scatter!(ax, first.(probe_rz), last.(probe_rz); color=:white, markersize=6)

        ax = Axis(fig[4, col]; xlabel="R [m]", ylabel=col == 1 ? "Z [m]" : "",
            aspect=DataAspect(), title=rich("Toroidal current J", subscript("ϕ")))
        heatmap!(ax, R_view, Z_view, J ./ 1e6;
            colormap=:balance, colorrange=(-J_max, J_max))
        if field isa AxisymmetricFieldReconstruction &&
                !isempty(field.element_currents)
            scatter!(ax, field.element_positions[:, 1], field.element_positions[:, 2];
                color=(:black, 0.3), markersize=4)
        end

        color, linestyle = profile_styles[col]
        label = col == 1 ? "Reference" : title
        lines!(ax_R, R_view, toroidal_current.(Ref(source), R_view, 0.0) ./ 1e6;
            color, linestyle, linewidth=3, label)
        lines!(ax_Z, Z_view, toroidal_current.(Ref(source), 0.95, Z_view) ./ 1e6;
            color, linestyle, linewidth=3)
    end
    Legend(fig[6, 1:3], ax_R; orientation=:horizontal, framevisible=false)
    Colorbar(fig[2, 4]; limits=(0, Bp_max), colormap=:cividis, label=Bp_label)
    Colorbar(fig[3, 4]; limits=(0, 0.05), colormap=:inferno, label="Relative error")
    Colorbar(fig[4, 4]; limits=(-J_max, J_max), colormap=:balance, label=J_label)
    colsize!(fig.layout, 1, Auto(1)); colsize!(fig.layout, 2, Auto(1))
    colsize!(fig.layout, 3, Auto(1))
    rowgap!(fig.layout, 1, 8)
    return fig
end
nothing # hide
```

### Uniform Probe Array

We first sample the field on a uniform $12 \times 11$ probe array.

```@example flare
uniform_rz = [(R, Z) for R in range(0.48, 1.22; length=12)
                      for Z in range(-0.25, 0.25; length=11)]
gaussian_fit = reconstruct(uniform_rz, :gaussian)
element_fit = reconstruct(uniform_rz, :elements)
println("Sum of fitted element currents: ",
    round(sum(element_fit.element_currents) / 1e3, digits=2), " kA")
fig4 = compare_reconstructions(uniform_rz, [
    "Gaussian RBF" => gaussian_fit,
    "Wendland RBF + current elements" => element_fit
])
```

Both methods recover the magnetic field to about a percent and the total current
to within a few tenths of a kiloampere. The largest field errors lie inside the
sheet, where $B_R$ reverses over a few centimeters between probe rows. The current
density is never measured, yet its location, shape and total are recovered from
the field alone; its peak is smoothed by the minimum-current regularization and
the finite probe spacing. The Gaussian kernel gives a slightly smoother current,
while current elements keep the default sparse Wendland basis and report the
plasma current directly through `element_currents`.

### Non-Uniform Probe Coverage

Probe access is rarely uniform in an experiment. Here, the upper half of the
domain keeps a dense $12 \times 8$ array, while the lower half has only three rows
of six probes. The reconstruction does not assume the up-down symmetry of the
reference field, so the poorly sampled half must be inferred from fewer data.

```@example flare
upper_rz = [(R, Z) for R in range(0.48, 1.22; length=12)
                    for Z in range(0.0, 0.25; length=8)]
lower_rz = [(R, Z) for R in range(0.48, 1.22; length=6)
                    for Z in range(-0.25, -0.25 / 3; length=3)]
nonuniform_rz = vcat(upper_rz, lower_rz)
fig5 = compare_reconstructions(nonuniform_rz, [
    "Gaussian RBF" => reconstruct(nonuniform_rz, :gaussian),
    "Wendland RBF + current elements" => reconstruct(nonuniform_rz, :elements)
])
```

The upper half is reconstructed as accurately as with the uniform array. The error
grows to several percent in the sparsely sampled lower half, concentrated in the
gap between the midplane and the first lower probe row, where the lower part of the
current sheet is not bracketed by measurements, and along the bottom edge between
the widely spaced probes. The vertical profile shows that the recovered current
sheet becomes asymmetric: it stays sharp on the well-sampled side but is smeared
toward the sparsely sampled side, and the total current drops by about 2%. Both
methods degrade in the same way, so the probe coverage, not the basis, limits the
accuracy here. When the configuration is known to be up-down symmetric, mirroring
the measurements, $B_R(R, -Z) = -B_R(R, Z)$ and $B_Z(R, -Z) = B_Z(R, Z)$, restores
the missing coverage.

$B_\phi$ is fit independently from its probe samples. Known exact boundary flux
values can be supplied with `boundary_positions` and `boundary_flux`, and measured
toroidal currents with `current_positions` and `current_density`. For a smooth
vacuum background without nearby coils, `vacuum_degree` adds low-order polynomial
solutions of $\Delta^*\psi = 0$ as an alternative to virtual filaments.
Reconstruction is restricted to the specified domain and does not assume zero
field at its boundary.
