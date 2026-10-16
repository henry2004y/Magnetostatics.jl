# [Spacecraft Grad-Shafranov Reconstruction](@id grad_shafranov_example)

This is **planar spacecraft GSR**, not the existing axisymmetric RBF flux fit.
It uses the equilibrium constraint

```math
\mathbf B=(A_y,-A_x,B_z(A)),\qquad
P_t(A)=p(A)+\frac{B_z(A)^2}{2\mu_0},\qquad
\nabla^2 A=-\mu_0 P_t'(A).
```

The toroidal Grad-Shafranov equation instead uses a cylindrical operator
``\Delta^\ast`` and radius-dependent source factors. Neither
[`reconstruct_axisymmetric_field`](@ref) nor [`reconstruct_field_3d`](@ref)
enforces pressure balance: those APIs reconstruct divergence-free fields from
spatial probes without assuming a magnetohydrostatic plasma.

Here a mock spacecraft crosses a static **Fadeev magnetic island** above its
center. We supply the invariant axis and crossing velocity and reconstruct a
2D strip from **121 magnetic-field and thermal-pressure samples along one line**.
The solver does not receive the island model or off-track truth.

This implementation follows the physical setup of
[Hau and Sonnerup (1999), DOI: 10.1029/1999JA900002](https://doi.org/10.1029/1999JA900002),
but uses truncated Chebyshev modes and second-order midpoint continuation,
not an exact reproduction of their numerical scheme.

## Mock spacecraft time series

Let

```math
f(x,y)=\cosh\left(\frac{y+y_0}{L}\right)
       -\epsilon\cos\left(\frac{x}{L}\right),\qquad
A=-B_0L\log f.
```

With equal longitudinal/transverse scales, this is an exact isotropic
equilibrium:

```math
p=p_0+\frac{B_0^2(1-\epsilon^2)}{2\mu_0 f^2},\qquad
B_z=B_g,\qquad
J_z=\frac{B_0(1-\epsilon^2)}{\mu_0L f^2}.
```

The island center is at ``(0,-y_0)``. We obtain its in-plane field from
[`FadeevIsland`](@ref), shifting by half a period and mapping its x-z plane
onto the reconstruction x-y plane. A constant guide field completes the 2.5D
geometry. Times are in seconds; the API uses SI units throughout.

```@example gsr
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie, Markdown

B0, L, epsilon = 20.0e-9, 1.0e6, 0.3
guide, offset, p0 = 5.0e-9, 0.2e6, 0.1e-9
island = FadeevIsland(B0, L, L, epsilon)
function truth(point)
    B = island(SVector(point[1] + π * L, 0.0, point[2] + offset))
    return SVector(B[1], B[3], guide)
end
denominator(x, y) = cosh((y + offset) / L) - epsilon * cos(x / L)
true_flux(x, y) = -B0 * L * log(denominator(x, y))
true_pressure(x, y) = p0 +
    B0^2 * (1 - epsilon^2) / (2Magnetostatics.μ₀ * denominator(x, y)^2)

times = collect(range(-40.0, 40.0; length = 121))
crossing_velocity = SVector(50.0e3, 0.0, 0.0)
invariant_axis = SVector(0.0, 0.0, 1.0)
track_x = crossing_velocity[1] .* times
positions = [SVector(x, 0.0, 0.0) for x in track_x]
observed = [truth(p)[k] for p in positions, k in 1:3]
observed .+= 0.02e-9 .* reshape(sin.(1:length(observed)), size(observed))
pressure = [true_pressure(x, 0.0) for x in track_x]
nothing # hide
```

The magnetic noise amplitude is **0.02 nT**, deterministic for reproducibility.
Pressure is noiseless here so that the pressure-branch comparison isolates
magnetic integration/noise error. For experimental data, include all species
contributing to isotropic thermal pressure; do not substitute density or
temperature alone.

`crossing_velocity` means **spacecraft velocity relative to the stationary
structure**, ``v_{\rm sc}-V_{\rm HT}``, not ``V_{\rm HT}`` itself. A stationary
spacecraft therefore uses the negative HT velocity. The solver projects this
velocity perpendicular to the supplied invariant axis and constructs a
right-handed frame. `origin` is the spacecraft position at the midpoint of
the time interval; its default is zero. Motion along the invariant axis is
irrelevant to this model.

## Reconstruct a narrow strip

```@example gsr
kwargs = (;
    crossing_velocity, invariant_axis,
    y = collect(range(-0.3L, 0.3L; length = 61)),
    spatial_degree = 10, pressure_degree = 6, axial_degree = 0,
)
reconstructed = reconstruct_grad_shafranov(times, observed, pressure; kwargs...)
reconstructed.diagnostics
```

The steps are:

1. Rotate the field into the supplied frame and integrate ``A(x,0)=-\int B_y\,dx``
   with the trapezoidal rule, setting the first sample's flux to zero.
2. Fit single-valued Chebyshev pressure and axial-field functions using normalized
   flux coordinates. The constant axial fit is appropriate for this mock only.
3. Fit along-track flux and ``A_y=B_x`` with a degree-10 spatial expansion.
4. Continue in both y directions with second-order midpoint steps. Project the
   nonlinear source onto the retained spatial modes at every step.
5. Interpolate flux with cubic Hermite segments in y and differentiate that same
   flux for the field and current. The representation is divergence-free, while
   the GS equation is satisfied only approximately. Current can have small jumps
   at interpolation knots because the flux is C1, not C2.

The source-function warning is intentional: some off-track flux values lie
outside the observed range. The diagnostic reports how much of the sampled
strip extrapolates; the warning triggers beyond 5% of the measured flux span.
This mock has an independently known smooth equilibrium, so we can check those
values. In real data, polynomial extrapolation can be unreliable. This initial
API does **not** fit exponential tails or infer an invariant axis automatically.

## Time series and pressure branches

Inbound and outbound samples of the same flux surfaces should agree in
``P_t(A)``. A small least-squares residual is a useful consistency check, but
does not independently establish stationarity, a correct axis, or force balance.

```@example gsr
track_flux = [flux_function(reconstructed, p) for p in positions]
measured_flux = zeros(length(times))
for i in 2:length(times)
    measured_flux[i] = measured_flux[i - 1] -
        (track_x[i] - track_x[i - 1]) * (observed[i, 2] + observed[i - 1, 2]) / 2
end
Pt = pressure + observed[:, 3].^2 / (2Magnetostatics.μ₀)
fitted_Pt = [
    thermal_pressure(reconstructed, p) + reconstructed(p)[3]^2 / (2Magnetostatics.μ₀)
        for p in positions
]

fig_track = Figure(; size = (1150, 750), fontsize = 18)
ax_B = Axis(fig_track[1, 1]; xlabel = "Time [s]", ylabel = "B [nT]")
for (k, label) in enumerate(("Bx", "By", "Bz"))
    lines!(ax_B, times, 1.0e9 .* observed[:, k]; label)
end
axislegend(ax_B; position = :rb)
ax_p = Axis(fig_track[2, 1]; xlabel = "Time [s]", ylabel = "Thermal pressure [nPa]")
lines!(ax_p, times, 1.0e9 .* pressure)
ax_branch = Axis(
    fig_track[1:2, 2]; xlabel = "A [mWb/m]", ylabel = "Pt [nPa]",
    title = "Two branches of Pt(A)"
)
midpoint = cld(length(times), 2)
scatter!(
    ax_branch, 1.0e3 .* measured_flux[1:midpoint], 1.0e9 .* Pt[1:midpoint];
    label = "Inbound", markersize = 8
)
scatter!(
    ax_branch, 1.0e3 .* measured_flux[midpoint:end], 1.0e9 .* Pt[midpoint:end];
    label = "Outbound", marker = :xcross, markersize = 9
)
order = sortperm(track_flux)
lines!(
    ax_branch, 1.0e3 .* track_flux[order], 1.0e9 .* fitted_Pt[order];
    label = "Fitted pressure function", color = :black
)
axislegend(ax_branch; position = :lt)
fig_track
```

## Reconstructed topology and held-out accuracy

Contours of ``A`` are the projected field lines, not full 3D guide-field lines.
The dashed line is the spacecraft trajectory; no measurements were supplied off
that line. Both topology panels use the **same flux levels and gauge**.

```@example gsr
xs = range(-1.7L, 1.7L; length = 101)
ys = range(-0.3L, 0.3L; length = 61)
gauge = true_flux(first(track_x), 0.0)
A_truth = [1.0e3 * (true_flux(x, y) - gauge) for x in xs, y in ys]
A_fit = [
    1.0e3 * flux_function(reconstructed, SVector(x, y, 0.0)) for x in xs, y in ys
]
error = [
    1.0e9 * norm(reconstructed(SVector(x, y, 0.0)) - truth(SVector(x, y, 0.0)))
        for x in xs, y in ys
]
base_levels = range(
    min(minimum(A_truth), minimum(A_fit)), max(maximum(A_truth), maximum(A_fit));
    length = 14
)
levels = sort(unique(vcat(collect(base_levels), maximum(A_truth) .- [0.02, 0.08, 0.2])))
fig_map = Figure(; size = (1450, 440), fontsize = 18)
for (column, data, title) in (
        (1, A_truth, "True flux topology"), (2, A_fit, "GSR from one crossing"),
    )
    ax = Axis(
        fig_map[1, column]; xlabel = "x [km]", ylabel = "y [km]", title
    )
    contour!(ax, xs ./ 1000, ys ./ 1000, data; levels, colormap = :viridis)
    hlines!(ax, [0.0]; color = :black, linestyle = :dash)
    scatter!(ax, [0.0], [-offset / 1000]; marker = :xcross, color = :black)
end
ax_error = Axis(
    fig_map[1, 3]; xlabel = "x [km]", ylabel = "y [km]", title = "Vector field error"
)
image_plot = heatmap!(ax_error, xs ./ 1000, ys ./ 1000, error; colormap = :magma)
Colorbar(fig_map[1, 4], image_plot; label = "|Bfit - B| [nT]")
fig_map
```

The display stretches y relative to x to make the narrow strip visible; axis
labels give the actual distances. Additional closely spaced flux levels reveal
small closed contours around the island center (black cross). A magnetic island
is a quasi-static equilibrium here, **not evidence of active reconnection**.

```@example gsr
heldout = [
    SVector(x, y, 0.0) for x in range(-1.7L, 1.7L; length = 17),
        y in (-0.23L, 0.17L, 0.27L)
]
expected = [truth(p)[k] for p in vec(heldout), k in 1:3]
predicted = [reconstructed(p)[k] for p in vec(heldout), k in 1:3]
relative_B_error = norm(predicted - expected) / norm(expected)
expected_J = [
    B0 * (1 - epsilon^2) /
        (Magnetostatics.μ₀ * L * denominator(p[1], p[2])^2) for p in heldout
]
predicted_J = [current_density(reconstructed, p)[3] for p in heldout]
relative_J_error = norm(predicted_J - expected_J) / norm(expected_J)
@assert relative_B_error < 0.03
@assert relative_J_error < 0.03
refined = reconstruct_grad_shafranov(
    times, observed, pressure; kwargs..., max_step = 2000.0
)
step_difference = norm(
    [refined(p)[k] for p in vec(heldout), k in 1:3] - predicted
) / norm(expected)
@assert step_difference < 0.005
println("Held-out vector field error: ", round(100relative_B_error; digits = 3), "%")
println("Held-out axial current error: ", round(100relative_J_error; digits = 3), "%")
println("Change on reducing the internal y step: ", step_difference)
println("GS residual [T/m]: ", reconstructed.diagnostics.gs_residual_rms)
```

For this mock the field/current errors are approximately **0.10% / 0.75%**.
These are synthetic validation errors, not uncertainty estimates. Independent
tests also check a polynomial equilibrium with varying ``B_z(A)``, arbitrary
frame rotation, exact pressure balance, and analytic current against the curl.

## Regularization, strip width and cost

An elliptic Cauchy problem is ill-posed: along-track noise at wavenumber
``k`` can grow like ``\exp(k|y|)``. Retaining more spatial modes may fit the
track better while making off-track reconstruction worse. Reducing the
integration step alone does not cure this instability.

```@example gsr
sensitivity = map((6, 10, 14)) do degree
    settings = merge(kwargs, (;
        spatial_degree = degree, y = collect(range(-0.6L, 0.6L; length = 121)),
    ))
    fit = reconstruct_grad_shafranov(times, observed, pressure; settings...)
    distances = (0.1L, 0.3L, 0.6L)
    errors = map(distances) do distance
        points = [SVector(x, distance, 0.0) for x in range(-1.7L, 1.7L; length = 31)]
        reference = [truth(p)[k] for p in points, k in 1:3]
        values = [fit(p)[k] for p in points, k in 1:3]
        norm(values - reference) / norm(reference)
    end
    timings = [
        @timed reconstruct_grad_shafranov(times, observed, pressure; settings...)
            for _ in 1:3
    ]
    fastest = timings[argmin([t.time for t in timings])]
    (; degree, errors, milliseconds = 1000fastest.time, bytes = fastest.bytes)
end
rows = [
    "| $(r.degree) | " *
    join(["$(round(100e; digits = 3))%" for e in r.errors], " | ") *
    " | $(round(r.milliseconds; digits = 1)) | $(r.bytes) |" for r in sensitivity
]
Markdown.parse(join(vcat([
    "| Spatial degree | y = 0.1L | y = 0.3L | y = 0.6L | Warm fit [ms] | Allocated bytes |",
    "|--:|--:|--:|--:|--:|--:|",
], rows), "\n"))
```

Warm timings exclude compilation but depend on hardware/load. The wider strip
also extrapolates the fitted pressure function farther: noise amplification and
source extrapolation both contribute to the error. The 0.2-track-length warning
is only a heuristic, **not a guaranteed safe width**. Choose spatial degree,
strip width and pressure-function complexity jointly and report sensitivity.

For this noise realization, degree 14 improves the error at `y = 0.1L` relative
to degree 10, but is worse at `y = 0.6L` (about 1.63% versus 0.87%). This is why
the best-looking track fit need not be the best continuation. The dedicated
`SUITE["grad_shafranov"]` benchmark group measures fitting and field, current,
flux and pressure evaluation separately.

This initial model deliberately omits HT/MVA/axis optimization, adaptive
regularization, exponential pressure tails, field-aligned-flow corrections,
pressure anisotropy and Hall/electron physics. It is intended for controlled
magnetohydrostatic reconstruction experiments, not automated interpretation of
active diffusion regions or arbitrary spacecraft intervals.
