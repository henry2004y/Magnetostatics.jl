# [3D Magnetic Field Reconstruction](@id reconstruction3d_example)

Reconstruct a continuous magnetic field inside a box from **36 vector probes**
and **54 normal-boundary samples**. A tilted, displaced external coil provides
the synthetic truth. The reconstruction receives measurements, not the coil
geometry or its analytical field.

The fitted representation is

```math
\mathbf B(\mathbf x)=\mathbf b_0+
\nabla\times\left[\ell\sum_j\phi_j(\mathbf x)\mathbf a_j\right],
```

so the field is divergence-free by construction. Sparse observations still do
not uniquely determine it: the basis and regularization select a smooth solution.
Boundary measurements below are uncertain observations, not exact wall conditions.

For other source terms and physical boundaries, see
[Reconstruction Across Sources and Boundaries](@ref reconstruction_cases_example).
The [inverse-method comparison](@ref inverse_methods) explains how this
general 3D model differs from axisymmetric flux reconstruction.

## Selected magnetic field

```@example reconstruction3d
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie, FieldTracer

domain = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
coil = CurrentLoop(
    1.0, 1.0e4, SVector(0.2, -0.1, 2.0), normalize(SVector(0.3, 0.2, 1.0))
)
truth = CurrentLoopAnalytic(coil)
```

The source is a **1 m radius loop carrying 10 kA**, centered at
``(0.2,-0.1,2.0)`` m with a tilted normal. The left panel locates it relative to
the 1 m reconstruction box; the right panel zooms into that box. Arrows point
along ``\mathbf B`` and use a common length scale, while their colors show
``|\mathbf B|`` in mT. The coil lies outside the measurement volume.

```@example reconstruction3d
box_edges = Point3f[]
for axis in 1:3, a in (-0.5, 0.5), b in (-0.5, 0.5)
    other_axes = filter(!=(axis), 1:3)
    for side in (-0.5, 0.5)
        point = zeros(3)
        point[axis] = side
        point[other_axes] .= (a, b)
        push!(box_edges, Point3f(point))
    end
end
coil_wire = discretize_loop(coil, 150)
field_points = vec([
    SVector(x, y, z) for x in (-0.3, 0.0, 0.3),
        y in (-0.3, 0.0, 0.3), z in (-0.3, 0.0, 0.3)
])
field_vectors = truth.(field_points)
field_strength = 1.0e3 .* norm.(field_vectors)
arrow_scale = 0.20 / maximum(norm.(field_vectors))

fig_setup = Figure(; size = (1150, 550), fontsize = 22)
source_axis = Axis3(
    fig_setup[1, 1]; title = "External tilted coil",
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
    aspect = :data, azimuth = 1.15pi, elevation = 0.18pi,
    xticks = [-1, 0, 1], yticks = [-1, 0, 1]
)
lines!(
    source_axis, Point3f.(coil_wire.points); color = :darkorange,
    linewidth = 4, label = "Coil: 10 kA"
)
linesegments!(
    source_axis, box_edges; color = :steelblue, linewidth = 2,
    label = "Reconstruction box"
)
Legend(fig_setup[2, 1], source_axis; orientation = :horizontal, framevisible = false)

field_axis = Axis3(
    fig_setup[1, 2]; title = "Analytical field inside the box",
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
    aspect = :data, azimuth = 1.15pi, elevation = 0.25pi
)
linesegments!(field_axis, box_edges; color = (:gray, 0.35), linewidth = 1)
field_arrows = arrows3d!(
    field_axis, Point3f.(field_points), Vec3f.(arrow_scale .* field_vectors);
    color = field_strength, colormap = :viridis,
    colorrange = extrema(field_strength), markerscale = 1,
    minshaftlength = 0, tiplength = 0.04, tipradius = 0.022, shaftradius = 0.004
)
limits!(field_axis, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
Colorbar(fig_setup[1, 3], field_arrows; label = "|B| [mT]", height = 300)
Label(fig_setup[2, 2], "Arrow direction: +B; length proportional to |B|")
fig_setup
```

Only this synthetic demonstration knows the source and analytical field. The
inverse fit below receives probe and boundary observations, not the coil model.

## Synthetic measurements

```@example reconstruction3d
points_to_matrix(points) =
    reduce(vcat, [reshape(collect(p), 1, 3) for p in points])
sample_field(field, positions) =
    points_to_matrix([field(SVector{3}(p)) for p in eachrow(positions)])

function box_surface(coordinates)
    points = SVector{3, Float64}[]
    normals = SVector{3, Float64}[]
    for axis in 1:3, side in (-0.5, 0.5), a in coordinates, b in coordinates
        other_axes = filter(!=(axis), 1:3)
        point = zeros(3)
        point[axis] = side
        point[other_axes] .= (a, b)
        push!(points, SVector{3}(point))
        push!(normals, SVector{3}(ntuple(k -> k == axis ? sign(side) : 0.0, 3)))
    end
    return points, normals
end

probe_positions = points_to_matrix(Iterators.product(
    (-0.45, -0.15, 0.15, 0.45), (-0.4, 0.0, 0.4), (-0.4, 0.0, 0.4)
))
probe_sigma = 1.0e-5
probe_B = sample_field(truth, probe_positions)
probe_B .+= probe_sigma .* reshape(sin.(1:length(probe_B)), size(probe_B))
probe_B[1, 2] = NaN
probe_B[12, 3] = NaN

boundary_points, boundary_directions = box_surface((-0.35, 0.0, 0.35))
boundary_positions = points_to_matrix(boundary_points)
boundary_normals = points_to_matrix(boundary_directions)
boundary_sigma = 2.0e-5
boundary_values = [
    dot(n, truth(p)) + boundary_sigma * cos(i)
        for (i, (p, n)) in enumerate(zip(boundary_points, boundary_directions))
]

println("Vector probes: ", size(probe_positions, 1))
println("Available probe components: ", count(isfinite, probe_B))
println("Normal-boundary samples: ", length(boundary_values))
```

Noise is deterministic to make the figures reproducible; its amplitude is set
by the supplied uncertainty scale. Two missing components demonstrate partial
probe measurements. In an experiment, use calibrated noise estimates and actual
boundary diagnostics rather than values generated from a known model.

## Fit the field

```@example reconstruction3d
basis_centers = points_to_matrix(Iterators.product(
    ntuple(_ -> range(-0.5, 0.5; length = 4), 3)...
))
reconstructed = reconstruct_field_3d(
    probe_positions, probe_B;
    domain, length_scale = 1.2, basis_centers, kernel = :gaussian,
    uncertainties = fill(probe_sigma, size(probe_B)),
    boundary_positions, boundary_normals, boundary_values,
    boundary_uncertainties = fill(boundary_sigma, length(boundary_values)),
    regularization = 1.0e-8, regularization_grid = (7, 7, 7),
    abstol = 1.0e-12, reltol = 1.0e-12
)
reconstructed.diagnostics
```

The global Gaussian basis is appropriate for this smooth external-coil field.
The default compact Wendland basis is useful for localized structure, but needs
support throughout the region of interest. Regularization penalizes correction
field and curl energy; it does **not** impose zero current exactly.

## Three-dimensional comparison

The left panel overlays analytical (solid blue) and reconstructed (dashed orange)
field lines traced from the **same six seeds**. Arrows indicate the direction of
the analytical field; gray markers locate the probes. These are magnetic field
lines, not particle trajectories. Nearly overlapping curves indicate agreement,
so a second panel explicitly shows the vector error at held-out volume locations
in microtesla. Both panels use the same camera and box limits.

FieldTracer.jl traces ``d\mathbf x/ds = \mathbf B/|\mathbf B|`` in both directions
with native fixed-step RK4. Its domain predicate is checked before every field
evaluation, including intermediate stages, so the reconstruction is never evaluated
outside its box. A step crossing the boundary is discarded: endpoints remain
inside, within one step of the boundary, rather than being clipped onto it.
We check termination explicitly so invalid fields or exhausted step limits cannot
silently produce incomplete curves.

```@example reconstruction3d
function trace_field_line(field, seed, box; ds = 0.01)
    isoutside(p, s) = any(p[k] < box[2k - 1] || p[k] > box[2k] for k in 1:3)
    result = trace(
        (p, s) -> field(p), collect(seed); ds, max_length = 4.0, isoutside
    )
    for branch in (result.backward, result.forward)
        branch.status == :outside || error("Unexpected tracing status: $(branch.status)")
    end
    return vcat(reverse(result.backward.points), result.forward.points[2:end])
end

seeds = vec([SVector(x, y, 0.0) for x in (-0.27, 0.0, 0.27), y in (-0.22, 0.22)])
true_paths = [trace_field_line(truth, seed, domain) for seed in seeds]
fit_paths = [trace_field_line(reconstructed, seed, domain) for seed in seeds]

fig_vectors = Figure(; size = (1150, 570), fontsize = 16)
geometry_axis = Axis3(
    fig_vectors[1, 1]; title = "Field geometry: identical seeds",
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
    aspect = :data, azimuth = 1.15pi, elevation = 0.25pi
)
for path in true_paths
    lines!(geometry_axis, Point3f.(path); color = :steelblue, linewidth = 3)
end
for path in fit_paths
    lines!(
        geometry_axis, Point3f.(path); color = :darkorange,
        linewidth = 2, linestyle = :dash
    )
end
arrows3d!(
    geometry_axis, Point3f.(seeds),
    Vec3f.([0.12 * normalize(truth(p)) for p in seeds]);
    color = :steelblue, markerscale = 1, minshaftlength = 0,
    tiplength = 0.05, tipradius = 0.028, shaftradius = 0.003
)
scatter!(
    geometry_axis, Point3f.(eachrow(probe_positions));
    color = (:gray, 0.45), markersize = 5
)
Legend(
    fig_vectors[2, 1],
    [
        LineElement(; color = :steelblue, linewidth = 3),
        LineElement(; color = :darkorange, linewidth = 2, linestyle = :dash),
        MarkerElement(; color = :gray, marker = :circle, markersize = 7)
    ],
    ["Analytical (+B arrows)", "Reconstructed", "Probes"];
    orientation = :horizontal, framevisible = false
)

error_points = [
    SVector(x, y, z) for x in (-0.37, -0.13, 0.11, 0.39),
        y in (-0.31, -0.07, 0.17, 0.37), z in (-0.23, 0.03, 0.29)
]
volume_error = [1.0e6 * norm(reconstructed(p) - truth(p)) for p in error_points]
error_axis = Axis3(
    fig_vectors[1, 2]; title = "Held-out vector error",
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
    aspect = :data, azimuth = 1.15pi, elevation = 0.25pi
)
error_plot = scatter!(
    error_axis, Point3f.(vec(error_points)); color = vec(volume_error),
    colormap = :inferno, colorrange = (0.0, maximum(volume_error)), markersize = 13
)
Colorbar(fig_vectors[1, 3], error_plot; label = "Vector error [uT]")
Label(fig_vectors[2, 2], "Error = norm(B reconstructed - B analytical)")
for ax in (geometry_axis, error_axis)
    limits!(ax, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
end

@assert all(
    all(domain[2k - 1] - 1e-12 <= p[k] <= domain[2k] + 1e-12 for k in 1:3)
        for path in vcat(true_paths, fit_paths) for p in path
)
@assert all(
    minimum(abs(p[k] - domain[j]) for k in 1:3 for j in (2k - 1, 2k)) <= 0.01
        for path in vcat(true_paths, fit_paths) for p in (first(path), last(path))
)
refined_paths = [
    trace_field_line(reconstructed, seed, domain; ds = 0.005) for seed in seeds
]
@assert maximum(
    norm(endpoint(coarse) - endpoint(fine))
        for (coarse, fine) in zip(fit_paths, refined_paths) for endpoint in (first, last)
) < 0.015
fig_vectors
```

The geometry panel intentionally does not amplify differences: the reconstructed
curves should nearly coincide with the analytical ones. The error panel makes
those small differences measurable without implying different field topology.
Field-line endpoints are also compared when the tracing step is halved. Since
FieldTracer does not locate boundary intersections, their separation is checked
against the sum of the coarse and fine step lengths (0.015 m), not a boundary-event
tolerance. Reduce `ds` for finer endpoint resolution; tracing accuracy is separate
from reconstruction accuracy.

## Magnitude and error on a slice

Compare ``|\mathbf B|`` at ``z=0.13`` m, which does not coincide with a probe
plane. All three panels use **millitesla (mT)**. The error panel measures
``|\mathbf B_{\rm reconstructed}-\mathbf B_{\rm analytical}|``, not just the
difference in field magnitudes. For scale, an error of 0.01 mT is 10 microtesla.

```@example reconstruction3d
xs = range(-0.5, 0.5; length = 41)
ys = range(-0.5, 0.5; length = 41)
slice_z = 0.13
true_magnitude = [1.0e3 * norm(truth(SVector(x, y, slice_z))) for x in xs, y in ys]
fit_magnitude = [
    1.0e3 * norm(reconstructed(SVector(x, y, slice_z))) for x in xs, y in ys
]
vector_error = [
    1.0e3 * norm(reconstructed(SVector(x, y, slice_z)) - truth(SVector(x, y, slice_z)))
        for x in xs, y in ys
]
field_range = extrema(vcat(vec(true_magnitude), vec(fit_magnitude)))
panel_size = 280
fig_slice = Figure(; fontsize = 24)
for (column, data, title) in (
        (1, true_magnitude, "Analytical |B|"),
        (2, fit_magnitude, "Reconstructed |B|"),
        (3, vector_error, "Vector error")
    )
    ax = Axis(
        fig_slice[1, 2column - 1]; title, xlabel = "x [m]", ylabel = "y [m]",
        aspect = DataAspect(), width = panel_size, height = panel_size,
        titlesize = 26
    )
    color_range = column == 3 ? (0.0, maximum(vector_error)) : field_range
    image = heatmap!(ax, xs, ys, data; colorrange = color_range, colormap = :viridis)
    Colorbar(
        fig_slice[1, 2column], image;
        label = column == 3 ? "Vector error [mT]" : "|B| [mT]",
        height = panel_size, valign = :center, width = 18
    )
end
resize_to_layout!(fig_slice)
fig_slice
```

## Held-out validation

Evaluate an independent volume grid and a new set of face locations. Neither
is passed to the fit. A finite-difference divergence check verifies the analytic
construction numerically; its residual depends on the differencing step.

```@example reconstruction3d
heldout_positions = points_to_matrix(Iterators.product(
    (-0.37, 0.11, 0.39), (-0.31, 0.17), (-0.23, 0.29)
))
heldout_truth = sample_field(truth, heldout_positions)
heldout_fit = sample_field(reconstructed, heldout_positions)
relative_error = norm(heldout_fit - heldout_truth) / norm(heldout_truth)
heldout_surface, heldout_normals = box_surface((-0.23, 0.21))
heldout_boundary_error = maximum(
    abs(dot(n, reconstructed(p) - truth(p)))
        for (p, n) in zip(heldout_surface, heldout_normals)
)
function numerical_divergence(field, p; h = 1.0e-5)
    return sum(1:3) do k
        step = SVector{3}(ntuple(j -> j == k ? h : 0.0, 3))
        (field(p + step)[k] - field(p - step)[k]) / (2h)
    end
end
maximum_divergence = maximum(
    abs(numerical_divergence(reconstructed, SVector{3}(p)))
        for p in eachrow(heldout_positions)
)
println("Held-out relative vector error: ", round(100relative_error; digits = 2), "%")
println("Held-out normal-boundary error [uT]: ",
    round(1.0e6 * heldout_boundary_error; digits = 2))
println("Maximum numerical |div B| [T/m]: ", maximum_divergence)
println("Weighted observation RMS: ", reconstructed.diagnostics.weighted_rms)

@assert relative_error < 0.03
@assert heldout_boundary_error < 1.0e-4
@assert maximum_divergence < 1.0e-8
```

The held-out error is a synthetic benchmark, not a guarantee for experimental
data. Repeat the fit with different basis resolutions, length scales, and
regularization strengths to assess sensitivity. Noisy normal-boundary data need
not have exactly zero integrated flux; weighted constraints accommodate this
without relaxing the divergence-free representation. Residuals and sampled
support coverage are not posterior uncertainty estimates.
