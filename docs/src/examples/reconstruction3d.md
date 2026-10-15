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

## Synthetic measurements

```@example reconstruction3d
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

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

domain = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
coil = CurrentLoop(
    1.0, 1.0e4, SVector(0.2, -0.1, 2.0), normalize(SVector(0.3, 0.2, 1.0))
)
truth = CurrentLoopAnalytic(coil)
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

Plot field vectors as centered line segments on a coarse 3D grid, with identical
scaling in both panels. Orange points mark probe positions. Segment length is
proportional to field magnitude; segments show orientation, not polarity.

```@example reconstruction3d
vector_points = [
    SVector(x, y, z) for x in (-0.35, 0.0, 0.35),
        y in (-0.35, 0.0, 0.35), z in (-0.35, 0.0, 0.35)
]
vector_scale = 0.22 / maximum(norm(truth(p)) for p in vector_points)
fig_vectors = Figure(; size = (1000, 480), fontsize = 16)
for (column, field, title) in (
        (1, truth, "Analytical coil field"),
        (2, reconstructed, "Reconstructed from sparse probes")
    )
    ax = Axis3(
        fig_vectors[1, column]; title, xlabel = "x [m]", ylabel = "y [m]",
        zlabel = "z [m]", aspect = :data
    )
    segments = Point3f[]
    for p in vector_points
        half_vector = vector_scale * field(p) / 2
        push!(segments, Point3f(p - half_vector), Point3f(p + half_vector))
    end
    linesegments!(ax, segments; color = :steelblue, linewidth = 3)
    scatter!(
        ax, Point3f.(eachrow(probe_positions)); color = :darkorange, markersize = 7
    )
    limits!(ax, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
end
fig_vectors
```

## Magnitude and error on a slice

Compare ``|\mathbf B|`` at ``z=0.13`` m, which does not coincide with a probe
plane. The error panel measures the **vector difference**, not just the
difference in field magnitudes.

```@example reconstruction3d
xs = range(-0.5, 0.5; length = 41)
ys = range(-0.5, 0.5; length = 41)
slice_z = 0.13
true_magnitude = [1.0e3 * norm(truth(SVector(x, y, slice_z))) for x in xs, y in ys]
fit_magnitude = [
    1.0e3 * norm(reconstructed(SVector(x, y, slice_z))) for x in xs, y in ys
]
vector_error = [
    1.0e6 * norm(reconstructed(SVector(x, y, slice_z)) - truth(SVector(x, y, slice_z)))
        for x in xs, y in ys
]
field_range = extrema(vcat(vec(true_magnitude), vec(fit_magnitude)))
fig_slice = Figure(; size = (1150, 400), fontsize = 16)
for (column, data, title) in (
        (1, true_magnitude, "Analytical |B|"),
        (2, fit_magnitude, "Reconstructed |B|"),
        (3, vector_error, "Vector reconstruction error")
    )
    ax = Axis(
        fig_slice[1, 2column - 1]; title, xlabel = "x [m]", ylabel = "y [m]",
        aspect = DataAspect()
    )
    color_range = column == 3 ? (0.0, maximum(vector_error)) : field_range
    image = heatmap!(ax, xs, ys, data; colorrange = color_range, colormap = :viridis)
    Colorbar(fig_slice[1, 2column], image; label = column == 3 ? "Error [uT]" : "|B| [mT]")
end
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
