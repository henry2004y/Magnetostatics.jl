# [Reconstruction Across Sources and Boundaries](@id reconstruction_cases_example)

The [Biot-Savart demo](@ref biot_savart_example) solves the **forward** problem:
given sources and boundaries, compute the field. Here those configurations
generate synthetic observations for the **inverse** problem. The reconstruction
receives positions, field samples and explicitly supplied boundary/current data,
not the source geometry. All fits below use the general Cartesian 3D model;
see [Choosing a reconstruction model](../api.md#Choosing-a-reconstruction-model)
for the distinction from axisymmetric flux reconstruction.

| Case | Synthetic truth | Information given to the fit |
|:--|:--|:--|
| Open loop | Biot-Savart sum over 150 tilted-loop segments | Vector probes only; no artificial condition on the box faces |
| Conducting wall | Finite wire plus its opposite-current image | Vector probes and exact zero normal field at 25 wall points |
| Conducting sphere | Uniform background plus the analytic image dipole | Exterior probes, known uniform background and 25 uncertain zero-normal samples on a spherical cap |
| Volume current | Interior of a uniform-current Z-pinch | Vector probes and 27 vector-current observations |

The sphere uses the analytic limit from the forward demo, rather than a
discretized BEM field, to avoid conflating BEM error with reconstruction error.
Its box includes a small part of the sphere: probes, validation points and plots
are restricted to the **exterior**. The fitted smooth continuation inside the
sphere is not the physical interior solution. The cap is not a closed surface,
and this example does not reconstruct a global multiply connected exterior.

## Problems and measurements

```@example reconstruction_cases
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie, Markdown

point_matrix(points) = [p[k] for p in vec(collect(points)), k in 1:3]
sample_field(field, positions) =
    point_matrix(field(SVector{3}(p)) for p in eachrow(positions))
function lattice(box, coordinates)
    return point_matrix(Iterators.product(
        ntuple(k -> [
            box[2k - 1] + t * (box[2k] - box[2k - 1]) for t in coordinates
        ], 3)...
    ))
end

box = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
wall_box = (-0.5, 0.5, -0.5, 0.5, 0.0, 0.6)
sphere_box = (-0.25, 0.25, -0.25, 0.25, 0.9, 1.4)
coil = discretize_loop(CurrentLoop(
    1.0, 1.0e4, SVector(0.2, -0.1, 2.0), normalize(SVector(0.3, 0.2, 1.0))
), 150)
wire = Wire([SVector(0.0, -10.0, 1.5), SVector(0.0, 10.0, 1.5)], 1.0e4)
background = UniformField(SVector(0.0, 0.0, 0.01))
image = Dipole(SVector(0.0, 0.0, -2π * 0.01 / Magnetostatics.μ₀))
pinch = ZPinch(1.0e4, 2.0)

wall_points = point_matrix(
    SVector(x, y, 0.0) for x in range(-0.4, 0.4; length = 5),
        y in range(-0.4, 0.4; length = 5)
)
sphere_points = point_matrix(
    SVector(x, y, sqrt(1 - x^2 - y^2))
        for x in range(-0.2, 0.2; length = 5), y in range(-0.2, 0.2; length = 5)
)
current_points = lattice(box, (0.1, 0.5, 0.9))
current_values = point_matrix(current_density(pinch, p) for p in eachrow(current_points))
cases = (
    (;
        name = "Open loop", truth = p -> solve(BiotSavart(), coil, p),
        domain = box, scale = 1.2, constraints = (;), exterior = false,
    ),
    (;
        name = "Conducting wall",
        truth = p -> solve(BiotSavart(ConductingWall(3, 0.0)), wire, p),
        domain = wall_box, scale = 0.8, exterior = false,
        constraints = (;
            boundary_positions = wall_points,
            boundary_normals = repeat([0.0 0.0 1.0], size(wall_points, 1)),
            boundary_values = zeros(size(wall_points, 1)),
        ),
    ),
    (;
        name = "Sphere exterior", truth = p -> background(p) + image(p),
        domain = sphere_box, scale = 0.8, exterior = true,
        constraints = (;
            reference_field = background,
            boundary_positions = sphere_points, boundary_normals = sphere_points,
            boundary_values = zeros(size(sphere_points, 1)),
            boundary_uncertainties = fill(1.0e-6, size(sphere_points, 1)),
        ),
    ),
    (;
        name = "Volume current", truth = pinch,
        domain = box, scale = 1.2, exterior = false,
        constraints = (;
            current_positions = current_points, current_density = current_values,
            current_uncertainties = fill(1.0, size(current_values)),
        ),
    ),
)
nothing # hide
```

All source currents are 10 kA; distances are in metres and the sphere radius
is 1 m in a 10 mT background. The Z-pinch radius is 2 m, so the entire fitted
box is inside the distributed-current region, away from the discontinuous edge.
There is no default current-free assumption. In particular, field-energy
regularization is not a substitute for current measurements.

## Fit and validate off the probe grid

```@example reconstruction_cases
results = map(cases) do problem
    positions = lattice(problem.domain, (0.05, 0.35, 0.65, 0.95))
    heldout = lattice(problem.domain, (0.17, 0.53, 0.81))
    if problem.exterior
        positions = positions[[norm(p) > 1.02 for p in eachrow(positions)], :]
        heldout = heldout[[norm(p) > 1.02 for p in eachrow(heldout)], :]
    end
    observed = sample_field(problem.truth, positions)
    fitted = reconstruct_field_3d(
        positions, observed;
        domain = problem.domain, length_scale = problem.scale, kernel = :gaussian,
        basis_centers = lattice(problem.domain, range(0, 1; length = 4)),
        uncertainties = fill(1.0e-5, size(observed)),
        regularization = 1.0e-8, regularization_grid = (5, 5, 5),
        abstol = 1.0e-13, reltol = 1.0e-13, problem.constraints...
    )
    expected = sample_field(problem.truth, heldout)
    error = norm(sample_field(fitted, heldout) - expected) / norm(expected)
    @assert error < (problem.exterior ? 0.05 : 0.03)
    (; problem, fitted, error, probes = size(positions, 1), heldout)
end

rows = [
    "| $(r.problem.name) | $(r.probes) | $(round(100r.error; digits = 2))% |"
        for r in results
]
Markdown.parse(join(vcat(
    ["| Case | Vector probes | Held-out vector error |", "|:--|--:|--:|"], rows
), "\n"))
```

These are noiseless synthetic values, weighted with finite uncertainty scales
to specify their relative influence; the scales are not estimated error bars.
The [tilted-coil demo](@ref reconstruction3d_example) separately includes noise
and missing components. The error above is the aggregate vector norm
``\|B_{\rm fit}-B_{\rm truth}\|_F/\|B_{\rm truth}\|_F``, not magnitude-only error
or a pointwise relative error near a magnetic null.

Check **new** surface points and volume currents, not just training residuals:

```@example reconstruction_cases
wall_fit = results[2].fitted
wall_error = maximum(
    abs(wall_fit(SVector(x, y, 0.0))[3])
        for x in (-0.27, 0.21), y in (-0.19, 0.31)
)
sphere_fit = results[3].fitted
cap_error = maximum(
    abs(dot(p, sphere_fit(p))) for p in (
        SVector(x, y, sqrt(1 - x^2 - y^2))
            for x in (-0.13, 0.17), y in (-0.11, 0.15)
    )
)
pinch_fit = results[4]
J_expected = point_matrix(
    current_density(pinch, p) for p in eachrow(pinch_fit.heldout)
)
J_fitted = point_matrix(
    current_density(pinch_fit.fitted, p) for p in eachrow(pinch_fit.heldout)
)
current_error = norm(J_fitted - J_expected) / norm(J_expected)
@assert wall_error < 1.0e-5
@assert cap_error < 1.0e-4
@assert current_error < 0.02
println("Held-out wall normal error [T]: ", wall_error)
println("Held-out spherical-cap normal error [T]: ", cap_error)
println("Held-out relative current error: ", current_error)
```

Exact wall samples do not enforce zero normal field **between** samples.
Similarly, the soft cap samples do not impose a global spherical boundary
condition. The held-out checks measure these limitations explicitly.

## Compare fields and errors

Top panels show the true field magnitude in the ``y=0`` slice. Bottom panels
show vector error on the same slice, with separate color scales in microtesla.
White regions in the sphere panels are masked physical interior, not zero error.

```@example reconstruction_cases
fig_cases = Figure(; size = (1500, 780), fontsize = 16)
for (column, result) in enumerate(results)
    problem, fitted = result.problem, result.fitted
    xs = range(problem.domain[1], problem.domain[2]; length = 41)
    zs = range(problem.domain[5], problem.domain[6]; length = 41)
    magnitude = [
        problem.exterior && hypot(x, z) < 1 ? NaN :
            1.0e3 * norm(problem.truth(SVector(x, 0.0, z))) for x in xs, z in zs
    ]
    error = [
        problem.exterior && hypot(x, z) < 1 ? NaN :
            1.0e6 * norm(fitted(SVector(x, 0.0, z)) -
                problem.truth(SVector(x, 0.0, z))) for x in xs, z in zs
    ]
    for (row, values, title, units) in (
            (1, magnitude, problem.name, "|B| [mT]"),
            (2, error, "Held-out slice error", "|Bfit - B| [uT]"),
        )
        ax = Axis(
            fig_cases[row, 2column - 1]; title, xlabel = "x [m]", ylabel = "z [m]",
            aspect = DataAspect()
        )
        image_plot = heatmap!(ax, xs, zs, values; colormap = :viridis)
        Colorbar(fig_cases[row, 2column], image_plot; label = units)
        if problem.exterior
            lines!(ax, xs, sqrt.(1 .- xs.^2); color = :white, linewidth = 2)
        elseif column == 2
            hlines!(ax, [0.0]; color = :white, linewidth = 2)
        end
    end
end
fig_cases
```

## Accuracy versus computational cost

The standalone sweep `benchmark/reconstruction_scaling.jl` varies kernel, center
count, probe count, length scale and regularization, reporting accuracy alongside
fit/evaluation time and allocations. See [performance guidance](../api.md#Reconstruction-performance).

On the external-coil sweep with deterministic 10 microtesla noise, increasing
Gaussian centers from 64 to 125 at 125 probes and scale 1.2 m changed held-out
field error only from about **1.74% to 1.73%**, while roughly doubling fit cost.
For that same 64-center fit the dimensionless spurious-curl measure was about
**0.18**: good field agreement does not imply equally accurate currents.
Compact Wendland fits with the same data could be much worse despite full
sampled support coverage. A few dozen probes cannot resolve arbitrary local
current structure. These results are problem-specific, not a universal ranking
of kernels.

The optimized evaluator sums basis contributions directly, and fitting assembles
sparse matrices from triplets instead of repeatedly concatenating sparse rows.
For a warmed, noiseless 27-probe coil fit, support scale 1.2 m and a 5-by-5-by-5
regularization grid, representative before/after measurements on Windows,
Julia 1.13.1, one thread, were:

| Kernel / centers | Fit before / after | Fit allocation before / after | Field query before / after |
|:--|--:|--:|--:|
| Gaussian / 27 | 174 / 5.9 ms | 257 / 10.8 MB | 9.1 / 0.63 us |
| Gaussian / 125 | 665 / 42.6 ms | 1159 / 45.2 MB | 32.1 / 2.94 us |
| Wendland / 27 | 58.3 / 5.3 ms | 209 / 9.2 MB | 7.1 / 1.25 us |
| Wendland / 125 | 552 / 42.0 ms | 1021 / 41.3 MB | 29.9 / 6.26 us |

Query allocations fell from roughly 22--84 kB to **32 bytes** per returned
field/current vector in this benchmark. Timings exclude compilation and vary
with hardware and load; rerun the benchmark rather than treating these as
performance guarantees. Direct summation is regression-tested against the sparse
operators for both kernels, including centers, support edges and box boundaries.
