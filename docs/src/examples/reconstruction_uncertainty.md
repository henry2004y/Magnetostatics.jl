# [Reconstruction Uncertainty and Hidden Topology](@id reconstruction_uncertainty_example)

Can a sparse-probe reconstruction reliably locate and classify a magnetic null?
This example separates **measurement-noise propagation within a fixed model**
from **ambiguity outside that model**. It does not introduce a Bayesian solver
or interpret ensemble frequencies as posterior probabilities.

The same executable benchmark is used by the reconstruction uncertainty tests.

## Controlled sparse-probe experiment

The domain is a two-metre cube, with **18 vector probes** and optional **24
normal-field boundary measurements**. The synthetic field has an isolated null
at the origin with two contracting and one expanding eigen-directions. It is
manufactured from the same eight-centre Gaussian vector-potential basis used
for fitting. This favorable, deliberately well-specified case isolates noise
propagation; it is not evidence that eight centres suffice for an unknown plasma.

We generate one nominal dataset with independent Gaussian component errors of
10 microtesla and normal-boundary errors of 5 microtesla. Boundary samples are
uncertain observations of the true nonzero normal field, **not conducting walls**.
The fit receives noisy observations and basis parameters, not the true coefficients
or null position.

```@example reconstruction_uncertainty
include("reconstruction_uncertainty.jl")
using .ReconstructionUncertaintyDemo
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

demo = ReconstructionUncertaintyDemo
experiment = demo.run_ensemble(; pairs = 64)
data = experiment.data
println("Vector probes: ", size(data.positions, 1))
println("Normal boundary samples: ", length(data.boundary))
println("Measurement sigma [uT]: ", 1.0e6 * experiment.sigma)
```

For each observation type, we draw 64 independent Gaussian perturbations and fit
both signs of each draw: **128 antithetic ensemble members**, not 128 independent
experiments. We repeat this with and without boundary measurements, perturbing
the uncertain boundary data as well as the probes. Paired draws make the ensemble
field mean equal to the nominal fit, providing a numerical linearity check.

The basis, length scale, and coefficient-ridge penalty are fixed across ensemble
members. A direct QR solve avoids confusing iterative convergence variation with
noise sensitivity. For the separate noise-scaling check, the code keeps
`regularization * sigma^2 = 1e-4` fixed so the effective unweighted estimator is
unchanged. This is a controlled sensitivity experiment, not an inferred prior
strength or an automatic regularization prescription.

## Check the propagated field uncertainty independently

For this linear, fixed-basis estimator, an independent dense design matrix gives
the exact covariance of its response to the assumed Gaussian errors:

```math
C_B = G K K^\mathsf{T} G^\mathsf{T},\qquad
K = (H^\mathsf{T}H+\lambda I)^{-1}H^\mathsf{T}.
```

Here `H` is the uncertainty-whitened measurement operator and `G` evaluates the
field at a new point. The implementation solves the augmented least-squares
problem with QR rather than explicitly forming an inverse. The covariance
describes **estimator variability due to measurement errors**, not regularization
bias, missing basis functions, uncertain geometry, or a posterior distribution.

```@example reconstruction_uncertainty
point = SVector(0.13, 0.17, -0.23)
covariances = [
    demo.prediction_covariance(data, point; sigma = experiment.sigma, constrained)
        for constrained in (false, true)
]
for (name, result, covariance) in zip(
        ("Probes only", "Probes + boundary"),
        (experiment.free, experiment.bounded), covariances
    )
    draws = hcat((field(point) for field in result.fields)...)
    centered = draws .- sum(draws; dims = 2) / size(draws, 2)
    empirical = centered * centered' / size(draws, 2)
    relative_difference = norm(empirical - covariance) / norm(covariance)
    println(name, ": covariance sampling discrepancy = ", relative_difference)
    @assert relative_difference < 0.3
    @assert isapprox(
        result.summary.average, demo.samples(result.nominal, data.heldout); atol = 1.0e-12
    )
end
@assert tr(covariances[2]) < tr(covariances[1])

fig_uncertainty = Figure(; size = (1100, 470), fontsize = 18)
ax_sigma = Axis(
    fig_uncertainty[1, 1]; title = "Conditional field standard deviation",
    xlabel = "Cartesian component", ylabel = "Standard deviation [uT]",
    xticks = (1:3, ["Bx", "By", "Bz"])
)
for (offset, covariance, label, color) in (
        (-0.16, covariances[1], "Probes only", :darkorange),
        (0.16, covariances[2], "Probes + boundary", :steelblue)
    )
    barplot!(
        ax_sigma, collect(1:3) .+ offset, 1.0e6 .* sqrt.(diag(covariance));
        width = 0.3, label, color
    )
end
Legend(fig_uncertainty[2, 1], ax_sigma; orientation = :horizontal, framevisible = false)

ax_null = Axis(
    fig_uncertainty[1, 2]; title = "Recovered nulls: x-z projection",
    xlabel = "x [m]", ylabel = "z [m]"
)
for (result, label, color) in (
        (experiment.free, "Probes only", :darkorange),
        (experiment.bounded, "Probes + boundary", :steelblue)
    )
    @assert all(n.status == :converged && n.positive == 1 for n in result.nulls)
    scatter!(
        ax_null, [n.point[1] for n in result.nulls], [n.point[3] for n in result.nulls];
        color = (color, 0.6), markersize = 8, label
    )
end
scatter!(ax_null, [0.0], [0.0]; color = :black, marker = :star5, markersize = 18,
    label = "True null")
Legend(fig_uncertainty[2, 2], ax_null; orientation = :horizontal, framevisible = false)
fig_uncertainty
```

In this experiment, boundary observations reduce aggregate field and null-location
spread. This is not a guarantee for every field component, geometry, or finite
ensemble. Nulls are found by a local Newton search with explicit failure statuses;
the search is not an exhaustive null census. The tests also halve the Jacobian
differencing step, check null residuals and eigenvalue signs, and verify that
doubling the noise doubles field spread when the estimator is held fixed.

The displayed cloud is a sensitivity visualization, **not a 95% confidence
region**. Even in a correct model, a finite cloud has no guaranteed coverage.
Ridge bias and the nonlinear map from field coefficients to nulls require separate
assessment. The held-out volume grid is never supplied to the reconstruction.

## Counterexample: identical data, different null structure

Now add a compactly supported curl field entirely between the probes and away
from the boundary. It is exactly zero at every probe and on the **entire cube
boundary**. Both fields remain divergence-free, but the central null changes
from one expanding eigen-direction to two, and at least two additional isolated
nulls appear along the x axis.
No force-balance condition, upper bound on current, or current-free prior is
imposed here: the ambiguity concerns sparse field and normal-boundary data alone.

```@example reconstruction_uncertainty
ambiguity = demo.ambiguity_case()
alternative = ambiguity.alternative
@assert demo.samples(alternative, data.positions) == demo.samples(data.truth, data.positions)
@assert all(
    dot(n, alternative(p)) == value
        for (p, n, value) in zip(
            eachrow(data.surface), eachrow(data.normals), data.boundary
        )
)
central = demo.locate_null(alternative)
@assert central.status == :converged && central.positive == 2
for root in (ambiguity.negative_root, ambiguity.positive_root)
    @assert root.status == :converged && root.residual < 1.0e-11
    @assert norm(data.truth(root.point)) > 1.0e-5
end
println("Additional null positions [m]: ",
    (ambiguity.negative_root.point, ambiguity.positive_root.point))

xs = range(-0.35, 0.35; length = 301)
fig_hidden = Figure(; size = (900, 430), fontsize = 18)
ax_hidden = Axis(
    fig_hidden[1, 1]; title = "An unobserved current changes the null structure",
    xlabel = "x [m], at y = z = 0", ylabel = "Bx [uT]"
)
lines!(ax_hidden, xs, [1.0e6 * data.truth(SVector(x, 0.0, 0.0))[1] for x in xs];
    label = "Baseline field", color = :steelblue, linewidth = 3)
lines!(ax_hidden, xs, [1.0e6 * alternative(SVector(x, 0.0, 0.0))[1] for x in xs];
    label = "Same observations + hidden current", color = :darkorange, linewidth = 3)
hlines!(ax_hidden, [0.0]; color = :gray, linestyle = :dash)
scatter!(
    ax_hidden,
    [ambiguity.negative_root.point[1], 0.0, ambiguity.positive_root.point[1]], zeros(3);
    color = :black, marker = :diamond, markersize = 12, label = "Verified nulls"
)
axislegend(ax_hidden; position = :lt)
fig_hidden
```

On this symmetry axis, `By = Bz = 0`, so these crossings of `Bx` are full nulls;
the tests verify all three components and local nondegeneracy. This plot does
not exclude additional nulls elsewhere.

The compact perturbation lies outside the fitted Gaussian basis. Identical data
therefore produce identical conditional ensembles whether the baseline or
alternative field is the actual truth. Near the hidden current, the two fields
differ by more than ten conditional RMS standard deviations in this benchmark.
The narrow ensemble cannot reveal what its model excludes.

## What this establishes

- **Achieved:** reproducible propagation of explicitly assumed measurement errors,
  an independent covariance check for the linear field estimator, numerical checks
  on null recovery, and a controlled improvement from additional boundary data.
- **Not established:** posterior probabilities, general coverage of topology,
  uniqueness, or robustness to all physically admissible currents.
- **Practical consequence:** combine noise ensembles with basis/prior sensitivity
  and searches for alternative fields. Measurements near the hidden current,
  independently justified current constraints, or a valid current-free assumption
  can discriminate the two fields. More samples of the same normal boundary
  cannot distinguish this interior perturbation.

See [Reconstructing topology from sparse measurements](@ref inverse_methods)
for the broader identifiability question and
[3D reconstruction](@ref inverse_3d) for the underlying fitting API.
