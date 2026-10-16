# [Three-Dimensional Reconstruction](@id inverse_3d)

[`reconstruct_field_3d`](@ref) fits a continuous Cartesian field from sparse
vector or directional measurements without assuming spatial symmetry.
See [Choosing an inverse method](@ref inverse_methods) before selecting this model.

## Representation

It fits Cartesian vector measurements using

```math
\mathbf B = \mathbf B_{\rm ref} + \mathbf b_0
  + \nabla\times\left[\ell\sum_j\phi_j(\mathbf x)\mathbf a_j\right].
```

The fitted correction is divergence-free throughout the evaluation box, not
only at probe locations. A supplied reference field must also be divergence-free.
The uniform mode captures background fields that localized bases represent poorly.

## Basic workflow

The following fits synthetic measurements of a tilted external coil:

```@example reconstruction_quickstart
using Magnetostatics, StaticArrays, LinearAlgebra

points_to_matrix(points) =
    reduce(vcat, [reshape(collect(p), 1, 3) for p in points])

coil = CurrentLoop(
    1.0, 1.0e4, SVector(0.2, -0.1, 2.0), normalize(SVector(0.3, 0.2, 1.0))
)
reference = CurrentLoopAnalytic(coil)
positions = points_to_matrix(Iterators.product(
    (-0.45, -0.15, 0.15, 0.45), (-0.4, 0.0, 0.4), (-0.4, 0.0, 0.4)
))
observed = points_to_matrix([reference(SVector{3}(p)) for p in eachrow(positions)])
centers = points_to_matrix(Iterators.product(
    ntuple(_ -> range(-0.5, 0.5; length = 4), 3)...
))

field = reconstruct_field_3d(
    positions, observed;
    domain = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
    length_scale = 1.2, basis_centers = centers, kernel = :gaussian,
    uncertainties = fill(1.0e-5, size(observed)),
    regularization = 1.0e-10, abstol = 1.0e-12, reltol = 1.0e-12
)
B = field(SVector(0.11, 0.17, -0.23))
J = current_density(field, SVector(0.11, 0.17, -0.23))
diagnostics = field.diagnostics
```

Replace synthetic observations with calibrated data in tesla and uncertainties
in tesla. Mark missing components with `NaN`. Single-axis probes can instead use
`directional_positions`, unit `directional_directions`, `directional_values`, and
`directional_uncertainties`; use empty `0 × 3` Cartesian matrices for
directional-only reconstruction.

## Constraints and regularization

Physical constraints are supplied independently:

- **Normal boundary field:** `boundary_positions`, unit `boundary_normals`, and
  `boundary_values` in tesla. Constraints are exact at those points unless
  `boundary_uncertainties` is supplied. No boundary condition is assumed by default.
- **Current measurements:** `current_positions`, `current_density`, and
  `current_uncertainties`, with three Cartesian components in A/m². These use
  analytic derivatives of the fitted potential, not finite-difference currents.
- **Known background:** `reference_field`, whose contribution is subtracted from
  observations before fitting a correction. Current observations/evaluation
  require a current-density method for this reference (null/uniform fields are
  supported directly).

By default, regularization penalizes the volume-averaged correction-field energy
and curl energy, with the curl term scaled by `length_scale²`. It does **not**
declare the volume current-free. `regularization_grid` controls the tensor
trapezoidal quadrature; refine it to check prior convergence.
`regularization_mode = :coefficients` selects a coefficient ridge instead.
Choose regularization jointly with measurement uncertainties and length scale.
The default unit uncertainties are placeholders, not calibrated noise estimates.

## Validation and limitations

Finite probes and normal boundary samples do not uniquely
determine the interior field. Check held-out measurements and sensitivity to
basis resolution, length scale, and regularization. Exact surface samples do not
enforce the condition between samples, and closed-surface data must satisfy zero
net magnetic flux physically. A conducting wall requires an explicit flux
assumption; no wall or open-exterior model is inferred from the evaluation box.
Periodic and multiply connected domains require additional global modes and are
not modeled by this initial box-based API.

`field.diagnostics` contains RMS residuals per observation type, weighted RMS,
maximum exact-boundary error, solver status and iterations, and sampled support
coverage. An absent observation type has an RMS of `nothing`. Support gaps generate
a warning: only uniform/reference contributions remain there. Sampled coverage
is not a guarantee of coverage everywhere, and residuals are not posterior
uncertainty estimates. Exact constraints currently use a dense null-space
reduction; use uncertainty-weighted boundary data for the sparse iterative path.
Small length scales increase default center and quadrature counts cubically.

## [Reconstruction performance](@id reconstruction_performance)

Fitting and evaluation have different costs. Fitting assembles sparse observation
and regularization operators; evaluating a fitted 3D field sums analytic basis
contributions directly without constructing those operators. Compact Wendland
evaluation visits neighboring spatial bins; Gaussian evaluation visits all centers.
For a compact basis at fixed center density and support radius, adding remote
centers should not increase the local summation cost. Larger support radii can
make even a compact design matrix effectively dense.

Run the accuracy/cost sweep from the repository root:

```powershell
julia --project=benchmark benchmark\reconstruction_scaling.jl
```

This prints CSV rows with basis/probe counts, kernel, support radius, regularization,
held-out vector error, a dimensionless spurious-current measure, fit time and
allocation, and field/current evaluation time and allocation. Compilation is
excluded by warming each configuration. Noise and held-out points are deterministic.
The existing `SUITE["reconstruction3d"]` benchmark group additionally separates
operator construction, matrix assembly, fitting and evaluation.

Increasing center count is not guaranteed to improve accuracy: sparse probes,
potential gauge modes and ill-conditioning can dominate. Compare held-out error
**and** derivative/current error before selecting the fastest acceptable model.
Refine probe coverage and physical priors, not just the basis. Exact boundary
constraints introduce a dense null-space solve and can dominate larger fits.
The [physical cases](@ref reconstruction_cases_example) also validate normal
fields at surface locations not supplied to the reconstruction.

## Examples and API

- [Tilted-coil reconstruction](@ref reconstruction3d_example): noisy probes,
  uncertain boundary samples, field-line tracing, and held-out validation.
- [Reconstruction across sources and boundaries](@ref reconstruction_cases_example):
  open loops, conducting walls, a conducting sphere, and distributed volume current.
- [`reconstruct_field_3d`](@ref) and [`FieldReconstruction3D`](@ref):
  keyword arguments and the returned field type.
