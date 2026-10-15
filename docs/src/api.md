# API Reference

## Types

```@docs
AbstractMagneticField
AbstractCurrentSource
AbstractSolver
Wire
CurrentLoop
HarrisSheet
AsymmetricHarrisSheet
ForceFreeHarrisSheet
BifurcatedHarrisSheet
FadeevIsland
Dipole
AnalyticalMagnetosphere
TorusKnot
TrefoilKnot
InfiniteWire
CurrentLoopAnalytic
UniformField
NullField
MagneticMirror
MagneticBottle
TokamakCoils
TokamakProfile
ZPinch
SolovevEquilibrium
TaylorState
ScrewPinch
ParkerSpiral
GuideFieldHarrisSheet
GoldHoyleFluxTube
FLARECoils
AxisymmetricFieldReconstruction
FieldReconstruction3D
```

## Boundary Conditions

```@docs
AbstractBoundary
OpenBoundary
PeriodicBoundary
ConductingWall
ConductingSphere
PrescribedBoundary
```

## Solvers

```@docs
BiotSavart
FFTSolver
VectorPotential
PoissonSolver
reconstruct_axisymmetric_field
reconstruct_field_3d
solve
vector_potential
```

### Three-dimensional reconstruction

`reconstruct_field_3d` is separate from the axisymmetric reconstruction API.
It fits Cartesian vector measurements using

```math
\mathbf B = \mathbf B_{\rm ref} + \mathbf b_0
  + \nabla\times\left[\ell\sum_j\phi_j(\mathbf x)\mathbf a_j\right].
```

The fitted correction is divergence-free throughout the evaluation box, not
only at probe locations. A supplied reference field must also be divergence-free.
The uniform mode captures background fields that localized bases represent poorly.

For example, the following fits synthetic measurements of a tilted external coil:

```julia
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

**Limitations:** finite probes and normal boundary samples do not uniquely
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

## Magnetic Geometry and Plasma Drifts

Differential geometry quantities and guiding-center drift velocity formulas. For full
orbit trajectory integration and particle pushers, see
[TestParticle.jl](https://github.com/henry2004y/TestParticle.jl).

```@docs
b_unit
field_magnitude
grad_b
magnetic_curvature
drift_velocities
mirror_ratio
loss_cone_angle
magnetic_helicity
poloidal_flux
```

## Surface Current (BEM)

```@docs
SurfaceCurrentMesh
sphere_mesh
compute_surface_current
```

## Utilities

```@docs
discretize_loop
discretize_knot
getB_loop
set_current_wire!
set_current_wire
compute_curl!
compute_curl
current_density
sph2cart
dipole_fieldline
image_source
draped_imf_field
```

## Configurations

```@docs
getB_mirror
getB_bottle
getB_tokamak_coil
getB_tokamak_profile
getB_zpinch
```
