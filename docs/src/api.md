# API Reference

For method selection and workflows, start with [Forward methods](@ref forward_methods)
or [Inverse methods](@ref inverse_methods). This page groups the public API by role.

## Core interfaces

```@docs
AbstractMagneticField
AbstractCurrentSource
AbstractSolver
solve
```

## Forward construction

### Current sources

```@docs
Wire
CurrentLoop
TorusKnot
TrefoilKnot
InfiniteWire
SurfaceCurrentMesh
```

### Analytical fields and configurations

```@docs
UniformField
NullField
Dipole
CurrentLoopAnalytic
HarrisSheet
AsymmetricHarrisSheet
ForceFreeHarrisSheet
BifurcatedHarrisSheet
GuideFieldHarrisSheet
FadeevIsland
GoldHoyleFluxTube
MagneticMirror
MagneticBottle
TokamakCoils
TokamakProfile
ZPinch
SolovevEquilibrium
TaylorState
ScrewPinch
FLARECoils
ParkerSpiral
AnalyticalMagnetosphere
getB_loop
getB_mirror
getB_bottle
getB_tokamak_coil
getB_tokamak_profile
getB_zpinch
draped_imf_field
```

### Numerical solvers and potentials

```@docs
BiotSavart
FFTSolver
VectorPotential
PoissonSolver
vector_potential
```

### Boundary conditions and surface currents

```@docs
AbstractBoundary
OpenBoundary
PeriodicBoundary
ConductingWall
ConductingSphere
PrescribedBoundary
sphere_mesh
compute_surface_current
image_source
```

## Inverse reconstruction

### Axisymmetric probe fitting

```@docs
AxisymmetricFieldReconstruction
reconstruct_axisymmetric_field
```

### Three-dimensional probe fitting

```@docs
FieldReconstruction3D
reconstruct_field_3d
```

### Spacecraft Grad-Shafranov reconstruction

```@docs
GradShafranovReconstruction
reconstruct_grad_shafranov
```

### Choosing a reconstruction model

See [Choosing an inverse method](@ref inverse_methods) for the comparison of
assumptions, measurement formats, constraints, and validation requirements.

### Three-dimensional reconstruction

See the [3D method guide](@ref inverse_3d) for a minimal fit, physical constraints,
regularization, diagnostics, and limitations.

### Reconstruction performance

See [Reconstruction performance](@ref reconstruction_performance) for the
accuracy/cost sweep and guidance on fitting and evaluation costs.

## Field evaluation and diagnostics

These operations apply to the field types documented by each method; not every
field supports every diagnostic.

```@docs
current_density
poloidal_flux
flux_function
thermal_pressure
```

### Magnetic geometry and plasma drifts

For full orbit trajectory integration and particle pushers, see
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
```

## Utilities

```@docs
discretize_loop
discretize_knot
set_current_wire!
set_current_wire
compute_curl!
compute_curl
sph2cart
dipole_fieldline
```
