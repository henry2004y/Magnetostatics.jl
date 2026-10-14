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
solve
vector_potential
```

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
