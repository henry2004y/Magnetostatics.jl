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

