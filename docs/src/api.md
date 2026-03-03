# API Reference

## Types

```@docs
AbstractMagneticField
AbstractCurrentSource
AbstractSolver
Wire
CurrentLoop
HarrisSheet
Dipole
CurrentLoopAnalytic
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
getB_loop
set_current_wire!
set_current_wire
compute_curl!
compute_curl
sph2cart
dipole_fieldline
image_source
```

## Configurations

```@docs
getB_mirror
getB_bottle
getB_tokamak_coil
getB_tokamak_profile
getB_zpinch
```

