# Magnetostatics.jl

A Julia package for computing magnetostatic fields from current sources.

## Overview

Magnetostatics.jl provides solvers and analytical models for computing magnetic fields produced by steady-state current distributions.

## Methods

### Analytical Models
Exact solutions for specific geometries:
- **`HarrisSheet`**: Harris current sheet model ($B_x(z) = B_0 \tanh(z/L)$).
- **`Dipole`**: Magnetic dipole field.
- **`CurrentLoop`**: Analytics for circular current loops using elliptic integrals.
- **`getB_mirror`**: Two-coil magnetic mirror configuration.
- **`getB_bottle`**: Magnetic bottle configuration (mirror + central coil).
- **`getB_tokamak_coil`**: Tokamak field from 16 toroidal coils + plasma current.
- **`getB_tokamak_profile`**: Tokamak field from safety factor $q$-profile.
- **`getB_zpinch`**: Z-pinch wire field.

### Numerical Solvers
General-purpose solvers for arbitrary geometries:
- **`BiotSavart`**: Numerical integration of the Biot-Savart law for arbitrary wire geometries.
- **`FFTSolver`**: Spectral method for computing **B** from a volumetric current density **J** on a uniform grid with periodic boundaries.
- **`VectorPotential`**: Computes the magnetic vector potential **A** for wires, current loops, and dipoles.

## Quick Start

```julia
using Magnetostatics
using StaticArrays

# Create a circular current loop and discretize it into wire segments
loop = CurrentLoop(1.0, 1.0, [0, 0, 0], [0, 0, 1])
wire = discretize_loop(1.0, 100, 1.0)

# Compute B at a point using the Biot-Savart solver
B = solve(BiotSavart(), wire, SVector(0.0, 0.0, 0.5))

# Compare with the analytical solution
B_exact = getB_loop(SVector(0.0, 0.0, 0.5), loop)
```

## Contents

```@contents
Pages = ["examples.md", "api.md"]
```
