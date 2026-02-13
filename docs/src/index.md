# Magnetostatics.jl

A Julia package for computing magnetostatic fields from current sources.

## Overview

Magnetostatics.jl provides solvers and analytical models for computing magnetic fields produced by steady-state current distributions. The package includes:

- **Biot-Savart solver** — numerical integration of the Biot-Savart law for arbitrary wire geometries.
- **FFT solver** — spectral method for computing **B** from a volumetric current density **J** on a uniform grid with periodic boundaries.
- **Vector Potential solver** — computes the magnetic vector potential **A** for wires, current loops, and dipoles.
- **Analytical fields** — closed-form models for Harris current sheets, magnetic dipoles, and circular current loops (via elliptic integrals).

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
Pages = ["api.md"]
```
