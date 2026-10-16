# Magnetostatics.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://henry2004y.github.io/Magnetostatics.jl/dev)
[![Coverage](https://codecov.io/gh/henry2004y/Magnetostatics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/henry2004y/Magnetostatics.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/refs/heads/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Julia package for forward construction of magnetostatic fields from sources
or potentials, and inverse reconstruction from measurements.

## Overview

Magnetostatics.jl provides analytical models and numerical solvers for fields
produced by steady-state currents or potentials, alongside methods for
reconstructing fields from sparse measurements. The package includes:

- **Biot-Savart solver** — numerical integration of the Biot-Savart law for arbitrary wire geometries.
- **FFT solver** — spectral method for computing **B** from a volumetric current density **J** on a uniform grid with periodic boundaries.
- **Vector Potential solver** — computes the magnetic vector potential **A** for wires, current loops, and dipoles.
- **3D reconstruction** — divergence-free field fitting from sparse Cartesian or
  directional probes, with normal-boundary constraints and current measurements.
- **Axisymmetric reconstruction** — poloidal-flux fitting from cylindrical probes,
  with vacuum, external-filament and plasma-current bases for axisymmetric systems.
- **Spacecraft Grad-Shafranov reconstruction** — planar magnetohydrostatic
  reconstruction from a magnetic/pressure time series in a supplied crossing frame.
- **Magnetosphere models** — planetary dipoles, magnetopause shielding, magnetotail, and interplanetary magnetic field (IMF) components.
- **Tokamaks** — analytical magnetic field models for toroidal/poloidal coils, central solenoids, and q-profile based configurations.
- **Analytical fields** — closed-form models for Harris current sheets (asymmetric, force-free), magnetic dipoles, mirrors, pinches, and circular current loops.

## Installation

To install Magnetostatics.jl, run the following command in the Julia REPL:

```julia
using Pkg
Pkg.add("Magnetostatics")
```

## Quick Start

```julia
using Magnetostatics

# Create a circular current loop and discretize it into wire segments
loop = CurrentLoop(1.0, 1.0, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
wire = discretize_loop(loop, 100)

# Compute B at a point using the Biot-Savart solver
B = solve(BiotSavart(), wire, [0.0, 0.0, 0.5])

# Compare with the analytical solution
B_exact = getB_loop([0.0, 0.0, 0.5], loop)
```

## Documentation

For more detailed information on the API and usage, please refer to the [documentation](https://henry2004y.github.io/Magnetostatics.jl/dev).

- [Forward methods](https://henry2004y.github.io/Magnetostatics.jl/dev/forward):
  choose analytical models and numerical solvers from known sources and boundaries.
- [Inverse methods](https://henry2004y.github.io/Magnetostatics.jl/dev/inverse/):
  compare reconstruction assumptions, measurement formats, constraints, and validation.
- [API reference](https://henry2004y.github.io/Magnetostatics.jl/dev/api):
  public interfaces grouped by forward construction, inverse reconstruction,
  and shared field diagnostics.

## Related Packages

- [FieldTracer.jl](https://github.com/henry2004y/FieldTracer.jl) — Field-line tracing
  for analytical and reconstructed fields, and numerical fields on grids. The
  examples use v0.2.4 or later; UniformStreamlines.jl remains useful for automatically
  placing evenly spaced streamlines.
- [TestParticle.jl](https://github.com/henry2004y/TestParticle.jl) — Charged test particle
  tracing in electromagnetic fields, including guiding-center drifts and relativistic pushers.
