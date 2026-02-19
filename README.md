# Magnetostatics.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://henry2004y.github.io/Magnetostatics.jl/dev)
[![Coverage](https://codecov.io/gh/henry2004y/Magnetostatics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/henry2004y/Magnetostatics.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Julia package for computing magnetostatic fields from current sources.

## Overview

Magnetostatics.jl provides solvers and analytical models for computing magnetic fields produced by steady-state current distributions. The package includes:

- **Biot-Savart solver** — numerical integration of the Biot-Savart law for arbitrary wire geometries.
- **FFT solver** — spectral method for computing **B** from a volumetric current density **J** on a uniform grid with periodic boundaries.
- **Vector Potential solver** — computes the magnetic vector potential **A** for wires, current loops, and dipoles.
- **Analytical fields** — closed-form models for Harris current sheets, magnetic dipoles, and circular current loops.

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
