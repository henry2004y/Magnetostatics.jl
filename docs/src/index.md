# Magnetostatics.jl

Magnetostatics.jl provides **forward construction** of magnetic fields from
specified sources or potentials, and **inverse reconstruction** from measurements.

## Start here

| Goal | Start with |
|:--|:--|
| Compute a field from known sources, currents, or model parameters | [Forward methods](@ref forward_methods) |
| Reconstruct a field from spatial probes or a spacecraft crossing | [Inverse methods](@ref inverse_methods) |
| Compare axisymmetric, general 3D, and spacecraft GSR assumptions | [Choosing an inverse method](@ref inverse_methods) |
| Compute field geometry or guiding-center drift diagnostics | [Magnetic geometry](examples/magnetic_geometry.md) |
| Look up types, keyword arguments, or supported operations | [API reference](api.md) |

Method guides explain assumptions, inputs, constraints, and validation.
Worked examples demonstrate complete physical configurations; the API reference
documents individual functions and types.

## Quick Start

```julia
using Magnetostatics
using StaticArrays

# Create a circular current loop and discretize it into wire segments
loop = CurrentLoop(1.0, 1.0, [0, 0, 0], [0, 0, 1])
wire = discretize_loop(loop, 100)

# Compute B at a point using the Biot-Savart solver
B = solve(BiotSavart(), wire, SVector(0.0, 0.0, 0.5))

# Compare with the analytical solution
B_exact = getB_loop(SVector(0.0, 0.0, 0.5), loop)
```

## Related Packages

- [FieldTracer.jl](https://github.com/henry2004y/FieldTracer.jl): Field-line tracing
  for analytical and reconstructed fields through a `(position, arc_length)` callback,
  and for numerical fields through component arrays on grids. See the
  [Parker spiral](@ref parker_spiral_example), [3D reconstruction](@ref reconstruction3d_example),
  and [FFT solver](examples/fft.md) demos, which use v0.2.4 or later.
  UniformStreamlines.jl is retained for automatic evenly spaced streamline
  placement in the other visualizations. The closed-form `dipole_fieldline`
  utility is an exact geometric construction, not a numerical tracer.
- [TestParticle.jl](https://github.com/henry2004y/TestParticle.jl): Charged test particle
  tracing and orbit dynamics (full-orbit and guiding-center pushers) in electromagnetic
  fields. Magnetic field models and numerical solutions from `Magnetostatics.jl` can be
  directly queried by `TestParticle.jl`.

## Contents

```@contents
Pages = [
    "forward.md",
    "inverse/index.md",
    "inverse/axisymmetric.md",
    "inverse/three_dimensional.md",
    "inverse/grad_shafranov.md",
    "api.md",
]
Depth = 2
```
