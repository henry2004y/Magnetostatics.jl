# Magnetostatics.jl

Magnetostatics.jl provides solvers and analytical models for computing magnetic fields produced by steady-state current distributions or potentials.

## Methods

### Analytical Models
Exact solutions for specific geometries:
- **`HarrisSheet`**: Harris current sheet model ($B_x(z) = B_0 \tanh(z/L)$).
- **`AsymmetricHarrisSheet`**: Asymmetric current sheet ($B_x(z) = \frac{B_1+B_2}{2} + \frac{B_1-B_2}{2} \tanh(z/L)$).
- **`ForceFreeHarrisSheet`**: Force-free sheet ($B_x(z) = B_0 \tanh(z/L)$, $B_y(z) = B_0 \text{sech}(z/L)$).
- **`BifurcatedHarrisSheet`**: Double-peaked current sheet ($B_x(z) = \frac{B_0}{2} [\tanh((z-d)/L) + \tanh((z+d)/L)]$).
- **`FadeevIsland`**: 2D magnetic island equilibrium ($A_y(x, z) = -B_0 L \ln[\cosh(z/L) + \epsilon \cos(x/L_x)]$).
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
- **`PoissonSolver`**: Finite-difference Poisson solver (∇²**A** = −μ₀**J**) on a Cartesian grid with Dirichlet boundary conditions, supporting anisotropic grid spacing via `LinearSolve.jl`.

### Field Reconstruction
- **[`reconstruct_field_3d`](@ref)**: Continuous divergence-free Cartesian field
  fitting from sparse vector or directional probes, with exact/weighted normal
  boundary samples, current measurements, and field-energy regularization.
  See the [API guide](api.md#Three-dimensional-reconstruction) for an example
  and the limitations of sparse-probe inference.
  The [3D reconstruction demo](@ref reconstruction3d_example) includes noisy
  probes, boundary measurements, visual comparisons, and held-out validation.

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
Pages = ["examples.md", "api.md"]
```
