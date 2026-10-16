# [Forward Field Construction](@id forward_methods)

The **forward problem** starts from a specified source, field model, or potential
and computes the magnetic field. Choose a representation that matches the known
geometry, current distribution, and boundary assumptions.

If the inputs are measured fields rather than known sources, start with
[Inverse reconstruction](@ref inverse_methods).

## Choosing a forward method

| Known input | Method | Assumptions and scope | Worked example |
|:--|:--|:--|:--|
| Parameters of a supported field model | Analytical fields such as [`Dipole`](@ref), [`CurrentLoopAnalytic`](@ref), and [`HarrisSheet`](@ref) | Geometry and physics are fixed by the selected model | [Dipole](@ref dipole_example), [current sheets](examples/current_sheet.md) |
| Discretized wire geometry and current | [`BiotSavart`](@ref) | Numerical line integration; accuracy depends on wire discretization | [Biot-Savart solver](@ref biot_savart_example) |
| Volumetric current on a uniform grid | [`FFTSolver`](@ref) | Periodic spectral solution by default; zero-padded free-space convolution with `OpenBoundary()` | [FFT solver](examples/fft.md) |
| Volumetric current on a Cartesian grid with prescribed vector-potential boundary values | [`PoissonSolver`](@ref) | Finite-difference solution of the vector-potential Poisson equation with Dirichlet boundaries | [Poisson solver](examples/poisson.md) |
| Wire, loop, or dipole source for a vector potential | [`VectorPotential`](@ref), [`vector_potential`](@ref) | Source-dependent potential representation | [Current loop](examples/current_loop.md), [Biot-Savart solver](@ref biot_savart_example) |

Boundary types are not interchangeable options for every solver. Check the
selected [`solve`](@ref) method before supplying a boundary. The
[Biot-Savart example](@ref biot_savart_example) covers image-source and surface-current
constructions for supported conducting boundaries.

## Basic workflow

1. Define the source geometry and current, or select an analytical field model.
2. Choose a solver consistent with the domain and boundary assumptions.
3. Evaluate the field at the required positions or on a grid.
4. Check discretization convergence and, where available, compare with an
   analytical solution before using the result as synthetic inverse data.

```@example forward_quickstart
using Magnetostatics, StaticArrays, LinearAlgebra

loop = CurrentLoop(1.0, 1.0, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
wire = discretize_loop(loop, 100)
position = SVector(0.0, 0.0, 0.5)
B = solve(BiotSavart(), wire, position)
B_exact = getB_loop(position, loop)
relative_error = norm(B - B_exact) / norm(B_exact)
```

Analytical field objects such as `CurrentLoopAnalytic(loop)` can instead be
called directly as `field(position)`. Numerical grid solvers have their own
input and output layouts; follow their worked examples rather than assuming
all solvers use pointwise evaluation.

## Example families

- **Sources and geometry:** [dipoles](examples/dipole.md),
  [current loops](examples/current_loop.md), [knotted wires](examples/trefoil_knot.md),
  and [magnetic mirrors](examples/magnetic_mirror.md).
- **Plasma and laboratory configurations:** [tokamaks](examples/tokamak.md),
  [Solov'ev equilibria](examples/solovev.md), [Taylor states](examples/taylor_state.md),
  [current sheets](examples/current_sheet.md), [FLARE coils](@ref flare_example),
  and [Z-pinches](examples/z_pinch.md).
- **Space and astrophysical fields:** [Parker spiral](examples/parker_spiral.md),
  [magnetospheres](examples/magnetosphere.md), and
  [empirical magnetospheres](examples/empirical_magnetosphere.md).

## From forward models to inverse validation

A forward model provides controlled synthetic observations and independent
validation points. Keep the generating model separate from the reconstruction:
only pass the data and physical constraints that would actually be available.
The [reconstruction cases](@ref reconstruction_cases_example) demonstrate this
separation across source and boundary configurations.
