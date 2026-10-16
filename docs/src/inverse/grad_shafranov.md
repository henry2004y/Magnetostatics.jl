# [Spacecraft Grad-Shafranov Reconstruction](@id inverse_gsr)

[`reconstruct_grad_shafranov`](@ref) reconstructs a planar magnetohydrostatic
strip from a single straight spacecraft crossing. Unlike the spatial-probe
methods, it assumes isotropic pressure balance and translational invariance.
It is not the cylindrical axisymmetric flux fit.

## Measurements and frame

- `times` contains strictly increasing sample times in seconds.
- `B_observed` contains Cartesian magnetic components in tesla.
- `pressure` contains nonnegative isotropic thermal pressure in pascals.
  Missing or non-finite observations are not supported.
- `invariant_axis` is the supplied direction of translational invariance.
- `crossing_velocity` is the constant spacecraft velocity **relative to the
  structure**, in m/s, expressed in the measurement frame.
- `y` specifies the increasing transverse strip coordinates in metres,
  including zero and points on both sides of the track.

The solver does not estimate the de Hoffmann-Teller frame, search for the
invariant axis, or determine whether the equilibrium assumptions are valid.
Those inputs must be established independently.

## Representation and numerical controls

In the local frame,

```math
\mathbf B=(A_y,-A_x,B_z(A)),\qquad
P_t(A)=p(A)+\frac{B_z(A)^2}{2\mu_0},\qquad
\nabla^2 A=-\mu_0 P_t'(A).
```

The solver integrates flux along the track, fits `P_t(A)` and `B_z(A)` using
Chebyshev polynomials, and continues the planar equation in both transverse
directions. `pressure_degree` and `axial_degree` control the source-function
fits; `spatial_degree` truncates the along-track representation and acts as
spatial regularization. `max_step` limits internal continuation steps.

## Evaluation and validation

The returned [`GradShafranovReconstruction`](@ref) evaluates Cartesian positions
and returns Cartesian components in the original measurement frame.
[`flux_function`](@ref), [`thermal_pressure`](@ref), and [`current_density`](@ref)
provide derived quantities. Its finite domain is a strip in local x and y;
local z is unrestricted.

Elliptic Cauchy continuation is ill-posed. Keep the strip narrow and vary its
width, polynomial degrees, and step size. Inspect `field.diagnostics` for track
residuals, source-function extrapolation, and sampled GS residuals. Validate
off-track whenever possible: a small track residual does not guarantee stable
continuation. Warnings about strip width or flux extrapolation are not
quantitative stability or uncertainty bounds.

The interpolated field is divergence-free, but the GS equation is satisfied
only approximately. Analytic current is obtained from the interpolated field's
curl, not just from the fitted source function.

## Example and API

The [mock spacecraft crossing](@ref grad_shafranov_example) reconstructs a Fadeev
island and compares field, pressure, current, and continuation sensitivity
against off-track truth. See [`reconstruct_grad_shafranov`](@ref) for the complete
input contract and [Choosing an inverse method](@ref inverse_methods) for models
that do not assume pressure balance.
