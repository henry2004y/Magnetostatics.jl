# [Axisymmetric Reconstruction](@id inverse_axisymmetric)

[`reconstruct_axisymmetric_field`](@ref) fits a continuous field from cylindrical
probe measurements under azimuthal symmetry about the z axis. It interpolates
poloidal flux and fits the toroidal field independently; it does **not** solve
magnetohydrostatic pressure balance.

## Measurements and representation

- `positions` is an `n x 2` matrix with rows `(R, Z)` in metres.
- `B_observed` is an `n x 3` matrix of `(B_R, B_phi, B_Z)` in tesla.
  Non-finite field components are treated as missing measurements.
- `domain = (Rmin, Rmax, Zmin, Zmax)` must have `Rmin > 0`.
- `uncertainties` supplies component-wise magnetic uncertainties in tesla.

The poloidal field derives from a radial-basis expansion of flux and is
divergence-free by construction. `kernel = :wendland_c2` uses compact support;
`:gaussian` uses a global basis. `length_scale`, basis centers, and regularization
control the spatial model.

## Constraints and source-informed modes

- `boundary_positions` and `boundary_flux` specify exact poloidal-flux samples;
  `boundary_uncertainties` makes them weighted observations instead.
- `current_positions`, `current_density`, and `current_uncertainties` specify
  toroidal current density in A/m². Zero-current samples are an explicit vacuum
  prior, not a default physical assumption.
- `vacuum_degree` adds current-free polynomial modes for smooth distant-coil fields.
- `filament_positions` adds virtual external coaxial filaments with fitted currents.
  The filaments must lie outside the reconstruction domain.
- `current_element_positions` adds softened toroidal current elements inside the
  plasma region, with width controlled by `current_element_width`.

Compactly supported flux bumps carry no net toroidal current. For a distributed
plasma current, use current-carrying elements or assess whether a Gaussian basis
places the return current outside the domain. Merely increasing compact-basis
resolution does not supply a missing net-current mode.

## Evaluation and validation

The returned [`AxisymmetricFieldReconstruction`](@ref) is called at **Cartesian**
positions and returns Cartesian field components, despite using cylindrical
measurements for fitting. [`poloidal_flux`](@ref) and [`current_density`](@ref)
provide derived quantities. Evaluation outside the cylindrical domain raises
`DomainError`.

Use held-out probes and compare both field and toroidal current density. Vary
the basis, length scale, regularization, and source-informed modes; agreement
at probes alone does not establish the inferred plasma-current distribution.
Exact flux constraints use a dense null-space reduction, while uncertain
boundary data permit the sparse iterative path.

## Example and API

The [FLARE example](@ref flare_example) starts from known coils and then fits
sparse probes, comparing vacuum, filament, and plasma-current bases.
See [`reconstruct_axisymmetric_field`](@ref) for all keywords and
[Choosing an inverse method](@ref inverse_methods) for alternatives when
azimuthal symmetry is not appropriate.
