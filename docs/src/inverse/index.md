# [Inverse Field Reconstruction](@id inverse_methods)

The **inverse problem** starts from measurements and reconstructs a field
consistent with those data and explicit physical assumptions. It is not simply
a forward solver run backwards: sparse data generally do not determine a unique
field or current distribution.

## Choosing an inverse method

The package currently provides three physically distinct models:

| Method | Required measurements | Physical assumption | Guide and example |
|:--|:--|:--|:--|
| [`reconstruct_axisymmetric_field`](@ref) | Spatial probes at `(R, Z)` with cylindrical `(B_R, B_phi, B_Z)` components | Azimuthal symmetry about the z axis; pressure balance is not enforced | [Axisymmetric guide](@ref inverse_axisymmetric), [FLARE](@ref flare_example) |
| [`reconstruct_field_3d`](@ref) | Cartesian vector and/or directional field probes | No spatial symmetry; a divergence-free representation, not pressure balance | [3D guide](@ref inverse_3d), [tilted coil](@ref reconstruction3d_example) |
| [`reconstruct_grad_shafranov`](@ref) | Magnetic-field and isotropic thermal-pressure time series, supplied invariant axis and crossing velocity | Planar translational invariance and magnetohydrostatic balance | [Spacecraft GSR guide](@ref inverse_gsr), [mock crossing](@ref grad_shafranov_example) |

An axisymmetric fit is not a cheaper general 3D model: applying it to a tilted
coil or non-axisymmetric boundary discards real structure. Spacecraft GSR is
not cylindrical Grad-Shafranov equilibrium reconstruction; it uses the planar
Laplacian and requires pressure data and a supplied crossing frame.

### Spatial-probe constraints

| | Axisymmetric | Three-dimensional |
|:--|:--|:--|
| Domain | `(Rmin, Rmax, Zmin, Zmax)`, with `Rmin > 0` | `(xmin, xmax, ymin, ymax, zmin, zmax)` |
| Representation | Poloidal-flux RBFs and optional vacuum, filament, and current-element modes; independent toroidal fit | Uniform field plus curl of a vector-potential RBF expansion |
| Boundary data | Poloidal **flux**, `boundary_flux` | Normal **field**, `boundary_values = n dot B` |
| Current data | Toroidal `J_phi` | Cartesian `(J_x, J_y, J_z)` |
| Background treatment | Vacuum/filament modes fitted from data | Optional fixed `reference_field` |

Both spatial-probe methods accept exact boundary samples or uncertainty-weighted
boundary observations. Exactness applies at supplied sample locations, not to
the entire boundary. Neither infers a wall condition from the reconstruction box.

## Shared workflow

1. **Select the physics first.** Establish symmetry, equilibrium validity, domain,
   and the measurements actually available before choosing a solver.
2. **Prepare consistent data.** Use metres, tesla, and the coordinate/component
   conventions of the selected method. Supply calibrated uncertainties where
   supported; default unit weights are not measurement-noise estimates.
3. **State constraints and priors.** Distinguish measured boundary/current data
   from assumed constraints, known backgrounds, and regularization.
4. **Fit and evaluate.** All three methods return callable fields evaluated at
   Cartesian positions in metres, returning Cartesian magnetic components in tesla.
   All support [`current_density`](@ref); additional quantities are method-specific.
5. **Validate independently.** Check held-out data, physical residuals, and
   sensitivity to basis resolution, regularization, and domain extent.

## Validation and interpretation

- Small training residuals do not establish uniqueness or off-probe accuracy.
- Check current density as well as magnetic field: differentiation can amplify
  poorly constrained structure.
- Compare with held-out probes and boundary locations that were not used to fit.
- Repeat fits with plausible priors and numerical resolutions. For GSR, also vary
  continuation width, spatial degree, and step size.
- Inspect available diagnostics and warnings. The 3D and GSR models expose
  `field.diagnostics`; these residuals are not posterior uncertainty estimates.
- Keep synthetic truth out of the fit unless it represents genuinely known
  information, such as an explicitly supplied background field.

See [reconstruction performance](@ref reconstruction_performance) for 3D
accuracy/cost measurements, and [physical reconstruction cases](@ref reconstruction_cases_example)
for validation across sources and boundaries.

## Adding inverse methods

New inverse methods should have their own method guide rather than being presented
as interchangeable options of an existing solver. Document:

1. Required measurements, units, coordinate frames, and missing-data support.
2. Physical assumptions, representation, domain, and identifiable quantities.
3. Boundary/current constraints, priors, regularization, and numerical controls.
4. Returned field evaluation, supported derived quantities, and diagnostics.
5. A reproducible example with independent validation and stated failure modes.

Add each method to the comparison above, the inverse-method navigation, and the
inverse API reference. Keep method explanations separate from worked examples
and API docstrings so the documentation can grow without a single solver catalogue.
