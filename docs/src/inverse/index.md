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
The [uncertainty and hidden-topology benchmark](@ref reconstruction_uncertainty_example)
checks conditional noise propagation and demonstrates why it cannot establish
topological uniqueness outside the fitted model.

## Reconstructing topology from sparse measurements

The scientific target is often not just the field at unmeasured points, but its
**topology**: magnetic nulls, separatrices, flux surfaces, and connections between
regions or boundaries. A smooth divergence-free fit is a starting point, not
evidence that these structures have been recovered.

In an axisymmetric field, critical points of poloidal flux identify poloidal
O-points and X-points. With a nonzero toroidal field, these are not full magnetic
nulls. In general 3D fields, connectivity requires field-line analysis; a contour
plot or small vector-field residual is not a substitute. Near nulls and
separatrices, small field errors can change inferred connections substantially.

### What the present methods establish

The existing methods cover flexible divergence-free spatial fitting,
source-informed axisymmetric fitting, and planar equilibrium continuation.
They can generate candidate topologies under their respective assumptions.
They do not currently quantify how uniquely sparse measurements determine
those topologies.

Even complete normal-field data on a closed boundary do not determine an
arbitrary current-carrying interior field. A divergence-free perturbation
localized away from that boundary and the probes can preserve those observations
while changing the interior field. Finite current samples reduce this freedom
but do not remove it everywhere. Regularization chooses among admissible fields;
it does not demonstrate that their topologies agree.

The useful question is therefore: **which additional, independently justified
information rules out topologically different fields that fit the data?**

### Which additional inverse approaches would help?

The following are candidate capabilities, not currently available solver options.
Some are extensions to existing fits rather than separate reconstruction methods.

| Candidate | When it adds information | Gap relative to current methods | Main limitation |
|:--|:--|:--|:--|
| Known-source or reduced-order parameter inference | Coil geometry, possible current paths, or a credible family of fields is known | Fit a small set of physically meaningful currents or parameters instead of a flexible field; the axisymmetric filament/current-element bases already provide part of this approach | Unmodeled plasma currents or incorrect source geometry can bias the inferred topology |
| Current-free, boundary-adapted reconstruction | The reconstruction region is independently known to contain no current | Enforce `curl(B) = 0` throughout the region, rather than penalizing current or sampling zero current at selected points | Cannot represent current sheets or plasma currents inside that region; holes and linked external currents require circulation modes |
| Force-free reconstruction | Current is approximately parallel to the field and pressure forces are negligible | Enforce `curl(B) = alpha B` together with zero divergence; this couples field and current in a way the current spatial fits do not | Sparse probes and normal boundary data alone generally do not determine `alpha`; low pressure alone does not establish force-free balance |
| Axisymmetric equilibrium or magnetohydrostatic reconstruction | Static force balance and pressure/current-profile information are credible | Enforce force balance rather than only interpolate flux; cylindrical equilibrium reconstruction is distinct from the existing planar spacecraft GSR | Requires additional pressure/profile and boundary information; equilibrium assumptions can fail in a driven reconnection layer |
| Uncertainty-aware reconstruction and experimental design | Measurement errors and model assumptions can be characterized | Explore alternative fields consistent with the data, assess topology stability, and identify informative additional probes | Uncertainty is conditional on the chosen model and priors; a narrow posterior in an inadequate basis can be misleading |

For a current-free, simply connected region, write `B = -grad(Phi)` with
`laplacian(Phi) = 0`. Compatible normal-field data over the **entire** closed
boundary determine the field, with the potential fixed only up to a constant.
Sparse boundary samples do not provide that uniqueness. In multiply connected
regions, a single-valued scalar potential may omit fields with nonzero circulation;
the relevant global modes or linked-current information must also be supplied.

A conducting boundary is likewise not automatically a zero-normal-field boundary.
Its magnetic condition depends on the assumed flux history and material model.
Specify what is physically known, including external or far-field conditions,
rather than treating the numerical box as a measured wall.

Force-free and magnetohydrostatic approaches should be conditional additions,
not default upgrades. An equilibrium constraint is useful only where its
neglected forces and time dependence are demonstrably small. The
[force-free review by Wiegelmann and Sakurai](https://arxiv.org/abs/1208.4693)
discusses the assumptions and boundary-data issues of that class of models.

### How to decide whether another method is needed

Before adding a solver, compare the existing fits on synthetic problems designed
to challenge **topology**, not only pointwise accuracy:

1. Choose target observables: poloidal X/O-point locations, 3D null locations
   and types, separatrix geometry, or boundary-to-boundary connection classes.
   Define the domain and tracing termination rules for connectivity comparisons.
2. Use known forward fields with and without the target structures. Include
   cases where probes miss a localized current region or a separatrix, and cases
   that deliberately violate an assumed symmetry or equilibrium.
3. Hold out field, current, and boundary samples. Vary noise, probe placement,
   basis resolution, regularization, and plausible boundary assumptions.
4. Search for alternative fits that agree with the measurements within their
   errors but have different target structures. Parameter sweeps or bootstrap
   fits are useful sensitivity checks, not by themselves calibrated posterior
   probabilities.
5. Repeat with an independently justified source, vacuum, or equilibrium
   constraint. Check whether it improves topology recovery without hiding
   model mismatch or degrading held-out agreement.

The [existing reconstruction cases](@ref reconstruction_cases_example) provide
forward truths and held-out validation. The
[uncertainty benchmark](@ref reconstruction_uncertainty_example) additionally checks
null stability and an explicit topological ambiguity, but these controlled tests
do not establish identifiability for a general experiment. A topology can be considered
supported
only to the extent that it persists across the admissible models, measurement
errors, and numerical resolutions examined. If incompatible topologies remain,
report the ambiguity rather than selecting one from its fit residual alone.

### Recommended priorities

For sparse-measurement topology inference, the first priority is **topology
diagnostics and uncertainty/sensitivity analysis around the existing methods**.
This establishes whether the missing ingredient is information, representation,
or physics; another flexible interpolator does not resolve a data null space.

The next method should follow the experimental regime:

- **Known coils and uncertain source strengths:** extend source-informed inference,
  including general 3D configurations where axisymmetry is inappropriate.
- **Verified current-free measurement volume:** prioritize a harmonic,
  geometry-aware reconstruction with explicit global circulation constraints
  where needed.
- **Current-carrying plasma:** retain the divergence-free fit as a baseline, then
  assess force-free or pressure-balanced reconstruction only with supporting
  physical measurements.

Thus additional methods are warranted where they encode defensible information
that the existing representations do not enforce. No method can guarantee a
unique topology from arbitrary sparse probes and boundary samples alone.
