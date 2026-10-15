"""
    FieldReconstruction3D <: AbstractMagneticField

Continuous Cartesian field represented as a uniform field plus the analytic curl
of a radial-basis vector potential. The fitted field is divergence-free by
construction. Coordinates are in metres and field components in tesla.
Evaluation outside `domain` throws `DomainError`.
"""
struct FieldReconstruction3D{T, F <: AbstractMagneticField, D} <: AbstractMagneticField
    centers::Matrix{T}
    coefficients::Vector{T}
    length_scale::T
    domain::NTuple{6, T}
    kernel::Symbol
    center_bins::Dict{NTuple{3, Int}, Vector{Int}}
    reference_field::F
    diagnostics::D
end

function _validate_xyz(positions, domain, name)
    size(positions, 2) == 3 ||
        throw(ArgumentError("$name must have columns (x, y, z)"))
    all(isfinite, positions) || throw(ArgumentError("$name must be finite"))
    for axis in 1:3
        lo, hi = domain[2axis - 1], domain[2axis]
        all(lo .<= positions[:, axis] .<= hi) ||
            throw(ArgumentError("$name must lie inside domain"))
    end
    return nothing
end

function _check_3d_position(field, pos)
    length(pos) == 3 || throw(ArgumentError("position must have three coordinates"))
    all(
        field.domain[2k - 1] <= pos[k] <= field.domain[2k] for k in 1:3
    ) || throw(DomainError(pos, "position is outside the reconstruction domain"))
    return SVector{3}(pos)
end

function _centers_3d(domain, length_scale, ::Type{T}) where {T}
    coordinates = ntuple(3) do k
        lo, hi = domain[2k - 1], domain[2k]
        n = max(2, ceil(Int, sqrt(T(3)) * (hi - lo) / length_scale) + 1)
        range(lo, hi; length = n)
    end
    points = Iterators.product(coordinates...)
    centers = Matrix{T}(undef, length(points), 3)
    for (i, p) in enumerate(points)
        centers[i, :] .= p
    end
    return centers
end

function _bins_3d(centers, domain, scale)
    bins = Dict{NTuple{3, Int}, Vector{Int}}()
    for i in axes(centers, 1)
        cell = ntuple(k -> floor(Int, (centers[i, k] - domain[2k - 1]) / scale), 3)
        push!(get!(bins, cell, Int[]), i)
    end
    return bins
end

function _candidates_3d(pos, domain, scale, bins)
    cell = ntuple(k -> floor(Int, (pos[k] - domain[2k - 1]) / scale), 3)
    candidates = Int[]
    for offset in Iterators.product(-1:1, -1:1, -1:1)
        append!(candidates, get(bins, cell .+ offset, Int[]))
    end
    return candidates
end

function _kernel_derivatives_3d(displacement, scale, kernel)
    T = eltype(displacement)
    r = norm(displacement)
    inv_scale2 = inv(scale^2)
    identity = one(SMatrix{3, 3, T, 9})
    if kernel === :gaussian
        value = exp(-r^2 * inv_scale2 / 2)
        gradient = -value * inv_scale2 * displacement
        hessian = value * (
            inv_scale2^2 * displacement * displacement' - inv_scale2 * identity
        )
    else
        q = r / scale
        q >= 1 && return zero(displacement), zero(identity)
        tangential = -20 * (1 - q)^3 * inv_scale2
        gradient = tangential * displacement
        hessian = tangential * identity
        if r > 0
            hessian += 60 * (1 - q)^2 * inv_scale2 / (scale * r) *
                displacement * displacement'
        end
    end
    return gradient, hessian
end

function _operators_3d(pos, centers, scale, kernel, domain, bins)
    T = promote_type(eltype(pos), eltype(centers), typeof(scale))
    ncoefficients = 3 + 3size(centers, 1)
    field_rows = [spzeros(T, ncoefficients) for _ in 1:3]
    curl_rows = [spzeros(T, ncoefficients) for _ in 1:3]
    for k in 1:3
        field_rows[k][k] = one(T)
    end
    candidates = kernel === :gaussian ? axes(centers, 1) :
        _candidates_3d(pos, domain, scale, bins)
    support_count = 0
    for j in candidates
        displacement = SVector{3, T}(pos) - SVector{3, T}(centers[j, :])
        kernel === :wendland_c2 && norm(displacement) >= scale && continue
        support_count += 1
        gradient, hessian = _kernel_derivatives_3d(displacement, scale, kernel)
        # Scaling A by the support radius gives all coefficients field units.
        curl_block = scale * (hessian - tr(hessian) * one(hessian))
        for k in 1:3
            direction = SVector{3, T}(ntuple(i -> T(i == k), 3))
            field_column = scale * cross(gradient, direction)
            column = 3 + 3(j - 1) + k
            for i in 1:3
                field_rows[i][column] = field_column[i]
                curl_rows[i][column] = curl_block[i, k]
            end
        end
    end
    return field_rows, curl_rows, support_count
end

function _operators_3d(field::FieldReconstruction3D, pos)
    return _operators_3d(
        pos, field.centers, field.length_scale, field.kernel,
        field.domain, field.center_bins
    )
end

function (field::FieldReconstruction3D)(pos)
    point = _check_3d_position(field, pos)
    rows, _, _ = _operators_3d(field, point)
    return _reference_value_3d(field.reference_field, point) +
        SVector{3}(map(row -> dot(row, field.coefficients), rows))
end

_reference_current(::NullField, pos) = zero(SVector{3, eltype(pos)})
_reference_current(::UniformField, pos) = zero(SVector{3, eltype(pos)})
_reference_current(field::AbstractMagneticField, pos) = current_density(field, pos)

function current_density(field::FieldReconstruction3D, pos)
    point = _check_3d_position(field, pos)
    _, rows, _ = _operators_3d(field, point)
    return _reference_value_3d(field.reference_field, point; current = true) +
        SVector{3}(map(row -> dot(row, field.coefficients) / μ₀, rows))
end

function _matrix_3d(rows, ncoefficients, ::Type{T}) where {T}
    isempty(rows) && return spzeros(T, 0, ncoefficients)
    return sparse(reduce(hcat, rows)')
end

function _reference_value_3d(field, point; current = false)
    value = current ? _reference_current(field, point) : field(point)
    length(value) == 3 && all(isfinite, value) ||
        throw(ArgumentError("reference field must return three finite components"))
    return SVector{3}(value)
end

function _add_observation_3d!(rows, targets, sigmas, groups, row, target, sigma, group)
    isnan(target) && return nothing
    isfinite(target) || throw(ArgumentError("observations cannot be infinite"))
    _append_fit_row!(rows, targets, sigmas, row, target, sigma)
    push!(groups, group)
    return nothing
end

function _validate_directions_3d(directions, positions, name)
    size(directions) == size(positions) ||
        throw(ArgumentError("$name must have the same size as its positions"))
    all(isfinite, directions) || throw(ArgumentError("$name must be finite"))
    unit_rows = all(
        isapprox(norm(row), 1; atol = 1.0e-8, rtol = 1.0e-8) for row in eachrow(directions)
    )
    unit_rows ||
        throw(ArgumentError("$name must contain unit vectors"))
    return nothing
end

function _quadrature_3d(domain, grid, ::Type{T}) where {T}
    length(grid) == 3 && all(n isa Integer && n >= 2 for n in grid) ||
        throw(ArgumentError("regularization_grid must contain three integers >= 2"))
    coordinates = ntuple(k -> range(domain[2k - 1], domain[2k]; length = grid[k]), 3)
    points = Matrix{T}(undef, prod(grid), 3)
    weights = Vector{T}(undef, prod(grid))
    for (i, index) in enumerate(CartesianIndices(Tuple(grid)))
        points[i, :] .= ntuple(k -> coordinates[k][index[k]], 3)
        weights[i] = prod(1:3) do k
            endpoint = index[k] == 1 || index[k] == grid[k]
            (endpoint ? T(0.5) : one(T)) / (grid[k] - 1)
        end
    end
    return points, weights
end

function _regularizer_3d(
        centers, scale, kernel, domain, bins, grid, mode, ::Type{T}
    ) where {T}
    ncoefficients = 3 + 3size(centers, 1)
    points, weights = _quadrature_3d(domain, grid, T)
    rows = SparseVector{T, Int}[]
    unsupported = 0
    for i in axes(points, 1)
        field_rows, curl_rows, count = _operators_3d(
            SVector{3, T}(points[i, :]), centers, scale, kernel, domain, bins
        )
        unsupported += count == 0
        if mode === :field
            factor = sqrt(weights[i])
            append!(rows, factor .* field_rows)
            append!(rows, (factor * scale) .* curl_rows)
        end
    end
    identity = sparse(T(1) * I, ncoefficients, ncoefficients)
    L = mode === :coefficients ? identity :
        vcat(_matrix_3d(rows, ncoefficients, T), sqrt(eps(T)) * identity)
    return L, unsupported, size(points, 1)
end

function _group_rms_3d(residual, groups, group)
    selected = residual[groups .== group]
    return isempty(selected) ? nothing : norm(selected) / sqrt(length(selected))
end

function _solve_3d(A, b, C, d, algorithm, abstol, reltol, maxiters)
    T = eltype(b)
    particular = zeros(T, size(A, 2))
    basis = nothing
    if !isempty(d)
        dense_C = Matrix(C)
        particular = pinv(dense_C) * d
        tolerance = 100eps(T) * max(one(T), norm(d, Inf))
        norm(C * particular - d, Inf) <= tolerance ||
            throw(ArgumentError("exact boundary constraints are inconsistent"))
        basis = nullspace(dense_C)
        size(basis, 2) == 0 && return particular, :constrained, 0
    end
    reduced = basis === nothing ? A : A * basis
    solution = LinearSolve.solve(
        LinearProblem(reduced, b - A * particular), algorithm; abstol, reltol, maxiters
    )
    LinearSolve.SciMLBase.successful_retcode(solution) ||
        error("3D reconstruction did not converge: $(solution.retcode)")
    coefficients = basis === nothing ? solution.u : particular + basis * solution.u
    all(isfinite, coefficients) || error("3D reconstruction produced non-finite values")
    return coefficients, solution.retcode, solution.iters
end

"""
    reconstruct_field_3d(positions, B_observed; domain, length_scale, ...)

Fit a divergence-free Cartesian magnetic field. `positions` and `B_observed`
are `n × 3` matrices. `NaN` field components denote missing measurements;
infinite values are rejected. `uncertainties` has the observation shape and is
in tesla. The box `domain` is `(xmin, xmax, ymin, ymax, zmin, zmax)`.

The model is `reference_field + b₀ + curl(A)`, with
`A = length_scale * sum(φⱼ * aⱼ)`. `:wendland_c2` (default) has support radius
`length_scale`; `:gaussian` is global. Default centers cover the box; custom
`basis_centers` must be an `m × 3` matrix inside it. `reference_field` defaults
to `NullField()` and must itself be divergence-free.

Additional directional probes use `directional_positions`, unit
`directional_directions` and scalar `directional_values` in tesla, with
`directional_uncertainties` (default: one tesla). Cartesian or directional
components marked `NaN` are skipped, including their uncertainty entries.

`boundary_positions`, unit `boundary_normals` and finite `boundary_values`
constrain `n ⋅ B`. Without `boundary_uncertainties` these are exact constraints
at the supplied points; otherwise they are weighted observations in tesla.
The evaluation box need not coincide with this physical surface. Exact
constraints use a dense null-space reduction, intended for modest problems.
They do not enforce a condition between samples or establish closed-surface
flux compatibility. Prefer soft constraints for noisy surface data.

`current_positions` and `current_density` are `m × 3` matrices; current
components and `current_uncertainties` (default: one A/m²) are in A/m².
Current observations and `current_density(field, pos)` require the reference
field to implement current density (zero is supported for null/uniform fields).

`regularization` is positive. With `regularization_mode = :field` (default),
it penalizes the box-averaged `|δB|² + length_scale² * |curl(δB)|²`, where
`δB = B - reference_field`, using normalized tensor trapezoidal quadrature.
The grid defaults to at least two intervals per length scale along each axis;
override it with `regularization_grid = (nx, ny, nz)`. A machine-epsilon
coefficient ridge stabilizes potential gauge modes. `:coefficients` instead
uses a coefficient ridge alone. All coefficients have field units.

Sparse least squares uses `algorithm = KrylovJL_LSMR()` with `abstol`,
`reltol`, and `maxiters` (default: 10,000). Nonconvergence throws an error.
`field.diagnostics` reports residuals in observation units, weighted
RMS, exact-constraint error, convergence and sampled RBF support coverage.
Uncovered locations retain only the uniform and reference fields and generate
a fitting warning. Coverage on the quadrature grid is not a proof of coverage
everywhere, nor are fit residuals uncertainty estimates. Sparse probes do not
determine a unique physical field; basis, length scale and regularization are
modeling assumptions, not additional measurements.
"""
function reconstruct_field_3d(
        positions::AbstractMatrix{<:Real}, B_observed::AbstractMatrix{<:Real};
        domain, length_scale, basis_centers = nothing, kernel = :wendland_c2,
        uncertainties = ones(size(B_observed)), regularization = 1.0e-8,
        reference_field::AbstractMagneticField = NullField(),
        algorithm = KrylovJL_LSMR(), abstol = 1.0e-10, reltol = 1.0e-10,
        maxiters::Integer = 10_000,
        regularization_mode = :field, regularization_grid = nothing,
        directional_positions = nothing, directional_directions = nothing,
        directional_values = nothing, directional_uncertainties = nothing,
        boundary_positions = nothing, boundary_normals = nothing,
        boundary_values = nothing, boundary_uncertainties = nothing,
        current_positions = nothing, current_density = nothing,
        current_uncertainties = nothing
    )
    length(domain) == 6 ||
        throw(ArgumentError("domain must contain three pairs of bounds"))
    T = promote_type(Float64, eltype(positions), eltype(B_observed), typeof(length_scale))
    box = NTuple{6, T}(Tuple(domain))
    all(isfinite, box) && all(box[2k] > box[2k - 1] for k in 1:3) ||
        throw(ArgumentError("domain bounds must be finite and strictly increasing"))
    for (name, value) in (
            ("length_scale", length_scale), ("regularization", regularization),
            ("abstol", abstol), ("reltol", reltol),
        )
        isfinite(value) && value > 0 ||
            throw(ArgumentError("$name must be finite and positive"))
    end
    kernel in (:wendland_c2, :gaussian) ||
        throw(ArgumentError("kernel must be :wendland_c2 or :gaussian"))
    regularization_mode in (:field, :coefficients) ||
        throw(ArgumentError("regularization_mode must be :field or :coefficients"))
    maxiters > 0 || throw(ArgumentError("maxiters must be positive"))
    _validate_xyz(positions, box, "positions")
    size(B_observed) == size(positions) ||
        throw(ArgumentError("B_observed must have the same size as positions"))
    size(uncertainties) == size(B_observed) ||
        throw(ArgumentError("uncertainties must have the same size as B_observed"))
    centers = basis_centers === nothing ? _centers_3d(box, T(length_scale), T) :
        Matrix{T}(basis_centers)
    _validate_xyz(centers, box, "basis_centers")
    size(centers, 1) > 0 || throw(ArgumentError("at least one center is required"))
    bins = _bins_3d(centers, box, T(length_scale))
    ncoefficients = 3 + 3size(centers, 1)
    rows, targets, sigmas = SparseVector{T, Int}[], T[], T[]
    groups = Symbol[]
    exact_rows, exact_targets = SparseVector{T, Int}[], T[]
    uncovered = 0
    for i in axes(positions, 1)
        point = SVector{3, T}(positions[i, :])
        operators, _, count = _operators_3d(
            point, centers, T(length_scale), kernel, box, bins
        )
        uncovered += count == 0
        reference = _reference_value_3d(reference_field, point)
        for k in 1:3
            target = B_observed[i, k]
            _add_observation_3d!(
                rows, targets, sigmas, groups, operators[k], T(target - reference[k]),
                T(uncertainties[i, k]), :probe
            )
        end
    end
    for (locations, directions, values, errors, group) in (
            (
                directional_positions, directional_directions, directional_values,
                directional_uncertainties, :directional,
            ),
            (
                boundary_positions, boundary_normals, boundary_values,
                boundary_uncertainties, :boundary,
            ),
        )
        supplied = (locations !== nothing, directions !== nothing, values !== nothing)
        if !any(supplied)
            errors === nothing ||
                throw(ArgumentError("$group uncertainties require observations"))
            continue
        end
        all(supplied) ||
            throw(ArgumentError("$group positions, directions and values are required"))
        _validate_xyz(locations, box, "$group positions")
        _validate_directions_3d(directions, locations, "$group directions")
        values isa AbstractVector && length(values) == size(locations, 1) ||
            throw(ArgumentError("$group values must have one entry per position"))
        group === :boundary && !all(isfinite, values) &&
            throw(ArgumentError("boundary values must be finite"))
        exact = group === :boundary && errors === nothing
        matching_errors = errors isa AbstractVector && length(errors) == length(values)
        errors === nothing || matching_errors ||
            throw(ArgumentError("$group uncertainties must match values"))
        for i in eachindex(values)
            point = SVector{3, T}(locations[i, :])
            direction = SVector{3, T}(directions[i, :])
            operators, _, count = _operators_3d(
                point, centers, T(length_scale), kernel, box, bins
            )
            uncovered += count == 0
            row = sum(direction[k] * operators[k] for k in 1:3)
            reference = _reference_value_3d(reference_field, point)
            target = T(values[i] - dot(direction, reference))
            if exact
                push!(exact_rows, row)
                push!(exact_targets, target)
            else
                sigma = errors === nothing ? one(T) : T(errors[i])
                _add_observation_3d!(
                    rows, targets, sigmas, groups, row, target, sigma, group
                )
            end
        end
    end
    if current_positions !== nothing || current_density !== nothing
        current_positions !== nothing && current_density !== nothing ||
            throw(ArgumentError("current_positions and current_density are required"))
        _validate_xyz(current_positions, box, "current_positions")
        size(current_density) == size(current_positions) ||
            throw(ArgumentError("current_density must match current_positions"))
        errors = current_uncertainties === nothing ? ones(T, size(current_density)) :
            current_uncertainties
        size(errors) == size(current_density) ||
            throw(ArgumentError("current_uncertainties must match current_density"))
        for i in axes(current_positions, 1)
            point = SVector{3, T}(current_positions[i, :])
            _, operators, count = _operators_3d(
                point, centers, T(length_scale), kernel, box, bins
            )
            uncovered += count == 0
            reference = _reference_value_3d(reference_field, point; current = true)
            for k in 1:3
                _add_observation_3d!(
                    rows, targets, sigmas, groups, operators[k] / T(μ₀),
                    T(current_density[i, k] - reference[k]), T(errors[i, k]), :current
                )
            end
        end
    elseif current_uncertainties !== nothing
        throw(ArgumentError("current_uncertainties require current observations"))
    end
    isempty(rows) && isempty(exact_rows) &&
        throw(ArgumentError("a finite observation or exact constraint is required"))
    grid = regularization_grid === nothing ?
        ntuple(k -> max(3, ceil(Int, 2(box[2k] - box[2k - 1]) / length_scale) + 1), 3) :
        regularization_grid
    L, unsupported_grid, grid_count = _regularizer_3d(
        centers, T(length_scale), kernel, box, bins, grid, regularization_mode, T
    )
    H = _matrix_3d(rows, ncoefficients, T)
    weighted = Diagonal(inv.(sigmas)) * H
    A = vcat(weighted, sqrt(T(regularization)) * L)
    b = vcat(targets ./ sigmas, zeros(T, size(L, 1)))
    C = _matrix_3d(exact_rows, ncoefficients, T)
    coefficients, retcode, iterations = _solve_3d(
        A, b, C, exact_targets, algorithm, abstol, reltol, maxiters
    )
    residual = H * coefficients - targets
    exact_error = isempty(exact_targets) ? zero(T) :
        norm(C * coefficients - exact_targets, Inf)
    tolerance = 100eps(T) * max(one(T), norm(exact_targets))
    exact_error <= tolerance || error("exact boundary constraints lost numerical accuracy")
    diagnostics = (
        weighted_rms = isempty(targets) ? nothing :
            norm(residual ./ sigmas) / sqrt(length(targets)),
        probe_rms = _group_rms_3d(residual, groups, :probe),
        directional_rms = _group_rms_3d(residual, groups, :directional),
        boundary_rms = _group_rms_3d(residual, groups, :boundary),
        current_rms = _group_rms_3d(residual, groups, :current),
        exact_boundary_max_error = exact_error,
        exact_constraint_count = length(exact_targets),
        coefficient_count = ncoefficients,
        observation_count = length(targets), solver_retcode = retcode,
        iterations = iterations, unsupported_observations = uncovered,
        unsupported_grid_points = unsupported_grid, grid_point_count = grid_count,
        sampled_support_fraction = 1 - unsupported_grid / grid_count,
    )
    if uncovered > 0 || unsupported_grid > 0
        message = "RBF support has gaps; only uniform/reference fields act there"
        @warn message unsupported_observations = uncovered grid_gaps = unsupported_grid
    end
    return FieldReconstruction3D(
        centers, coefficients, T(length_scale), box, kernel, bins,
        reference_field, diagnostics
    )
end
