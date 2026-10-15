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
    return field.reference_field(point) +
        SVector{3}(map(row -> dot(row, field.coefficients), rows))
end

_reference_current(::NullField, pos) = zero(SVector{3, eltype(pos)})
_reference_current(::UniformField, pos) = zero(SVector{3, eltype(pos)})
_reference_current(field::AbstractMagneticField, pos) = current_density(field, pos)

function current_density(field::FieldReconstruction3D, pos)
    point = _check_3d_position(field, pos)
    _, rows, _ = _operators_3d(field, point)
    return _reference_current(field.reference_field, point) +
        SVector{3}(map(row -> dot(row, field.coefficients) / μ₀, rows))
end

function _matrix_3d(rows, ncoefficients, ::Type{T}) where {T}
    isempty(rows) && return spzeros(T, 0, ncoefficients)
    return sparse(reduce(hcat, rows)')
end

function _solve_3d(A, b, C, d, algorithm, abstol, reltol)
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
        LinearProblem(reduced, b - A * particular), algorithm; abstol, reltol
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

`regularization` is a positive coefficient ridge penalty. Sparse least squares
uses `algorithm = KrylovJL_LSMR()` with `abstol` and `reltol`. Sparse probes do
not determine a unique physical field; length scale, basis and regularization
are modeling assumptions, not additional measurements.
"""
function reconstruct_field_3d(
        positions::AbstractMatrix{<:Real}, B_observed::AbstractMatrix{<:Real};
        domain, length_scale, basis_centers = nothing, kernel = :wendland_c2,
        uncertainties = ones(size(B_observed)), regularization = 1.0e-8,
        reference_field::AbstractMagneticField = NullField(),
        algorithm = KrylovJL_LSMR(), abstol = 1.0e-10, reltol = 1.0e-10
    )
    length(domain) == 6 ||
        throw(ArgumentError("domain must contain three pairs of bounds"))
    T = promote_type(Float64, eltype(positions), eltype(B_observed), typeof(length_scale))
    box = NTuple{6, T}(Tuple(domain))
    all(isfinite, box) && all(box[2k] > box[2k - 1] for k in 1:3) ||
        throw(ArgumentError("domain bounds must be finite and strictly increasing"))
    for (name, value) in (
            ("length_scale", length_scale), ("regularization", regularization),
            ("abstol", abstol), ("reltol", reltol)
        )
        isfinite(value) && value > 0 ||
            throw(ArgumentError("$name must be finite and positive"))
    end
    kernel in (:wendland_c2, :gaussian) ||
        throw(ArgumentError("kernel must be :wendland_c2 or :gaussian"))
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
    uncovered = 0
    for i in axes(positions, 1)
        point = SVector{3, T}(positions[i, :])
        operators, _, count = _operators_3d(
            point, centers, T(length_scale), kernel, box, bins
        )
        uncovered += count == 0
        reference = reference_field(point)
        length(reference) == 3 && all(isfinite, reference) ||
            throw(ArgumentError("reference_field must return three finite components"))
        for k in 1:3
            target = B_observed[i, k]
            isnan(target) && continue
            isfinite(target) || throw(ArgumentError("observations cannot be infinite"))
            _append_fit_row!(
                rows, targets, sigmas, operators[k], T(target - reference[k]),
                T(uncertainties[i, k])
            )
        end
    end
    isempty(rows) && throw(ArgumentError("at least one finite observation is required"))
    H = _matrix_3d(rows, ncoefficients, T)
    weighted = Diagonal(inv.(sigmas)) * H
    A = vcat(weighted, sqrt(T(regularization)) * sparse(I, ncoefficients, ncoefficients))
    b = vcat(targets ./ sigmas, zeros(T, ncoefficients))
    coefficients, retcode, iterations = _solve_3d(
        A, b, spzeros(T, 0, ncoefficients), T[], algorithm, abstol, reltol
    )
    residual = H * coefficients - targets
    diagnostics = (
        weighted_rms = norm(residual ./ sigmas) / sqrt(length(targets)),
        observation_count = length(targets), solver_retcode = retcode,
        iterations = iterations, unsupported_probes = uncovered,
    )
    uncovered > 0 && @warn "RBF support does not cover all probes" unsupported = uncovered
    return FieldReconstruction3D(
        centers, coefficients, T(length_scale), box, kernel, bins,
        reference_field, diagnostics
    )
end
