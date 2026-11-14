"""
    FieldReconstruction3D <: AbstractMagneticField

Continuous Cartesian field represented as a uniform field plus the analytic curl
of a radial-basis vector potential. The fitted field is divergence-free by
construction. Coordinates are in metres and field components in tesla.
Evaluation outside `domain` throws `DomainError`.
Unlike [`AxisymmetricFieldReconstruction`](@ref), this model makes no symmetry
assumption and fits Cartesian probes in a three-dimensional box.
"""
struct _GridBins3D
    origin::NTuple{3, Int}
    dims::NTuple{3, Int}
    offsets::Vector{Int}
    indices::Vector{Int}
end

struct FieldReconstruction3D{T, F <: AbstractMagneticField, D} <: AbstractMagneticField
    centers::Matrix{T}
    coefficients::Vector{T}
    length_scale::T
    domain::NTuple{6, T}
    kernel::Symbol
    center_bins::Dict{NTuple{3, Int}, Vector{Int}}
    grid_bins::_GridBins3D
    reference_field::F
    diagnostics::D
end

function FieldReconstruction3D(
        centers::Matrix{T}, coefficients::Vector{T}, length_scale::T,
        domain::NTuple{6, T}, kernel::Symbol,
        center_bins::Dict{NTuple{3, Int}, Vector{Int}},
        reference_field::F, diagnostics::D
    ) where {T, F <: AbstractMagneticField, D}
    return FieldReconstruction3D{T, F, D}(
        centers, coefficients, length_scale, domain, kernel,
        center_bins, _to_grid_bins_3d(center_bins), reference_field, diagnostics
    )
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

function _to_grid_bins_3d(bins::Dict{NTuple{3, Int}, Vector{Int}})
    isempty(bins) && return _GridBins3D((0, 0, 0), (0, 0, 0), Int[1], Int[])
    min_x = typemax(Int); min_y = typemax(Int); min_z = typemax(Int)
    max_x = typemin(Int); max_y = typemin(Int); max_z = typemin(Int)
    total = 0
    for (cell, vec) in bins
        cx, cy, cz = cell
        cx < min_x && (min_x = cx); cx > max_x && (max_x = cx)
        cy < min_y && (min_y = cy); cy > max_y && (max_y = cy)
        cz < min_z && (min_z = cz); cz > max_z && (max_z = cz)
        total += length(vec)
    end
    origin = (min_x, min_y, min_z)
    dims = (max_x - min_x + 1, max_y - min_y + 1, max_z - min_z + 1)
    nx, ny, _ = dims
    nx_ny = nx * ny
    ncells = nx_ny * dims[3]
    counts = zeros(Int, ncells)
    for (cell, vec) in bins
        cidx = (cell[3] - min_z) * nx_ny + (cell[2] - min_y) * nx + (cell[1] - min_x) + 1
        @inbounds counts[cidx] += length(vec)
    end
    offsets = Vector{Int}(undef, ncells + 1)
    offsets[1] = 1
    @inbounds for c in 1:ncells
        offsets[c + 1] = offsets[c] + counts[c]
    end
    indices = Vector{Int}(undef, total)
    fill!(counts, 0)
    for (cell, vec) in bins
        cidx = (cell[3] - min_z) * nx_ny + (cell[2] - min_y) * nx + (cell[1] - min_x) + 1
        pos = offsets[cidx] + counts[cidx]
        @inbounds for v in vec
            indices[pos] = v
            pos += 1
        end
        counts[cidx] += length(vec)
    end
    return _GridBins3D(origin, dims, offsets, indices)
end

@inline function _check_3d_position(field, pos)
    length(pos) == 3 || throw(ArgumentError("position must have three coordinates"))
    @inbounds begin
        d = field.domain
        (d[1] <= pos[1] <= d[2] && d[3] <= pos[2] <= d[4] && d[5] <= pos[3] <= d[6]) ||
            throw(DomainError(pos, "position is outside the reconstruction domain"))
        return SVector{3}(pos[1], pos[2], pos[3])
    end
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

function _foreach_kernel_3d(
        f::F, pos, centers, scale, kernel, domain, bins
    ) where {F}
    if kernel === :gaussian
        for j in axes(centers, 1)
            _kernel_term_3d(f, j, pos, centers, scale, kernel)
        end
        return nothing
    end
    cell = ntuple(k -> floor(Int, (pos[k] - domain[2k - 1]) / scale), 3)
    for offset in Iterators.product(-1:1, -1:1, -1:1)
        bin = get(bins, cell .+ offset, nothing)
        bin === nothing && continue
        for j in bin
            _kernel_term_3d(f, j, pos, centers, scale, kernel)
        end
    end
    return nothing
end

@inline function _kernel_term_3d(f::F, j, pos, centers, scale, kernel) where {F}
    T = promote_type(eltype(pos), eltype(centers), typeof(scale))
    displacement = SVector{3, T}(pos) -
        SVector{3, T}(ntuple(k -> centers[j, k], 3))
    kernel === :wendland_c2 && norm(displacement) >= scale && return nothing
    gradient, hessian = _kernel_derivatives_3d(displacement, scale, kernel)
    f(j, gradient, hessian)
    return nothing
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
    field_indices = ntuple(k -> [k], 3)
    field_values = ntuple(_ -> T[1], 3)
    curl_indices = ntuple(_ -> Int[], 3)
    curl_values = ntuple(_ -> T[], 3)
    support_count = Ref(0)
    _foreach_kernel_3d(pos, centers, scale, kernel, domain, bins) do j, gradient, hessian
        support_count[] += 1
        # Scaling A by the support radius gives all coefficients field units.
        curl_block = scale * (hessian - tr(hessian) * one(hessian))
        for k in 1:3
            direction = SVector{3, T}(ntuple(i -> T(i == k), 3))
            field_column = scale * cross(gradient, direction)
            column = 3 + 3(j - 1) + k
            for i in 1:3
                if !iszero(field_column[i])
                    push!(field_indices[i], column)
                    push!(field_values[i], field_column[i])
                end
                if !iszero(curl_block[i, k])
                    push!(curl_indices[i], column)
                    push!(curl_values[i], curl_block[i, k])
                end
            end
        end
    end
    field_rows = [sparsevec(field_indices[k], field_values[k], ncoefficients) for k in 1:3]
    curl_rows = [sparsevec(curl_indices[k], curl_values[k], ncoefficients) for k in 1:3]
    return field_rows, curl_rows, support_count[]
end

function _operators_3d(field::FieldReconstruction3D, pos)
    return _operators_3d(
        pos, field.centers, field.length_scale, field.kernel,
        field.domain, field.center_bins
    )
end

function (field::FieldReconstruction3D)(pos)
    point = _check_3d_position(field, pos)
    return _reference_value_3d(field.reference_field, point) +
        _correction_3d(field, point, Val(false))
end

@inline function _correction_3d(
        field::FieldReconstruction3D, point::SVector{3}, ::Val{current}
    ) where {current}
    T = promote_type(eltype(point), eltype(field.coefficients))
    coefficients = field.coefficients
    centers = field.centers
    scale = field.length_scale
    kernel = field.kernel
    px, py, pz = T(point[1]), T(point[2]), T(point[3])
    pvec = SVector{3, T}(px, py, pz)
    @inbounds acc = current ? zero(SVector{3, T}) :
        SVector{3, T}(T(coefficients[1]), T(coefficients[2]), T(coefficients[3]))
    if kernel === :gaussian
        @inbounds for j in axes(centers, 1)
            cj = SVector{3, T}(T(centers[j, 1]), T(centers[j, 2]), T(centers[j, 3]))
            disp = pvec - cj
            gradient, hessian = _kernel_derivatives_3d(disp, scale, kernel)
            off = 3 + 3 * (j - 1)
            a = SVector{3, T}(
                T(coefficients[off + 1]), T(coefficients[off + 2]), T(coefficients[off + 3])
            )
            acc += current ? scale * (hessian * a - tr(hessian) * a) :
                scale * cross(gradient, a)
        end
    else
        gb = field.grid_bins
        if !isempty(gb.indices)
            d = field.domain
            ix0 = floor(Int, (px - d[1]) / scale) - gb.origin[1] + 1
            iy0 = floor(Int, (py - d[3]) / scale) - gb.origin[2] + 1
            iz0 = floor(Int, (pz - d[5]) / scale) - gb.origin[3] + 1
            nx, ny, nz = gb.dims
            nx_ny = nx * ny
            @inbounds for iz in max(1, iz0 - 1):min(nz, iz0 + 1)
                zoff = (iz - 1) * nx_ny
                for iy in max(1, iy0 - 1):min(ny, iy0 + 1)
                    yzoff = zoff + (iy - 1) * nx
                    for ix in max(1, ix0 - 1):min(nx, ix0 + 1)
                        cidx = yzoff + ix
                        for ptr in gb.offsets[cidx]:(gb.offsets[cidx + 1] - 1)
                            j = gb.indices[ptr]
                            disp = pvec - SVector{3, T}(
                                T(centers[j, 1]), T(centers[j, 2]), T(centers[j, 3])
                            )
                            norm(disp) >= scale && continue
                            gradient, hessian = _kernel_derivatives_3d(disp, scale, kernel)
                            off = 3 + 3 * (j - 1)
                            a = SVector{3, T}(
                                T(coefficients[off + 1]), T(coefficients[off + 2]),
                                T(coefficients[off + 3])
                            )
                            acc += current ? scale * (hessian * a - tr(hessian) * a) :
                                scale * cross(gradient, a)
                        end
                    end
                end
            end
        end
    end
    return current ? acc / T(μ₀) : acc
end

_reference_value_3d(::NullField, point::SVector{3}; current = false) = zero(point)

@inline function _reference_value_3d(field, point; current = false)
    value = current ? _reference_current(field, point) : field(point)
    length(value) == 3 && all(isfinite, value) ||
        throw(ArgumentError("reference field must return three finite components"))
    return SVector{3}(value)
end

_reference_current(::NullField, pos) = zero(SVector{3, eltype(pos)})
_reference_current(::UniformField, pos) = zero(SVector{3, eltype(pos)})
_reference_current(field::AbstractMagneticField, pos) = current_density(field, pos)

function current_density(field::FieldReconstruction3D, pos)
    point = _check_3d_position(field, pos)
    return _reference_value_3d(field.reference_field, point; current = true) +
        _correction_3d(field, point, Val(true))
end

function _matrix_3d(rows, ncoefficients, ::Type{T}) where {T}
    nentries = sum(nnz, rows; init = 0)
    I, J, V = Vector{Int}(undef, nentries), Vector{Int}(undef, nentries),
        Vector{T}(undef, nentries)
    entry = 0
    for (i, row) in enumerate(rows)
        indices, values = findnz(row)
        for k in eachindex(indices)
            entry += 1
            I[entry], J[entry], V[entry] = i, indices[k], values[k]
        end
    end
    return sparse(I, J, V, length(rows), ncoefficients)
end

mutable struct _TripletRows3D{T}
    nrows::Int
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{T}
end

_TripletRows3D{T}() where {T} = _TripletRows3D{T}(0, Int[], Int[], T[])

@inline function _push_triplet!(rows::_TripletRows3D{T}, i::Int, j::Int, v::T) where {T}
    iszero(v) && return nothing
    push!(rows.I, i)
    push!(rows.J, j)
    push!(rows.V, v)
    return nothing
end

function _sparse_matrix(rows::_TripletRows3D{T}, ncoefficients::Int) where {T}
    return sparse(rows.I, rows.J, rows.V, rows.nrows, ncoefficients)
end

struct _KernelScratch3D{T}
    indices::Vector{Int}
    grads::Vector{SVector{3, T}}
    hessians::Vector{SMatrix{3, 3, T, 9}}
end

function _KernelScratch3D{T}() where {T}
    return _KernelScratch3D{T}(Int[], SVector{3, T}[], SMatrix{3, 3, T, 9}[])
end

function _collect_kernel_3d!(
        scratch::_KernelScratch3D{T}, point::SVector{3, T},
        centers, scale::T, kernel::Symbol, domain, gb::_GridBins3D
    ) where {T}
    empty!(scratch.indices)
    empty!(scratch.grads)
    empty!(scratch.hessians)
    if kernel === :gaussian
        @inbounds for j in axes(centers, 1)
            cj = SVector{3, T}(T(centers[j, 1]), T(centers[j, 2]), T(centers[j, 3]))
            disp = point - cj
            gradient, hessian = _kernel_derivatives_3d(disp, scale, kernel)
            push!(scratch.indices, j)
            push!(scratch.grads, gradient)
            push!(scratch.hessians, hessian)
        end
    elseif !isempty(gb.indices)
        ix0 = floor(Int, (point[1] - domain[1]) / scale) - gb.origin[1] + 1
        iy0 = floor(Int, (point[2] - domain[3]) / scale) - gb.origin[2] + 1
        iz0 = floor(Int, (point[3] - domain[5]) / scale) - gb.origin[3] + 1
        nx, ny, nz = gb.dims
        nx_ny = nx * ny
        @inbounds for iz in max(1, iz0 - 1):min(nz, iz0 + 1)
            zoff = (iz - 1) * nx_ny
            for iy in max(1, iy0 - 1):min(ny, iy0 + 1)
                yzoff = zoff + (iy - 1) * nx
                for ix in max(1, ix0 - 1):min(nx, ix0 + 1)
                    cidx = yzoff + ix
                    for ptr in gb.offsets[cidx]:(gb.offsets[cidx + 1] - 1)
                        j = gb.indices[ptr]
                        disp = point - SVector{3, T}(
                            T(centers[j, 1]), T(centers[j, 2]), T(centers[j, 3])
                        )
                        norm(disp) >= scale && continue
                        gradient, hessian = _kernel_derivatives_3d(disp, scale, kernel)
                        push!(scratch.indices, j)
                        push!(scratch.grads, gradient)
                        push!(scratch.hessians, hessian)
                    end
                end
            end
        end
    end
    return length(scratch.indices)
end

function _push_field_component_row!(
        rows::_TripletRows3D{T}, m::Int, scale::T, scratch::_KernelScratch3D{T},
        weight::T = one(T)
    ) where {T}
    row = (rows.nrows += 1)
    _push_triplet!(rows, row, m, weight)
    factor = weight * scale
    @inbounds for idx in eachindex(scratch.indices)
        off = 3 + 3 * (scratch.indices[idx] - 1)
        g = scratch.grads[idx]
        if m == 1
            _push_triplet!(rows, row, off + 2, -factor * g[3])
            _push_triplet!(rows, row, off + 3, factor * g[2])
        elseif m == 2
            _push_triplet!(rows, row, off + 1, factor * g[3])
            _push_triplet!(rows, row, off + 3, -factor * g[1])
        else
            _push_triplet!(rows, row, off + 1, -factor * g[2])
            _push_triplet!(rows, row, off + 2, factor * g[1])
        end
    end
    return row
end

function _push_directional_row!(
        rows::_TripletRows3D{T}, direction::SVector{3, T}, scale::T,
        scratch::_KernelScratch3D{T}, weight::T = one(T)
    ) where {T}
    row = (rows.nrows += 1)
    for k in 1:3
        _push_triplet!(rows, row, k, weight * direction[k])
    end
    factor = weight * scale
    @inbounds for idx in eachindex(scratch.indices)
        off = 3 + 3 * (scratch.indices[idx] - 1)
        w = factor * cross(direction, scratch.grads[idx])
        for k in 1:3
            _push_triplet!(rows, row, off + k, w[k])
        end
    end
    return row
end

function _push_curl_component_row!(
        rows::_TripletRows3D{T}, m::Int, scale::T, scratch::_KernelScratch3D{T},
        factor::T = one(T)
    ) where {T}
    row = (rows.nrows += 1)
    coeff = factor * scale
    @inbounds for idx in eachindex(scratch.indices)
        off = 3 + 3 * (scratch.indices[idx] - 1)
        h = scratch.hessians[idx]
        tr_h = tr(h)
        for k in 1:3
            val = coeff * (m == k ? h[m, m] - tr_h : h[m, k])
            _push_triplet!(rows, row, off + k, val)
        end
    end
    return row
end

function _add_observation_3d!(rows, targets, sigmas, groups, row, target, sigma, group)
    isnan(target) && return nothing
    isfinite(target) || throw(ArgumentError("observations cannot be infinite"))
    isfinite(sigma) && sigma > 0 ||
        throw(ArgumentError("uncertainties must be finite and positive"))
    push!(rows, row)
    push!(targets, target)
    push!(sigmas, sigma)
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
        centers, scale, kernel, domain, grid_bins::_GridBins3D,
        grid, mode, ::Type{T}
    ) where {T}
    ncoefficients = 3 + 3 * size(centers, 1)
    points, weights = _quadrature_3d(domain, grid, T)
    reg_rows = _TripletRows3D{T}()
    scratch = _KernelScratch3D{T}()
    unsupported = 0
    scale_T = T(scale)
    @inbounds for i in axes(points, 1)
        pt = SVector{3, T}(points[i, 1], points[i, 2], points[i, 3])
        count = _collect_kernel_3d!(
            scratch, pt, centers, scale_T, kernel, domain, grid_bins
        )
        unsupported += count == 0
        if mode === :field
            factor = sqrt(weights[i])
            for m in 1:3
                _push_field_component_row!(reg_rows, m, scale_T, scratch, factor)
            end
            curl_factor = factor * scale_T
            for m in 1:3
                _push_curl_component_row!(reg_rows, m, scale_T, scratch, curl_factor)
            end
        end
    end
    ridge_weight = mode === :coefficients ? one(T) : sqrt(eps(T))
    for k in 1:ncoefficients
        row = (reg_rows.nrows += 1)
        _push_triplet!(reg_rows, row, k, ridge_weight)
    end
    return _sparse_matrix(reg_rows, ncoefficients), unsupported, size(points, 1)
end

function _regularizer_3d(
        centers, scale, kernel, domain, bins::Dict{NTuple{3, Int}, Vector{Int}},
        grid, mode, ::Type{T}
    ) where {T}
    return _regularizer_3d(
        centers, scale, kernel, domain, _to_grid_bins_3d(bins), grid, mode, T
    )
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
    grid_bins = _to_grid_bins_3d(bins)
    scale_T = T(length_scale)
    ncoefficients = 3 + 3size(centers, 1)
    H_rows = _TripletRows3D{T}()
    exact_rows = _TripletRows3D{T}()
    targets, sigmas = T[], T[]
    groups = Symbol[]
    exact_targets = T[]
    scratch = _KernelScratch3D{T}()
    uncovered = 0
    inv_mu0 = inv(T(μ₀))
    @inbounds for i in axes(positions, 1)
        point = SVector{3, T}(positions[i, 1], positions[i, 2], positions[i, 3])
        count = _collect_kernel_3d!(
            scratch, point, centers, scale_T, kernel, box, grid_bins
        )
        uncovered += count == 0
        reference = _reference_value_3d(reference_field, point)
        for k in 1:3
            raw = B_observed[i, k]
            isnan(raw) && continue
            target = T(raw - reference[k])
            isfinite(target) || throw(ArgumentError("observations cannot be infinite"))
            sigma = T(uncertainties[i, k])
            isfinite(sigma) && sigma > 0 ||
                throw(ArgumentError("uncertainties must be finite and positive"))
            _push_field_component_row!(H_rows, k, scale_T, scratch)
            push!(targets, target)
            push!(sigmas, sigma)
            push!(groups, :probe)
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
            point = SVector{3, T}(locations[i, 1], locations[i, 2], locations[i, 3])
            direction = SVector{3, T}(directions[i, 1], directions[i, 2], directions[i, 3])
            count = _collect_kernel_3d!(
                scratch, point, centers, scale_T, kernel, box, grid_bins
            )
            uncovered += count == 0
            reference = _reference_value_3d(reference_field, point)
            raw = values[i]
            isnan(raw) && continue
            target = T(raw - dot(direction, reference))
            isfinite(target) || throw(ArgumentError("observations cannot be infinite"))
            if exact
                _push_directional_row!(exact_rows, direction, scale_T, scratch)
                push!(exact_targets, target)
            else
                sigma = errors === nothing ? one(T) : T(errors[i])
                isfinite(sigma) && sigma > 0 ||
                    throw(ArgumentError("uncertainties must be finite and positive"))
                _push_directional_row!(H_rows, direction, scale_T, scratch)
                push!(targets, target)
                push!(sigmas, sigma)
                push!(groups, group)
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
            point = SVector{3, T}(
                current_positions[i, 1], current_positions[i, 2], current_positions[i, 3]
            )
            count = _collect_kernel_3d!(
                scratch, point, centers, scale_T, kernel, box, grid_bins
            )
            uncovered += count == 0
            reference = _reference_value_3d(reference_field, point; current = true)
            for k in 1:3
                raw = current_density[i, k]
                isnan(raw) && continue
                target = T(raw - reference[k])
                isfinite(target) || throw(ArgumentError("observations cannot be infinite"))
                sigma = T(errors[i, k])
                isfinite(sigma) && sigma > 0 ||
                    throw(ArgumentError("uncertainties must be finite and positive"))
                _push_curl_component_row!(H_rows, k, scale_T, scratch, inv_mu0)
                push!(targets, target)
                push!(sigmas, sigma)
                push!(groups, :current)
            end
        end
    elseif current_uncertainties !== nothing
        throw(ArgumentError("current_uncertainties require current observations"))
    end
    H_rows.nrows == 0 && exact_rows.nrows == 0 &&
        throw(ArgumentError("a finite observation or exact constraint is required"))
    grid = regularization_grid === nothing ?
        ntuple(k -> max(3, ceil(Int, 2(box[2k] - box[2k - 1]) / length_scale) + 1), 3) :
        regularization_grid
    L, unsupported_grid, grid_count = _regularizer_3d(
        centers, scale_T, kernel, box, grid_bins, grid, regularization_mode, T
    )
    H = _sparse_matrix(H_rows, ncoefficients)
    weighted = Diagonal(inv.(sigmas)) * H
    A = vcat(weighted, sqrt(T(regularization)) * L)
    b = vcat(targets ./ sigmas, zeros(T, size(L, 1)))
    C = _sparse_matrix(exact_rows, ncoefficients)
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
        centers, coefficients, scale_T, box, kernel, bins, grid_bins,
        reference_field, diagnostics
    )
end
