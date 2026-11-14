struct _GridBins3D
    origin::NTuple{3, Int}
    dims::NTuple{3, Int}
    offsets::Vector{Int}
    indices::Vector{Int}
end

"""
    FieldReconstruction3D <: AbstractMagneticField

Continuous Cartesian field represented as a uniform field `b₀`, the analytic curl of a
radial-basis vector potential, optional current-free solid harmonics (`vacuum_degree`),
virtual external 3D dipoles (`dipole_positions`), and parametric template fields
(`source_fields`). With `model = :div_free` (default) the fitted field is divergence-free
by construction; with `model = :curl_free` without `source_fields`, it is both
divergence-free and curl-free. Coordinates are in metres and field components in tesla.
Evaluation outside `domain` throws `DomainError`.
"""
struct FieldReconstruction3D{T, F <: AbstractMagneticField, D, SF} <: AbstractMagneticField
    centers::Matrix{T}
    coefficients::Vector{T}
    length_scale::T
    domain::NTuple{6, T}
    kernel::Symbol
    center_bins::Dict{NTuple{3, Int}, Vector{Int}}
    grid_bins::_GridBins3D
    reference_field::F
    diagnostics::D
    model::Symbol
    vacuum_coefficients::Vector{T}
    vacuum_degree::Int
    vacuum_center::SVector{3, T}
    vacuum_scale::T
    dipole_positions::Matrix{T}
    dipole_moments::Vector{SVector{3, T}}
    source_fields::SF
    source_weights::Vector{T}
end

function FieldReconstruction3D(
        centers::Matrix{T}, coefficients::Vector{T}, length_scale::T,
        domain::NTuple{6, T}, kernel::Symbol,
        center_bins::Dict{NTuple{3, Int}, Vector{Int}},
        reference_field::F, diagnostics::D
    ) where {T, F <: AbstractMagneticField, D}
    return FieldReconstruction3D{T, F, D, Tuple{}}(
        centers, coefficients, length_scale, domain, kernel,
        center_bins, _to_grid_bins_3d(center_bins), reference_field, diagnostics,
        :div_free, T[], -1, zero(SVector{3, T}), one(T),
        zeros(T, 0, 3), SVector{3, T}[], (), T[]
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

function _harmonic_count(vac_deg::Int)
    return vac_deg < 1 ? 0 : (vac_deg + 2)^2 - 4
end

@inline function _harmonic_norm(l::Int, m::Int, ::Type{T}) where {T <: AbstractFloat}
    p = one(T)
    for k in (l - m + 1):(l + m)
        p *= T(k)
    end
    odd_df = one(T)
    for k in 1:2:(2m - 1)
        odd_df *= T(k)
    end
    return odd_df / (T(l) * sqrt(p))
end

@inline function _foreach_harmonic_grad_3d(
        f::F, point::SVector{3, T}, center::SVector{3, T}, scale::T, vac_deg::Int
    ) where {F, T}
    vac_deg < 1 && return nothing
    inv_s = inv(scale)
    x = (point[1] - center[1]) * inv_s
    y = (point[2] - center[2]) * inv_s
    z = (point[3] - center[3]) * inv_s
    r2 = x * x + y * y + z * z
    L_pot = vac_deg + 1
    idx = 0
    c_mm = one(T)
    cx_mm = zero(T)
    cy_mm = zero(T)
    s_mm = zero(T)
    sx_mm = zero(T)
    sy_mm = zero(T)
    for m in 0:L_pot
        if m >= 1
            c_next = x * c_mm - y * s_mm
            cx_next = c_mm + x * cx_mm - y * sx_mm
            cy_next = x * cy_mm - s_mm - y * sy_mm
            s_next = y * c_mm + x * s_mm
            sx_next = y * cx_mm + s_mm + x * sx_mm
            sy_next = c_mm + y * cy_mm + x * sy_mm
            c_mm, cx_mm, cy_mm = c_next, cx_next, cy_next
            s_mm, sx_mm, sy_mm = s_next, sx_next, sy_next
        end
        if m >= 2
            nmm = _harmonic_norm(m, m, T)
            idx += 1
            f(idx, SVector{3, T}(cx_mm * nmm, cy_mm * nmm, zero(T)))
            idx += 1
            f(idx, SVector{3, T}(sx_mm * nmm, sy_mm * nmm, zero(T)))
        end
        c1, cx1, cy1, cz1 = c_mm, cx_mm, cy_mm, zero(T)
        c2, cx2, cy2, cz2 = zero(T), zero(T), zero(T), zero(T)
        s1, sx1, sy1, sz1 = s_mm, sx_mm, sy_mm, zero(T)
        s2, sx2, sy2, sz2 = zero(T), zero(T), zero(T), zero(T)
        for l in (m + 1):L_pot
            inv_lm = inv(T(l - m))
            a = T(2l - 1)
            b = T(l + m - 1)
            c0 = (a * z * c1 - b * r2 * c2) * inv_lm
            cx0 = (a * z * cx1 - b * (2 * x * c2 + r2 * cx2)) * inv_lm
            cy0 = (a * z * cy1 - b * (2 * y * c2 + r2 * cy2)) * inv_lm
            cz0 = (a * (c1 + z * cz1) - b * (2 * z * c2 + r2 * cz2)) * inv_lm
            nlm = _harmonic_norm(l, m, T)
            if l >= 2
                idx += 1
                f(idx, SVector{3, T}(cx0 * nlm, cy0 * nlm, cz0 * nlm))
            end
            c2, cx2, cy2, cz2 = c1, cx1, cy1, cz1
            c1, cx1, cy1, cz1 = c0, cx0, cy0, cz0
            if m >= 1
                s0 = (a * z * s1 - b * r2 * s2) * inv_lm
                sx0 = (a * z * sx1 - b * (2 * x * s2 + r2 * sx2)) * inv_lm
                sy0 = (a * z * sy1 - b * (2 * y * s2 + r2 * sy2)) * inv_lm
                sz0 = (a * (s1 + z * sz1) - b * (2 * z * s2 + r2 * sz2)) * inv_lm
                if l >= 2
                    idx += 1
                    f(idx, SVector{3, T}(sx0 * nlm, sy0 * nlm, sz0 * nlm))
                end
                s2, sx2, sy2, sz2 = s1, sx1, sy1, sz1
                s1, sx1, sy1, sz1 = s0, sx0, sy0, sz0
            end
        end
    end
    return nothing
end

@inline function _eval_harmonic_field_3d(
        vac_c::Vector{Tc}, point::SVector{3, T}, center::SVector{3, T},
        scale::T, vac_deg::Int
    ) where {Tc, T}
    vac_deg < 1 && return zero(SVector{3, T})
    inv_s = inv(scale)
    x = (point[1] - center[1]) * inv_s
    y = (point[2] - center[2]) * inv_s
    z = (point[3] - center[3]) * inv_s
    r2 = x * x + y * y + z * z
    L_pot = vac_deg + 1
    idx = 0
    bx = zero(T)
    by = zero(T)
    bz = zero(T)
    c_mm = one(T)
    cx_mm = zero(T)
    cy_mm = zero(T)
    s_mm = zero(T)
    sx_mm = zero(T)
    sy_mm = zero(T)
    @inbounds for m in 0:L_pot
        if m >= 1
            c_next = x * c_mm - y * s_mm
            cx_next = c_mm + x * cx_mm - y * sx_mm
            cy_next = x * cy_mm - s_mm - y * sy_mm
            s_next = y * c_mm + x * s_mm
            sx_next = y * cx_mm + s_mm + x * sx_mm
            sy_next = c_mm + y * cy_mm + x * sy_mm
            c_mm, cx_mm, cy_mm = c_next, cx_next, cy_next
            s_mm, sx_mm, sy_mm = s_next, sx_next, sy_next
        end
        if m >= 2
            nmm = _harmonic_norm(m, m, T)
            idx += 1
            w1 = T(vac_c[idx]) * nmm
            bx += w1 * cx_mm
            by += w1 * cy_mm
            idx += 1
            w2 = T(vac_c[idx]) * nmm
            bx += w2 * sx_mm
            by += w2 * sy_mm
        end
        c1, cx1, cy1, cz1 = c_mm, cx_mm, cy_mm, zero(T)
        c2, cx2, cy2, cz2 = zero(T), zero(T), zero(T), zero(T)
        s1, sx1, sy1, sz1 = s_mm, sx_mm, sy_mm, zero(T)
        s2, sx2, sy2, sz2 = zero(T), zero(T), zero(T), zero(T)
        for l in (m + 1):L_pot
            inv_lm = inv(T(l - m))
            a = T(2l - 1)
            b = T(l + m - 1)
            c0 = (a * z * c1 - b * r2 * c2) * inv_lm
            cx0 = (a * z * cx1 - b * (2 * x * c2 + r2 * cx2)) * inv_lm
            cy0 = (a * z * cy1 - b * (2 * y * c2 + r2 * cy2)) * inv_lm
            cz0 = (a * (c1 + z * cz1) - b * (2 * z * c2 + r2 * cz2)) * inv_lm
            nlm = _harmonic_norm(l, m, T)
            if l >= 2
                idx += 1
                w = T(vac_c[idx]) * nlm
                bx += w * cx0
                by += w * cy0
                bz += w * cz0
            end
            c2, cx2, cy2, cz2 = c1, cx1, cy1, cz1
            c1, cx1, cy1, cz1 = c0, cx0, cy0, cz0
            if m >= 1
                s0 = (a * z * s1 - b * r2 * s2) * inv_lm
                sx0 = (a * z * sx1 - b * (2 * x * s2 + r2 * sx2)) * inv_lm
                sy0 = (a * z * sy1 - b * (2 * y * s2 + r2 * sy2)) * inv_lm
                sz0 = (a * (s1 + z * sz1) - b * (2 * z * s2 + r2 * sz2)) * inv_lm
                if l >= 2
                    idx += 1
                    ws = T(vac_c[idx]) * nlm
                    bx += ws * sx0
                    by += ws * sy0
                    bz += ws * sz0
                end
                s2, sx2, sy2, sz2 = s1, sx1, sy1, sz1
                s1, sx1, sy1, sz1 = s0, sx0, sy0, sz0
            end
        end
    end
    return SVector{3, T}(bx, by, bz)
end

@inline function _dipole_scale(domain::NTuple{6, T}, rd::SVector{3, T}) where {T}
    dx = max(zero(T), domain[1] - rd[1], rd[1] - domain[2])
    dy = max(zero(T), domain[3] - rd[2], rd[2] - domain[4])
    dz = max(zero(T), domain[5] - rd[3], rd[3] - domain[6])
    return max(sqrt(dx * dx + dy * dy + dz * dz), eps(T))
end

@inline function _dipole_basis_matrix(
        point::SVector{3, T}, rd::SVector{3, T}, d_scale::T
    ) where {T}
    disp = point - rd
    inv_r = inv(norm(disp))
    n = disp * inv_r
    fac = (d_scale * inv_r)^3
    eye = one(SMatrix{3, 3, T, 9})
    return fac * (T(3) * (n * transpose(n)) - eye)
end

@inline function _eval_source_term(::Val{false}, sf, point::SVector{3, T}) where {T}
    return SVector{3, T}(sf(point))
end

@inline function _eval_source_term(::Val{true}, sf, point::SVector{3, T}) where {T}
    return SVector{3, T}(_reference_current(sf, point))
end

@inline function _eval_sources_3d(
        ::Tuple{}, weights, offset::Int, point::SVector{3, T}, ::Val{current}
    ) where {T, current}
    return zero(SVector{3, T})
end

@inline function _eval_sources_3d(
        sources::Tuple{Any, Vararg{Any}}, weights, offset::Int,
        point::SVector{3, T}, ::Val{current}
    ) where {T, current}
    term = weights[offset] * _eval_source_term(Val(current), sources[1], point)
    return term + _eval_sources_3d(
        Base.tail(sources), weights, offset + 1, point, Val(current)
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
    if !isempty(centers)
        if kernel === :gaussian
            @inbounds for j in axes(centers, 1)
                cj = SVector{3, T}(T(centers[j, 1]), T(centers[j, 2]), T(centers[j, 3]))
                disp = pvec - cj
                gradient, hessian = _kernel_derivatives_3d(disp, scale, kernel)
                off = 3 + 3 * (j - 1)
                a = SVector{3, T}(
                    T(coefficients[off + 1]), T(coefficients[off + 2]),
                    T(coefficients[off + 3])
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
                                gradient, hessian = _kernel_derivatives_3d(
                                    disp, scale, kernel
                                )
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
    end
    current && (acc /= T(μ₀))
    if !current
        vac_deg = field.vacuum_degree
        if vac_deg >= 1
            acc += _eval_harmonic_field_3d(
                field.vacuum_coefficients, pvec, SVector{3, T}(field.vacuum_center),
                T(field.vacuum_scale), vac_deg
            )
        end
        dips = field.dipole_positions
        moms = field.dipole_moments
        if !isempty(moms)
            @inbounds for j in eachindex(moms)
                rd = SVector{3, T}(T(dips[j, 1]), T(dips[j, 2]), T(dips[j, 3]))
                disp = pvec - rd
                rmag = norm(disp)
                n = disp / rmag
                m_j = SVector{3, T}(moms[j])
                acc += (T(μ0_4π) / rmag^3) * (T(3) * dot(m_j, n) * n - m_j)
            end
        end
    end
    if !isempty(field.source_fields)
        acc += _eval_sources_3d(
            field.source_fields, field.source_weights, 1, pvec, Val(current)
        )
    end
    return acc
end

_reference_value_3d(::NullField, point::SVector{3}; current = false) = zero(point)

@inline function _reference_value_3d(field, point; current = false)
    value = current ? _reference_current(field, point) : field(point)
    length(value) == 3 && all(isfinite, value) ||
        throw(ArgumentError("reference field must return three finite components"))
    return SVector{3}(value)
end

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

function _push_ext_field_terms!(
        rows::_TripletRows3D{T}, row::Int, point::SVector{3, T},
        m::Int, weight::T, ext
    ) where {T}
    if ext.vac_deg >= 1
        @inbounds _foreach_harmonic_grad_3d(
            point, ext.vac_center, ext.vac_scale, ext.vac_deg
        ) do idx, grad
            _push_triplet!(rows, row, ext.offset_vac + idx, weight * grad[m])
        end
    end
    if !isempty(ext.dipoles)
        @inbounds for j in axes(ext.dipoles, 1)
            rd = SVector{3, T}(ext.dipoles[j, 1], ext.dipoles[j, 2], ext.dipoles[j, 3])
            M_j = _dipole_basis_matrix(point, rd, ext.dip_scales[j])
            off = ext.offset_dip + 3 * (j - 1)
            for k in 1:3
                _push_triplet!(rows, row, off + k, weight * M_j[m, k])
            end
        end
    end
    for (idx, sf) in enumerate(ext.source_fields)
        b_sf = _eval_source_term(Val(false), sf, point)
        _push_triplet!(rows, row, ext.offset_src + idx, weight * b_sf[m])
    end
    return nothing
end

function _push_ext_dir_terms!(
        rows::_TripletRows3D{T}, row::Int, point::SVector{3, T},
        direction::SVector{3, T}, weight::T, ext
    ) where {T}
    if ext.vac_deg >= 1
        @inbounds _foreach_harmonic_grad_3d(
            point, ext.vac_center, ext.vac_scale, ext.vac_deg
        ) do idx, grad
            _push_triplet!(rows, row, ext.offset_vac + idx, weight * dot(direction, grad))
        end
    end
    if !isempty(ext.dipoles)
        @inbounds for j in axes(ext.dipoles, 1)
            rd = SVector{3, T}(ext.dipoles[j, 1], ext.dipoles[j, 2], ext.dipoles[j, 3])
            M_j = _dipole_basis_matrix(point, rd, ext.dip_scales[j])
            w = weight * (M_j * direction)
            off = ext.offset_dip + 3 * (j - 1)
            for k in 1:3
                _push_triplet!(rows, row, off + k, w[k])
            end
        end
    end
    for (idx, sf) in enumerate(ext.source_fields)
        b_sf = _eval_source_term(Val(false), sf, point)
        _push_triplet!(rows, row, ext.offset_src + idx, weight * dot(direction, b_sf))
    end
    return nothing
end

function _push_ext_curl_terms!(
        rows::_TripletRows3D{T}, row::Int, point::SVector{3, T},
        m::Int, factor::T, ext
    ) where {T}
    for (idx, sf) in enumerate(ext.source_fields)
        j_sf = _eval_source_term(Val(true), sf, point)
        _push_triplet!(rows, row, ext.offset_src + idx, factor * j_sf[m])
    end
    return nothing
end

function _regularizer_3d(
        centers, scale, kernel, domain, grid_bins::_GridBins3D,
        grid, mode, ::Type{T}, ext = nothing
    ) where {T}
    ncenters = size(centers, 1)
    ncoefficients = ext === nothing ? 3 + 3 * ncenters : ext.ncoefficients
    points, weights = _quadrature_3d(domain, grid, T)
    reg_rows = _TripletRows3D{T}()
    scratch = _KernelScratch3D{T}()
    unsupported = 0
    scale_T = T(scale)
    has_curl = ncenters > 0 || (ext !== nothing && !isempty(ext.source_fields))
    @inbounds for i in axes(points, 1)
        pt = SVector{3, T}(points[i, 1], points[i, 2], points[i, 3])
        if ncenters > 0
            count = _collect_kernel_3d!(
                scratch, pt, centers, scale_T, kernel, domain, grid_bins
            )
            unsupported += count == 0
        end
        if mode === :field
            factor = sqrt(weights[i])
            for m in 1:3
                row = _push_field_component_row!(reg_rows, m, scale_T, scratch, factor)
                ext === nothing || _push_ext_field_terms!(reg_rows, row, pt, m, factor, ext)
            end
            if has_curl
                curl_factor = factor * scale_T
                for m in 1:3
                    row = _push_curl_component_row!(
                        reg_rows, m, scale_T, scratch, curl_factor
                    )
                    ext === nothing || _push_ext_curl_terms!(
                        reg_rows, row, pt, m, curl_factor * T(μ₀), ext
                    )
                end
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
        grid, mode, ::Type{T}, ext = nothing
    ) where {T}
    return _regularizer_3d(
        centers, scale, kernel, domain, _to_grid_bins_3d(bins), grid, mode, T, ext
    )
end

function _group_rms_3d(residual, groups, group)
    selected = residual[groups .== group]
    return isempty(selected) ? nothing : norm(selected) / sqrt(length(selected))
end

function _solve_3d(
        A, b, C, d, n_obs::Int, store_covariance::Bool,
        algorithm, abstol, reltol, maxiters
    )
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
        if size(basis, 2) == 0
            cov = store_covariance ? zeros(T, size(A, 2), size(A, 2)) : nothing
            return particular, :constrained, 0, cov
        end
    end
    reduced = basis === nothing ? A : A * basis
    solution = LinearSolve.solve(
        LinearProblem(reduced, b - A * particular), algorithm; abstol, reltol, maxiters
    )
    LinearSolve.SciMLBase.successful_retcode(solution) ||
        error("3D reconstruction did not converge: $(solution.retcode)")
    coefficients = basis === nothing ? solution.u : particular + basis * solution.u
    all(isfinite, coefficients) || error("3D reconstruction produced non-finite values")
    cov = store_covariance ? _qr_covariance(A, n_obs, basis, T) : nothing
    return coefficients, solution.retcode, solution.iters, cov
end

"""
    reconstruct_field_3d(positions, B_observed; domain, length_scale = nothing, ...)

Fit a Cartesian magnetic field. `positions` and `B_observed` are `n × 3` matrices. `NaN`
field components denote missing measurements; infinite values are rejected.
`uncertainties` has the observation shape and is in tesla. `domain` is
`(xmin, xmax, ymin, ymax, zmin, zmax)`.

With `model = :div_free` (default), the fitted field is divergence-free by construction:
`reference_field + b₀ + curl(A) + ∇Φ_vac + ∑ B_dip + ∑ wₖ Sₖ`, where `A` uses vector
potential RBFs with support radius `length_scale` (`:wendland_c2` default or `:gaussian`).
With `model = :curl_free`, no vector potential RBFs are included (so `length_scale` may be
omitted) and `vacuum_degree` defaults to `3` unless `dipole_positions` or `source_fields`
are supplied, yielding a strictly current-free reconstruction (when `source_fields` is
empty).

Optional basis extensions:
- `vacuum_degree::Union{Nothing, Integer} = nothing`: add current-free real solid harmonic
  scalar potential modes of field polynomial degree `1:vacuum_degree` (with uniform `b₀`
  covering degree 0).
- `dipole_positions = nothing`: `m × 3` matrix of virtual external 3D dipole positions lying
  strictly outside `domain`; fitted dipole moments are stored in `field.dipole_moments`.
- `source_fields = ()`: tuple or vector of `AbstractMagneticField` templates whose scalar
  amplitudes are fitted linearly and stored in `field.source_weights`.
- `store_covariance::Bool = false`: compute and store the frequentist parameter covariance
  in `field.diagnostics.coefficient_covariance` for use with [`field_covariance`](@ref).

Directional probes, boundary normal flux (`n ⋅ B`), current density observations, and
regularization options behave as described in the 3D inverse guide.
"""
function reconstruct_field_3d(
        positions::AbstractMatrix{<:Real}, B_observed::AbstractMatrix{<:Real};
        domain, length_scale = nothing, basis_centers = nothing, kernel = :wendland_c2,
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
        current_uncertainties = nothing,
        model::Symbol = :div_free,
        vacuum_degree::Union{Nothing, Integer} = nothing,
        dipole_positions = nothing,
        source_fields = (),
        store_covariance::Bool = false
    )
    length(domain) == 6 ||
        throw(ArgumentError("domain must contain three pairs of bounds"))
    model in (:div_free, :curl_free) ||
        throw(ArgumentError("model must be :div_free or :curl_free"))
    if model === :div_free && length_scale === nothing
        throw(ArgumentError("length_scale is required when model = :div_free"))
    end
    if model === :curl_free && basis_centers !== nothing
        throw(ArgumentError("basis_centers cannot be supplied when model = :curl_free"))
    end
    source_fields isa Union{Tuple, AbstractVector} ||
        throw(ArgumentError("source_fields must be a Tuple or AbstractVector"))
    src_tuple = source_fields isa Tuple ? source_fields : Tuple(source_fields)
    all(sf -> sf isa AbstractMagneticField, src_tuple) ||
        throw(ArgumentError("all entries in source_fields must be AbstractMagneticField"))
    if model === :curl_free && current_positions !== nothing && isempty(src_tuple)
        throw(
            ArgumentError(
                "curl_free model without source_fields cannot fit current observations"
            )
        )
    end
    scale_type = length_scale === nothing ? Float64 : typeof(length_scale)
    T = promote_type(Float64, eltype(positions), eltype(B_observed), scale_type)
    box = NTuple{6, T}(Tuple(domain))
    all(isfinite, box) && all(box[2k] > box[2k - 1] for k in 1:3) ||
        throw(ArgumentError("domain bounds must be finite and strictly increasing"))
    scale_T = length_scale === nothing ?
        minimum(box[2k] - box[2k - 1] for k in 1:3) / 2 : T(length_scale)
    for (name, value) in (
            ("length_scale", scale_T), ("regularization", regularization),
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

    vac_deg = if model === :curl_free && vacuum_degree === nothing &&
            dipole_positions === nothing && isempty(src_tuple)
        3
    elseif vacuum_degree === nothing
        -1
    else
        Int(vacuum_degree)
    end
    vacuum_degree === nothing || vac_deg >= 1 ||
        throw(ArgumentError("vacuum_degree must be >= 1"))
    vac_center = SVector{3, T}(
        (box[1] + box[2]) / 2, (box[3] + box[4]) / 2, (box[5] + box[6]) / 2
    )
    vac_scale = hypot(
        (box[2] - box[1]) / 2, (box[4] - box[3]) / 2, (box[6] - box[5]) / 2
    )
    n_vac = _harmonic_count(vac_deg)

    dipoles = dipole_positions === nothing ? zeros(T, 0, 3) : Matrix{T}(dipole_positions)
    size(dipoles, 2) == 3 ||
        throw(ArgumentError("dipole_positions must have columns (x, y, z)"))
    all(isfinite, dipoles) || throw(ArgumentError("dipole_positions must be finite"))
    n_dip = size(dipoles, 1)
    dip_scales = Vector{T}(undef, n_dip)
    for j in 1:n_dip
        rd = SVector{3, T}(dipoles[j, 1], dipoles[j, 2], dipoles[j, 3])
        inside = (box[1] <= rd[1] <= box[2]) &&
            (box[3] <= rd[2] <= box[4]) &&
            (box[5] <= rd[3] <= box[6])
        !inside || throw(ArgumentError("dipole_positions must lie outside the domain"))
        dip_scales[j] = _dipole_scale(box, rd)
    end

    centers = if model === :curl_free
        zeros(T, 0, 3)
    elseif basis_centers === nothing
        _centers_3d(box, scale_T, T)
    else
        Matrix{T}(basis_centers)
    end
    if model === :div_free
        _validate_xyz(centers, box, "basis_centers")
        size(centers, 1) > 0 || throw(ArgumentError("at least one center is required"))
    end
    ncenters = size(centers, 1)
    bins = ncenters == 0 ? Dict{NTuple{3, Int}, Vector{Int}}() :
        _bins_3d(centers, box, scale_T)
    grid_bins = _to_grid_bins_3d(bins)

    n_rbf = 3 + 3 * ncenters
    offset_vac = n_rbf
    offset_dip = offset_vac + n_vac
    offset_src = offset_dip + 3 * n_dip
    n_src = length(src_tuple)
    ncoefficients = offset_src + n_src
    ext = (;
        vac_deg, vac_center, vac_scale, dipoles, dip_scales,
        source_fields = src_tuple, offset_vac, offset_dip, offset_src, ncoefficients,
    )

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
        if ncenters > 0
            count = _collect_kernel_3d!(
                scratch, point, centers, scale_T, kernel, box, grid_bins
            )
            uncovered += count == 0
        end
        reference = _reference_value_3d(reference_field, point)
        for k in 1:3
            raw = B_observed[i, k]
            isnan(raw) && continue
            target = T(raw - reference[k])
            isfinite(target) || throw(ArgumentError("observations cannot be infinite"))
            sigma = T(uncertainties[i, k])
            isfinite(sigma) && sigma > 0 ||
                throw(ArgumentError("uncertainties must be finite and positive"))
            row = _push_field_component_row!(H_rows, k, scale_T, scratch)
            _push_ext_field_terms!(H_rows, row, point, k, one(T), ext)
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
            if ncenters > 0
                count = _collect_kernel_3d!(
                    scratch, point, centers, scale_T, kernel, box, grid_bins
                )
                uncovered += count == 0
            end
            reference = _reference_value_3d(reference_field, point)
            raw = values[i]
            isnan(raw) && continue
            target = T(raw - dot(direction, reference))
            isfinite(target) || throw(ArgumentError("observations cannot be infinite"))
            if exact
                row = _push_directional_row!(exact_rows, direction, scale_T, scratch)
                _push_ext_dir_terms!(exact_rows, row, point, direction, one(T), ext)
                push!(exact_targets, target)
            else
                sigma = errors === nothing ? one(T) : T(errors[i])
                isfinite(sigma) && sigma > 0 ||
                    throw(ArgumentError("uncertainties must be finite and positive"))
                row = _push_directional_row!(H_rows, direction, scale_T, scratch)
                _push_ext_dir_terms!(H_rows, row, point, direction, one(T), ext)
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
            if ncenters > 0
                count = _collect_kernel_3d!(
                    scratch, point, centers, scale_T, kernel, box, grid_bins
                )
                uncovered += count == 0
            end
            reference = _reference_value_3d(reference_field, point; current = true)
            for k in 1:3
                raw = current_density[i, k]
                isnan(raw) && continue
                target = T(raw - reference[k])
                isfinite(target) || throw(ArgumentError("observations cannot be infinite"))
                sigma = T(errors[i, k])
                isfinite(sigma) && sigma > 0 ||
                    throw(ArgumentError("uncertainties must be finite and positive"))
                row = _push_curl_component_row!(H_rows, k, scale_T, scratch, inv_mu0)
                _push_ext_curl_terms!(H_rows, row, point, k, one(T), ext)
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
        ntuple(k -> max(3, ceil(Int, 2(box[2k] - box[2k - 1]) / scale_T) + 1), 3) :
        regularization_grid
    L, unsupported_grid, grid_count = _regularizer_3d(
        centers, scale_T, kernel, box, grid_bins, grid, regularization_mode, T, ext
    )
    H = _sparse_matrix(H_rows, ncoefficients)
    weighted = Diagonal(inv.(sigmas)) * H
    A = vcat(weighted, sqrt(T(regularization)) * L)
    b = vcat(targets ./ sigmas, zeros(T, size(L, 1)))
    C = _sparse_matrix(exact_rows, ncoefficients)
    all_coefficients, retcode, iterations, covariance = _solve_3d(
        A, b, C, exact_targets, H_rows.nrows, store_covariance,
        algorithm, abstol, reltol, maxiters
    )
    residual = H * all_coefficients - targets
    exact_error = isempty(exact_targets) ? zero(T) :
        norm(C * all_coefficients - exact_targets, Inf)
    tolerance = 100eps(T) * max(one(T), norm(exact_targets))
    exact_error <= tolerance || error("exact boundary constraints lost numerical accuracy")

    rbf_coefficients = all_coefficients[1:n_rbf]
    vacuum_coefficients = all_coefficients[(offset_vac + 1):offset_dip]
    dipole_moments = Vector{SVector{3, T}}(undef, n_dip)
    for j in 1:n_dip
        off = offset_dip + 3 * (j - 1)
        c_j = SVector{3, T}(
            all_coefficients[off + 1], all_coefficients[off + 2], all_coefficients[off + 3]
        )
        dipole_moments[j] = c_j * (dip_scales[j]^3 / T(μ0_4π))
    end
    source_weights = all_coefficients[(offset_src + 1):ncoefficients]

    eff_uncovered = ncenters == 0 ? 0 : uncovered
    eff_unsupported_grid = ncenters == 0 ? 0 : unsupported_grid
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
        iterations = iterations, unsupported_observations = eff_uncovered,
        unsupported_grid_points = eff_unsupported_grid, grid_point_count = grid_count,
        sampled_support_fraction = ncenters == 0 ? one(T) :
            1 - eff_unsupported_grid / grid_count,
        coefficient_covariance = covariance,
    )
    if eff_uncovered > 0 || eff_unsupported_grid > 0
        message = "RBF support has gaps; only uniform/reference fields act there"
        @warn message unsupported_observations = eff_uncovered grid_gaps =
            eff_unsupported_grid
    end
    F = typeof(reference_field)
    D = typeof(diagnostics)
    SF = typeof(src_tuple)
    return FieldReconstruction3D{T, F, D, SF}(
        centers, rbf_coefficients, scale_T, box, kernel, bins, grid_bins,
        reference_field, diagnostics, model, vacuum_coefficients, vac_deg,
        vac_center, vac_scale, dipoles, dipole_moments, src_tuple, source_weights
    )
end
