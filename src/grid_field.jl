"""
    GridMagneticField{T, A, Rx, Ry, Rz} <: AbstractMagneticField

Callable magnetic field backed by a `(3, Nx, Ny, Nz)` array on regular Cartesian grid
ranges `(xs, ys, zs)`. Supports zero-allocation `:linear` (trilinear) and `:cubic`
(tricubic Catmull-Rom) interpolation, `:throw` / `:clamp` / `:zero` out-of-bounds handling,
and finite-difference [`current_density`](@ref).
"""
struct GridMagneticField{
        T <: AbstractFloat,
        A <: AbstractArray{T, 4},
        Rx <: AbstractRange,
        Ry <: AbstractRange,
        Rz <: AbstractRange,
    } <: AbstractMagneticField
    B::A
    xs::Rx
    ys::Ry
    zs::Rz
    interpolation::Symbol
    extrapolate::Symbol
    dx::T
    dy::T
    dz::T
    inv_dx::T
    inv_dy::T
    inv_dz::T
end

function GridMagneticField(
        B::AbstractArray{T, 4},
        xs::AbstractRange,
        ys::AbstractRange,
        zs::AbstractRange;
        interpolation::Symbol = :linear,
        extrapolate::Symbol = :throw
    ) where {T}
    size(B, 1) == 3 || throw(ArgumentError("B must have size (3, Nx, Ny, Nz)"))
    nx, ny, nz = length(xs), length(ys), length(zs)
    size(B) == (3, nx, ny, nz) ||
        throw(ArgumentError("B size $(size(B)) does not match grid ($nx, $ny, $nz)"))
    nx >= 2 && ny >= 2 && nz >= 2 ||
        throw(ArgumentError("Grid must have at least 2 points along each axis"))
    interpolation in (:linear, :cubic) ||
        throw(ArgumentError("interpolation must be :linear or :cubic"))
    extrapolate in (:throw, :clamp, :zero) ||
        throw(ArgumentError("extrapolate must be :throw, :clamp, or :zero"))

    FT = T <: AbstractFloat ? T : Float64
    B_arr = B isa AbstractArray{FT, 4} ? B : FT.(B)
    dx = FT(step(xs))
    dy = FT(step(ys))
    dz = FT(step(zs))
    dx > 0 && dy > 0 && dz > 0 || throw(ArgumentError("Grid spacing must be positive"))
    return GridMagneticField(
        B_arr, xs, ys, zs, interpolation, extrapolate,
        dx, dy, dz, inv(dx), inv(dy), inv(dz)
    )
end

@inline function _grid_coords(
        field::GridMagneticField{T}, r
    ) where {T}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    x = T(r[1])
    y = T(r[2])
    z = T(r[3])
    xmin, xmax = T(first(field.xs)), T(last(field.xs))
    ymin, ymax = T(first(field.ys)), T(last(field.ys))
    zmin, zmax = T(first(field.zs)), T(last(field.zs))

    in_bounds = (xmin <= x <= xmax) && (ymin <= y <= ymax) && (zmin <= z <= zmax)
    if !in_bounds
        if field.extrapolate === :throw
            throw(
                DomainError(
                    SVector(x, y, z),
                    "Evaluation point lies outside GridMagneticField domain."
                )
            )
        elseif field.extrapolate === :zero
            return (false, 1, 1, 1, zero(T), zero(T), zero(T))
        end
        x = clamp(x, xmin, xmax)
        y = clamp(y, ymin, ymax)
        z = clamp(z, zmin, zmax)
    end

    nx, ny, nz = length(field.xs), length(field.ys), length(field.zs)
    fx = clamp((x - xmin) * field.inv_dx + one(T), one(T), T(nx))
    fy = clamp((y - ymin) * field.inv_dy + one(T), one(T), T(ny))
    fz = clamp((z - zmin) * field.inv_dz + one(T), one(T), T(nz))

    ix = min(floor(Int, fx), nx - 1)
    iy = min(floor(Int, fy), ny - 1)
    iz = min(floor(Int, fz), nz - 1)
    tx = fx - ix
    ty = fy - iy
    tz = fz - iz
    return (true, ix, iy, iz, tx, ty, tz)
end

@inline function _grid_b(B::AbstractArray{T, 4}, i::Int, j::Int, k::Int) where {T}
    @inbounds return SVector{3, T}(B[1, i, j, k], B[2, i, j, k], B[3, i, j, k])
end

@inline function _catmull_rom_weights(t::T) where {T}
    t2 = t * t
    t3 = t2 * t
    w0 = T(0.5) * (-t + 2 * t2 - t3)
    w1 = T(0.5) * (2 - 5 * t2 + 3 * t3)
    w2 = T(0.5) * (t + 4 * t2 - 3 * t3)
    w3 = T(0.5) * (-t2 + t3)
    return (w0, w1, w2, w3)
end

function (field::GridMagneticField{T})(r) where {T}
    ok, ix, iy, iz, tx, ty, tz = _grid_coords(field, r)
    ok || return zero(SVector{3, T})
    B = field.B
    if field.interpolation === :linear
        wx0, wx1 = one(T) - tx, tx
        wy0, wy1 = one(T) - ty, ty
        wz0, wz1 = one(T) - tz, tz
        c00 = wx0 * _grid_b(B, ix, iy, iz) + wx1 * _grid_b(B, ix + 1, iy, iz)
        c10 = wx0 * _grid_b(B, ix, iy + 1, iz) + wx1 * _grid_b(B, ix + 1, iy + 1, iz)
        c01 = wx0 * _grid_b(B, ix, iy, iz + 1) + wx1 * _grid_b(B, ix + 1, iy, iz + 1)
        c11 = wx0 * _grid_b(B, ix, iy + 1, iz + 1) +
            wx1 * _grid_b(B, ix + 1, iy + 1, iz + 1)
        c0 = wy0 * c00 + wy1 * c10
        c1 = wy0 * c01 + wy1 * c11
        return wz0 * c0 + wz1 * c1
    end

    nx, ny, nz = size(B, 2), size(B, 3), size(B, 4)
    wx = _catmull_rom_weights(tx)
    wy = _catmull_rom_weights(ty)
    wz = _catmull_rom_weights(tz)
    ixm = max(1, ix - 1)
    ixp = min(nx, ix + 2)
    iym = max(1, iy - 1)
    iyp = min(ny, iy + 2)
    izm = max(1, iz - 1)
    izp = min(nz, iz + 2)
    is = (ixm, ix, ix + 1, ixp)
    js = (iym, iy, iy + 1, iyp)
    ks = (izm, iz, iz + 1, izp)

    acc = zero(SVector{3, T})
    @inbounds for a in 1:4
        ka, wka = ks[a], wz[a]
        acc_j = zero(SVector{3, T})
        for b in 1:4
            jb, wjb = js[b], wy[b]
            acc_i = zero(SVector{3, T})
            for c in 1:4
                acc_i += wx[c] * _grid_b(B, is[c], jb, ka)
            end
            acc_j += wjb * acc_i
        end
        acc += wka * acc_j
    end
    return acc
end

@inline function _nodal_j(field::GridMagneticField{T}, i::Int, j::Int, k::Int) where {T}
    B = field.B
    nx, ny, nz = size(B, 2), size(B, 3), size(B, 4)
    @inbounds begin
        dBx_dy, dBz_dy = if j == 1
            (
                (B[1, i, 2, k] - B[1, i, 1, k]) * field.inv_dy,
                (B[3, i, 2, k] - B[3, i, 1, k]) * field.inv_dy,
            )
        elseif j == ny
            (
                (B[1, i, ny, k] - B[1, i, ny - 1, k]) * field.inv_dy,
                (B[3, i, ny, k] - B[3, i, ny - 1, k]) * field.inv_dy,
            )
        else
            s = T(0.5) * field.inv_dy
            (
                (B[1, i, j + 1, k] - B[1, i, j - 1, k]) * s,
                (B[3, i, j + 1, k] - B[3, i, j - 1, k]) * s,
            )
        end

        dBy_dx, dBz_dx = if i == 1
            (
                (B[2, 2, j, k] - B[2, 1, j, k]) * field.inv_dx,
                (B[3, 2, j, k] - B[3, 1, j, k]) * field.inv_dx,
            )
        elseif i == nx
            (
                (B[2, nx, j, k] - B[2, nx - 1, j, k]) * field.inv_dx,
                (B[3, nx, j, k] - B[3, nx - 1, j, k]) * field.inv_dx,
            )
        else
            s = T(0.5) * field.inv_dx
            (
                (B[2, i + 1, j, k] - B[2, i - 1, j, k]) * s,
                (B[3, i + 1, j, k] - B[3, i - 1, j, k]) * s,
            )
        end

        dBx_dz, dBy_dz = if k == 1
            (
                (B[1, i, j, 2] - B[1, i, j, 1]) * field.inv_dz,
                (B[2, i, j, 2] - B[2, i, j, 1]) * field.inv_dz,
            )
        elseif k == nz
            (
                (B[1, i, j, nz] - B[1, i, j, nz - 1]) * field.inv_dz,
                (B[2, i, j, nz] - B[2, i, j, nz - 1]) * field.inv_dz,
            )
        else
            s = T(0.5) * field.inv_dz
            (
                (B[1, i, j, k + 1] - B[1, i, j, k - 1]) * s,
                (B[2, i, j, k + 1] - B[2, i, j, k - 1]) * s,
            )
        end
    end
    inv_mu = inv(T(μ₀))
    return SVector{3, T}(
        (dBz_dy - dBy_dz) * inv_mu,
        (dBx_dz - dBz_dx) * inv_mu,
        (dBy_dx - dBx_dy) * inv_mu
    )
end

function current_density(field::GridMagneticField{T}, r) where {T}
    ok, ix, iy, iz, tx, ty, tz = _grid_coords(field, r)
    ok || return zero(SVector{3, T})
    wx0, wx1 = one(T) - tx, tx
    wy0, wy1 = one(T) - ty, ty
    wz0, wz1 = one(T) - tz, tz
    c00 = wx0 * _nodal_j(field, ix, iy, iz) + wx1 * _nodal_j(field, ix + 1, iy, iz)
    c10 = wx0 * _nodal_j(field, ix, iy + 1, iz) +
        wx1 * _nodal_j(field, ix + 1, iy + 1, iz)
    c01 = wx0 * _nodal_j(field, ix, iy, iz + 1) +
        wx1 * _nodal_j(field, ix + 1, iy, iz + 1)
    c11 = wx0 * _nodal_j(field, ix, iy + 1, iz + 1) +
        wx1 * _nodal_j(field, ix + 1, iy + 1, iz + 1)
    c0 = wy0 * c00 + wy1 * c10
    c1 = wy0 * c01 + wy1 * c11
    return wz0 * c0 + wz1 * c1
end
