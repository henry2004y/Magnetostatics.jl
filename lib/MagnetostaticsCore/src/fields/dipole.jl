"""
    Dipole{T} <: AbstractMagneticField

Magnetic field of a point magnetic dipole moment `M` located at `center`
with optional core regularization `r_min` and scaling prefactor `k`
(defaulting to `μ₀ / 4π`; pass `k = 1.0` for normalized code units).

# Fields
- `M::SVector{3, T}`: Magnetic dipole moment [A·m²].
- `center::SVector{3, T}`: Position of the dipole [m].
- `r_min::T`: Minimum distance clamp radius (0 disables softening).
- `k::T`: Prefactor scaling `B = k * (3(M⋅r̂)r̂ - M) / r³`.
"""
struct Dipole{T} <: AbstractMagneticField
    M::SVector{3, T}
    center::SVector{3, T}
    r_min::T
    k::T

    function Dipole{T}(
            M::SVector{3, T}, center::SVector{3, T}, r_min::T, k::T
        ) where {T}
        r_min >= zero(T) || throw(ArgumentError("r_min must be non-negative"))
        return new{T}(M, center, r_min, k)
    end
end

function Dipole(
        M::AbstractVector;
        center::AbstractVector = SVector(0, 0, 0),
        r_min::Real = 0,
        k::Union{Real, Nothing} = nothing
    )
    T = k === nothing ?
        float(promote_type(eltype(M), eltype(center), typeof(r_min))) :
        float(promote_type(eltype(M), eltype(center), typeof(r_min), typeof(k)))
    k_val = k === nothing ? T(μ0_4π) : T(k)
    return Dipole{T}(
        SVector{3, T}(M[1], M[2], M[3]),
        SVector{3, T}(center[1], center[2], center[3]),
        T(r_min),
        k_val
    )
end

function Dipole(
        M::AbstractVector,
        center::AbstractVector;
        r_min::Real = 0,
        k::Union{Real, Nothing} = nothing
    )
    T = k === nothing ?
        float(promote_type(eltype(M), eltype(center), typeof(r_min))) :
        float(promote_type(eltype(M), eltype(center), typeof(r_min), typeof(k)))
    k_val = k === nothing ? T(μ0_4π) : T(k)
    return Dipole{T}(
        SVector{3, T}(M[1], M[2], M[3]),
        SVector{3, T}(center[1], center[2], center[3]),
        T(r_min),
        k_val
    )
end

function Dipole(;
        M::AbstractVector,
        center::AbstractVector = SVector(0, 0, 0),
        r_min::Real = 0,
        k::Union{Real, Nothing} = nothing
    )
    return Dipole(M, center; r_min = r_min, k = k)
end

@inline function (field::Dipole{T})(x, y, z) where {T}
    FT = float(promote_type(T, typeof(x), typeof(y), typeof(z)))
    r = SVector{3, FT}(
        x - field.center[1], y - field.center[2], z - field.center[3]
    )
    r_sq = dot(r, r)
    if field.r_min > zero(T)
        r_sq_eff = max(r_sq, FT(field.r_min)^2)
        inv_r = inv(sqrt(r_sq_eff))
        inv_r3 = inv_r / r_sq_eff
        inv_r5 = inv_r3 / r_sq_eff
        return FT(field.k) * (3 * dot(field.M, r) * r * inv_r5 - field.M * inv_r3)
    end
    r_sq < 1.0e-20 && return zero(SVector{3, FT})
    r_mag = sqrt(r_sq)
    n = r / r_mag
    return (FT(field.k) / (r_sq * r_mag)) * (3 * dot(field.M, n) * n - field.M)
end

@inline function (field::Dipole)(r)
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    @inbounds return field(r[1], r[2], r[3])
end

@inline function current_density(::Dipole{T}, r) where {T}
    return zero(SVector{3, _promote_float(T, r)})
end

@inline function vector_potential(field::Dipole{T}, r) where {T}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    @inbounds r_vec = SVector{3, FT}(
        r[1] - field.center[1], r[2] - field.center[2], r[3] - field.center[3]
    )
    r_sq = dot(r_vec, r_vec)
    if field.r_min > zero(T)
        r_sq_eff = max(r_sq, FT(field.r_min)^2)
        r_mag = sqrt(r_sq_eff)
        return (FT(field.k) / (r_sq_eff * r_mag)) * cross(field.M, r_vec)
    end
    r_sq < 1.0e-20 && return zero(SVector{3, FT})
    r_mag = sqrt(r_sq)
    return (FT(field.k) / (r_sq * r_mag)) * cross(field.M, r_vec)
end

@inline function field_and_potential(field::Dipole{T}, r) where {T}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    @inbounds r_vec = SVector{3, FT}(
        r[1] - field.center[1], r[2] - field.center[2], r[3] - field.center[3]
    )
    r_sq = dot(r_vec, r_vec)
    if field.r_min > zero(T)
        r_sq_eff = max(r_sq, FT(field.r_min)^2)
        inv_r = inv(sqrt(r_sq_eff))
        inv_r3 = inv_r / r_sq_eff
        inv_r5 = inv_r3 / r_sq_eff
        B = FT(field.k) * (3 * dot(field.M, r_vec) * r_vec * inv_r5 - field.M * inv_r3)
        A = (FT(field.k) * inv_r3) * cross(field.M, r_vec)
        return (B, A)
    end
    if r_sq < 1.0e-20
        z = zero(SVector{3, FT})
        return (z, z)
    end
    inv_r = inv(sqrt(r_sq))
    inv_r3 = inv_r / r_sq
    n = r_vec * inv_r
    B = (FT(field.k) * inv_r3) * (3 * dot(field.M, n) * n - field.M)
    A = (FT(field.k) * inv_r3) * cross(field.M, r_vec)
    return (B, A)
end

"""
    LineDipole2D{T} <: AbstractMagneticField

2D in-plane line dipole magnetic field with exact zero 2D divergence:
`B = k * (2 (M ⋅ r) r - r² M) / r⁴` in the `(x, y)` plane and `B_z = 0`.
Vector potential is `A = (0, 0, k * (M_x r_y - M_y r_x) / r²)`, satisfying `∇ × A = B`.

# Fields
- `M::SVector{2, T}`: 2D line dipole moment vector.
- `center::SVector{2, T}`: In-plane center `(x₀, y₀)`.
- `r_min::T`: Minimum distance clamp radius (0 disables softening).
- `k::T`: Prefactor (defaults to `μ₀ / 2π`; pass `k = 1.0` for normalized code units).
"""
struct LineDipole2D{T} <: AbstractMagneticField
    M::SVector{2, T}
    center::SVector{2, T}
    r_min::T
    k::T

    function LineDipole2D{T}(
            M::SVector{2, T}, center::SVector{2, T}, r_min::T, k::T
        ) where {T}
        r_min >= zero(T) || throw(ArgumentError("r_min must be non-negative"))
        return new{T}(M, center, r_min, k)
    end
end

function LineDipole2D(
        M::AbstractVector,
        center::AbstractVector = SVector(0, 0);
        r_min::Real = 0,
        k::Union{Real, Nothing} = nothing
    )
    T = k === nothing ?
        float(promote_type(eltype(M), eltype(center), typeof(r_min))) :
        float(promote_type(eltype(M), eltype(center), typeof(r_min), typeof(k)))
    k_val = k === nothing ? T(2 * μ0_4π) : T(k)
    return LineDipole2D{T}(
        SVector{2, T}(M[1], M[2]),
        SVector{2, T}(center[1], center[2]),
        T(r_min),
        k_val
    )
end

function LineDipole2D(;
        M::AbstractVector,
        center::AbstractVector = SVector(0, 0),
        r_min::Real = 0,
        k::Union{Real, Nothing} = nothing
    )
    return LineDipole2D(M, center; r_min = r_min, k = k)
end

@inline function (field::LineDipole2D{T})(x, y, z = zero(T)) where {T}
    FT = float(promote_type(T, typeof(x), typeof(y), typeof(z)))
    rx = x - field.center[1]
    ry = y - field.center[2]
    r2 = rx * rx + ry * ry
    r2_eff = if field.r_min > zero(T)
        max(r2, FT(field.r_min)^2)
    else
        r2 < 1.0e-20 && return zero(SVector{3, FT})
        r2
    end
    mx, my = field.M[1], field.M[2]
    mdotr = mx * rx + my * ry
    inv_r4 = inv(r2_eff * r2_eff)
    bx = FT(field.k) * (2 * mdotr * rx - r2_eff * mx) * inv_r4
    by = FT(field.k) * (2 * mdotr * ry - r2_eff * my) * inv_r4
    return SVector{3, FT}(bx, by, zero(FT))
end

@inline function (field::LineDipole2D)(r)
    @boundscheck length(r) >= 2 ||
        throw(ArgumentError("Input must have at least 2 elements."))
    @inbounds return field(r[1], r[2])
end

@inline function current_density(::LineDipole2D{T}, r) where {T}
    return zero(SVector{3, _promote_float(T, r)})
end

@inline function vector_potential(field::LineDipole2D{T}, r) where {T}
    @boundscheck length(r) >= 2 ||
        throw(ArgumentError("Input must have at least 2 elements."))
    FT = _promote_float(T, r)
    @inbounds rx = FT(r[1]) - field.center[1]
    @inbounds ry = FT(r[2]) - field.center[2]
    r2 = rx * rx + ry * ry
    r2_eff = if field.r_min > zero(T)
        max(r2, FT(field.r_min)^2)
    else
        r2 < 1.0e-20 && return zero(SVector{3, FT})
        r2
    end
    az = FT(field.k) * (field.M[1] * ry - field.M[2] * rx) / r2_eff
    return SVector{3, FT}(zero(FT), zero(FT), az)
end

@inline function field_and_potential(field::LineDipole2D{T}, r) where {T}
    @boundscheck length(r) >= 2 ||
        throw(ArgumentError("Input must have at least 2 elements."))
    FT = _promote_float(T, r)
    @inbounds rx = FT(r[1]) - field.center[1]
    @inbounds ry = FT(r[2]) - field.center[2]
    r2 = rx * rx + ry * ry
    r2_eff = if field.r_min > zero(T)
        max(r2, FT(field.r_min)^2)
    else
        if r2 < 1.0e-20
            z = zero(SVector{3, FT})
            return (z, z)
        end
        r2
    end
    mx, my = field.M[1], field.M[2]
    mdotr = mx * rx + my * ry
    inv_r2 = inv(r2_eff)
    inv_r4 = inv_r2 * inv_r2
    kf = FT(field.k)
    B = SVector{3, FT}(
        kf * (2 * mdotr * rx - r2_eff * mx) * inv_r4,
        kf * (2 * mdotr * ry - r2_eff * my) * inv_r4,
        zero(FT)
    )
    A = SVector{3, FT}(zero(FT), zero(FT), kf * (mx * ry - my * rx) * inv_r2)
    return (B, A)
end
