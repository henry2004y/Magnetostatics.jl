"""
    InfiniteWire{T} <: AbstractMagneticField

Magnetic field of an infinite straight wire.

# Fields
- `I::T`: Current in the wire [A].
- `r::T`: Radius of the wire [m].
- `center::SVector{3, T}`: A point on the wire [m].
- `direction::SVector{3, T}`: Unit vector pointing in the direction of the current.
"""
struct InfiniteWire{T} <: AbstractMagneticField
    I::T
    r::T
    center::SVector{3, T}
    direction::SVector{3, T}

    function InfiniteWire(
            I, r,
            center = SVector(0.0, 0.0, 0.0), direction = SVector(0.0, 0.0, 1.0)
        )
        T = float(promote_type(typeof(I), typeof(r), eltype(center), eltype(direction)))
        dir_hat = normalize(SVector{3, T}(direction))
        return new{T}(T(I), T(r), SVector{3, T}(center), dir_hat)
    end
end

@inline function (field::InfiniteWire{T})(pos) where {T}
    @boundscheck length(pos) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, pos)
    (; I, r, center, direction) = field
    @inbounds r_pos = SVector{3, FT}(pos[1], pos[2], pos[3]) - center

    z_local = dot(r_pos, direction)
    rho_vec = r_pos - z_local * direction
    rho_sq = dot(rho_vec, rho_vec)
    azimuthal_vec = cross(direction, rho_vec)

    if r > zero(T) && rho_sq <= r^2
        return (FT(μ₀) * I / (2 * FT(π) * r^2)) * azimuthal_vec
    end
    rho_sq < 1.0e-30 && return zero(SVector{3, FT})
    return (FT(μ₀) * I / (2 * FT(π) * rho_sq)) * azimuthal_vec
end

@inline function (field::InfiniteWire)(x, y, z)
    return field(SVector(x, y, z))
end

function current_density(field::InfiniteWire{T}, pos) where {T}
    @boundscheck length(pos) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, pos)
    (; I, r, center, direction) = field
    @inbounds r_pos = SVector{3, FT}(pos[1], pos[2], pos[3]) - center
    z_local = dot(r_pos, direction)
    rho_vec = r_pos - z_local * direction
    rho_sq = dot(rho_vec, rho_vec)
    if rho_sq <= r^2 && r > zero(T)
        return SVector{3, FT}((I / (FT(π) * r^2)) * direction)
    end
    return zero(SVector{3, FT})
end

"""
    CurrentLoopAnalytic{T} <: AbstractMagneticField

Analytical magnetic field of a circular current loop.
"""
struct CurrentLoopAnalytic{T} <: AbstractMagneticField
    loop::CurrentLoop{T}
end

"""
    getB_loop(r, loop::CurrentLoop)

Calculate the magnetic field `B` [T] at 3D point `r` from a `CurrentLoop`.
"""
@inline function getB_loop(r, loop::CurrentLoop)
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    @inbounds return getB_loop(r[1], r[2], r[3], loop)
end

@inline function getB_loop(x, y, z, loop::CurrentLoop{T}) where {T}
    FT = float(promote_type(T, typeof(x), typeof(y), typeof(z)))
    (; radius, current, center, normal) = loop
    n_hat = normal

    r_rel = SVector{3, FT}(x - center[1], y - center[2], z - center[3])
    z_local = dot(r_rel, n_hat)
    rho_vec = r_rel - z_local * n_hat
    rho_sq = dot(rho_vec, rho_vec)

    if rho_sq < (FT(1.0e-8) * radius)^2
        d2 = radius^2 + z_local^2
        d2 < FT(1.0e-30) && return zero(SVector{3, FT})
        b_z = FT(μ₀) * current * radius^2 / (2 * d2 * sqrt(d2))
        b_perp_fac = (3 * z_local / (2 * d2)) * b_z
        return b_perp_fac * rho_vec + b_z * n_hat
    end

    rho = sqrt(rho_sq)
    denom_sq = (radius + rho)^2 + z_local^2
    denom_diff_sq = (radius - rho)^2 + z_local^2
    denom_diff_sq < (FT(1.0e-12) * radius)^2 && return zero(SVector{3, FT})
    k_sq = 4 * radius * rho / denom_sq

    K_val, E_val = _ellipke(k_sq)

    factor = FT(μ₀) * current / (2 * FT(π) * sqrt(denom_sq))
    B_z_local = factor * (K_val + (radius^2 - rho_sq - z_local^2) / denom_diff_sq * E_val)
    B_rho_local = factor * (z_local / rho) *
        (-K_val + (radius^2 + rho_sq + z_local^2) / denom_diff_sq * E_val)

    rho_hat = rho_vec / rho
    return B_rho_local * rho_hat + B_z_local * n_hat
end

# Make structs callable
@inline function (field::CurrentLoopAnalytic)(r)
    return getB_loop(r, field.loop)
end

@inline function current_density(::CurrentLoopAnalytic{T}, pos) where {T}
    return zero(SVector{3, _promote_float(T, pos)})
end
