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
        T = promote_type(typeof(I), typeof(r), eltype(center), eltype(direction))
        dir_hat = normalize(SVector{3, T}(direction))
        return new{T}(T(I), T(r), SVector{3, T}(center), dir_hat)
    end
end

@inline function (field::InfiniteWire)(pos)
    @boundscheck length(pos) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    T = eltype(pos)
    r_pos = SVector{3, T}(pos[1], pos[2], pos[3]) - field.center

    z_local = dot(r_pos, field.direction)
    rho_vec = r_pos - z_local * field.direction
    rho = norm(rho_vec)

    rho < 1.0e-15 && return @SVector zeros(T, 3)

    azimuthal_vec = cross(field.direction, rho_vec)

    B_vec = if rho <= field.r
        (μ₀ * field.I / (2π * field.r^2)) * azimuthal_vec
    else
        (μ₀ * field.I / (2π * rho^2)) * azimuthal_vec
    end

    return B_vec
end

@inline function (field::InfiniteWire)(x, y, z)
    return field(SVector(x, y, z))
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
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    @inbounds return getB_loop(r[1], r[2], r[3], loop)
end

@inline function getB_loop(x, y, z, loop::CurrentLoop)
    (; radius, current, center, normal) = loop

    n_hat = normal

    # Relative position from center, r_rel = r - center
    rx, ry, rz = x - center[1], y - center[2], z - center[3]
    r_rel = SVector(rx, ry, rz)

    # Project r_rel onto the normal vector to get z component in local coordinates
    z_local = dot(r_rel, n_hat)

    # Vector component perpendicular to n (rho vector)
    rho_vec = r_rel - z_local * n_hat
    rho = norm(rho_vec)

    # Handle the singularity on the wire itself (rho = a, z = 0)
    # and the center (rho = 0) separately if needed.
    if rho < 1.0e-10 * radius # On the axis
        # B is purely in n direction
        # B = μ0 I a^2 / (2 (a^2 + z^2)^(3/2))
        B_mag = μ₀ * current * radius^2 / (2 * (radius^2 + z_local^2)^1.5)
        return B_mag * n_hat
    end

    # Cylindrical components calculation
    # B_rho and B_z in local coordinates
    denom_sq = (radius + rho)^2 + z_local^2
    k_sq = 4 * radius * rho / denom_sq

    K_val = ellipk(k_sq)
    E_val = ellipe(k_sq)

    # Common factor
    factor = μ₀ * current / (2 * π * sqrt(denom_sq))

    # B_z (local)
    denom_diff_sq = (radius - rho)^2 + z_local^2
    B_z_local = factor * (K_val + (radius^2 - rho^2 - z_local^2) / denom_diff_sq * E_val)

    # B_rho (local)
    if abs(z_local) < 1.0e-15
        B_rho_local = 0.0
    else
        B_rho_local = factor * (z_local / rho) *
            (-K_val + (radius^2 + rho^2 + z_local^2) / denom_diff_sq * E_val)
    end

    # Transform back to global Cartesian coordinates
    rho_hat = rho_vec / rho

    B = B_rho_local * rho_hat + B_z_local * n_hat

    return B
end

# Make structs callable
@inline function (field::CurrentLoopAnalytic)(r)
    return getB_loop(r, field.loop)
end
