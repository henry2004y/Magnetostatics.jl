"""
    HarrisSheet{T} <: AbstractMagneticField

Magnetic field of a Harris current sheet: B(z) = B0 * tanh(z/L) * x_hat.

# Fields
- `B0::T`: Asymptotic magnetic field strength.
- `L::T`: Half-width of the current sheet.
"""
struct HarrisSheet{T} <: AbstractMagneticField
    B0::T
    L::T
end

function (field::HarrisSheet)(r::SVector{3, T}) where {T}
    return SVector(field.B0 * tanh(r[3] / field.L), zero(T), zero(T))
end

"""
    Dipole{T} <: AbstractMagneticField

Magnetic field of a dipole moment M at the origin.

# Fields
- `M::SVector{3, T}`: Magnetic dipole moment.
"""
struct Dipole{T} <: AbstractMagneticField
    M::SVector{3, T}
end

function (field::Dipole)(r::SVector{3, T}) where {T}
    r_mag = norm(r)
    if r_mag < 1.0e-10
        return @SVector zeros(T, 3) # Singularity at origin
    end

    n = r / r_mag
    μ0_4π = 1.0e-7

    return (μ0_4π / r_mag^3) * (3 * dot(field.M, n) * n - field.M)
end

"""
    CurrentLoopAnalytic{T} <: AbstractMagneticField

Analytical magnetic field of a circular current loop.
"""
struct CurrentLoopAnalytic{T} <: AbstractMagneticField
    loop::CurrentLoop{T}
end

const μ₀ = 4π * 1.0e-7

"""
    getB_loop(r, loop::CurrentLoop)

Calculate the magnetic field `B` [T] at point `r` from a `CurrentLoop`.
"""
function getB_loop(r::AbstractVector, loop::CurrentLoop)
    (; radius, current, center, normal) = loop

    n_hat = normal

    # Relative position from center
    r_rel = r - center

    # Project r_rel onto the normal vector to get z component in local coordinates
    z_local = dot(r_rel, n_hat)

    # Vector component perpendicular to n (rho vector)
    rho_vec = r_rel - z_local * n_hat
    rho = norm(rho_vec)

    # Handle the singularity on the wire itself (rho = a, z = 0)
    # and the center (rho = 0) separately if needed.
    if rho < 1.0e-10 * radius # On the axis
        # B is purely in n direction
        # B = \mu_0 I a^2 / (2 (a^2 + z^2)^(3/2))
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

"""
    getB_loop(r, R, a, I, n)

Calculate the magnetic field `B` [T] at point `r` from a current loop with current `I` [A],
radius `a` [m], centered at `R`, and normal vector `n`.
"""
function getB_loop(
        r::AbstractVector, R::AbstractVector, a, I, n::AbstractVector
    )
    loop = CurrentLoop(a, I, R, n)
    return getB_loop(r, loop)
end

# Make structs callable
function (field::CurrentLoopAnalytic)(r)
    return getB_loop(r, field.loop)
end

function (solver::BiotSavart)(source::AbstractCurrentSource, r)
    return solve(solver, source, r)
end
