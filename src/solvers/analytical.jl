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

function (field::HarrisSheet)(r)
    T = eltype(r)
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

function (field::Dipole)(r)
    T = eltype(r)
    r_mag = norm(r)
    if r_mag < 1.0e-10
        return @SVector zeros(T, 3) # Singularity at origin
    end

    n = r / r_mag

    return (μ0_4π / r_mag^3) * (3 * dot(field.M, n) * n - field.M)
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
function getB_loop(r, loop::CurrentLoop)
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
function (field::CurrentLoopAnalytic)(r)
    return getB_loop(r, field.loop)
end

function (solver::BiotSavart)(source::AbstractCurrentSource, r)
    return solve(solver, source, r)
end

"""
    getB_mirror(x, y, z, distance, a, I1) :: SVector{3}
    getB_mirror(r, distance, a, I1) :: SVector{3}

Get magnetic field at `[x, y, z]` from a magnetic mirror generated from two coils.

# Arguments

  - `r`: location, vector of length 3 [m]
  - `x,y,z`: location [m]
  - `distance`: distance between solenoids in [m].
  - `a`: radius of each side coil in [m].
  - `I1`: current in the solenoid times number of windings in side coils.
"""
function getB_mirror(x, y, z, distance, a, I1)
    return getB_mirror(SVector(x, y, z), distance, a, I1)
end

function getB_mirror(r, distance, a, I1)
    cl1 = CurrentLoop(a, I1, SVector(0.0, 0.0, -0.5 * distance), SVector(0.0, 0.0, 1.0))
    cl2 = CurrentLoop(a, I1, SVector(0.0, 0.0, 0.5 * distance), SVector(0.0, 0.0, 1.0))
    B1 = getB_loop(r, cl1)
    B2 = getB_loop(r, cl2)

    return B1 + B2
end

"""
    getB_bottle(x, y, z, distance, a, b, I1, I2) :: SVector{3}
    getB_bottle(r, distance, a, b, I1, I2) :: SVector{3}

Get magnetic field from a magnetic bottle.
Reference: [wiki](https://en.wikipedia.org/wiki/Magnetic_mirror#Magnetic_bottles)

# Arguments

  - `r`: location, vector of length 3 [m]
  - `x,y,z`: location [m]
  - `distance::Float`: distance between solenoids in [m].
  - `a::Float`: radius of each side coil in [m].
  - `b::Float`: radius of central coil in [m].
  - `I1::Float`: current in the solenoid times number of windings in side coils in [A].
  - `I2::Float`: current in the central solenoid times number of windings in the
    central loop in [A].
"""
function getB_bottle(x, y, z, distance, a, b, I1, I2)
    return getB_bottle(SVector(x, y, z), distance, a, b, I1, I2)
end

function getB_bottle(r, distance, a, b, I1, I2)
    B = getB_mirror(r, distance, a, I1)

    # Central loop
    cl3 = CurrentLoop(b, I2, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))
    B3 = getB_loop(r, cl3)

    return B + B3
end

"""
    getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma) -> SVector{3}
    getB_tokamak_coil(r, a, b, ICoils, IPlasma) -> SVector{3}

Get the magnetic field from a Tokamak topology consists of 16 coils.
Original: [Tokamak-Fusion-Reactor](https://github.com/BoschSamuel/Simulation-of-a-Tokamak-Fusion-Reactor/blob/master/Simulation2.m)

# Arguments

  - `r`: location, vector of length 3 [m]
  - `x,y,z`: location [m]
  - `a`: radius of each coil in [m].
  - `b`: radius of central region in [m].
  - `ICoil`: current in the coil times number of windings in [A].
  - `IPlasma`: current of the plasma in [A].
"""
function getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma)
    return getB_tokamak_coil(SVector(x, y, z), a, b, ICoils, IPlasma)
end

function getB_tokamak_coil(r, a, b, ICoils, IPlasma)
    x, y, z = r
    a *= 2

    Bx, By, Bz = 0.0, 0.0, 0.0

    # magnetic field of the coils
    for i in 0:15
        θ = π / 16 + i * π / 8 # angle between the i-th coil and the x-axis

        # Coil center and normal
        # Center at (R_major * cos(θ), R_major * sin(θ), 0)
        # R_major = b + a
        R_major = b + a
        center = SVector(R_major * cos(θ), R_major * sin(θ), 0.0)

        # Normal is toroidal direction (perpendicular to poloidal plane)
        # Poloidal plane is at angle θ. Normal is (-sin(θ), cos(θ), 0)
        normal = SVector(-sin(θ), cos(θ), 0.0)

        cl = CurrentLoop(a, ICoils, center, normal)
        B_coil = getB_loop(r, cl)

        Bx += B_coil[1]
        By += B_coil[2]
        Bz += B_coil[3]
    end

    # magnetic field of the plasma current
    σ = a / 3 # parameter of the Gauss curve
    ϕ = atan(y, x)
    # distance to center of plasma ring
    # Plasma ring radius R_p = a + b
    R_p = a + b
    distance_val = √(z^2 + (x - R_p * cos(ϕ))^2 + (y - R_p * sin(ϕ))^2)

    if distance_val > 0.0001
        I2_r_plasma = IPlasma * erf(distance_val / (σ * √2))

        # Plasma is a horizontal loop at z=0, radius R_p
        cl_plasma = CurrentLoop(
            R_p, I2_r_plasma,
            SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0)
        )
        B_plasma = getB_loop(r, cl_plasma)

        Bx += B_plasma[1]
        By += B_plasma[2]
        Bz += B_plasma[3]
    end

    return SVector(Bx, By, Bz)
end

"""
    getB_tokamak_profile(x, y, z, q_profile, a, R₀, Bζ0) :: SVector{3}
    getB_tokamak_profile(r, q_profile, a, R₀, Bζ0) :: SVector{3}

Reconstruct the magnetic field distribution from a safe factor(q) profile.
Reference: Tokamak, 4th Edition, John Wesson.

# Arguments

  - `r`: location, vector of length 3 [m]
  - `x,y,z`: location [m]
  - `q_profile::Function`: profile of q. The variable of this function must be the normalized radius.
  - `a`: minor radius [m].
  - `R₀`: major radius [m].
  - `Bζ0`: toroidal magnetic field on axis [T].
"""
function getB_tokamak_profile(x, y, z, q_profile, a, R₀, Bζ0)
    return getB_tokamak_profile(SVector(x, y, z), q_profile, a, R₀, Bζ0)
end

function getB_tokamak_profile(pos, q_profile, a, R₀, Bζ0)
    x, y, z = pos
    R = √(x^2 + y^2)
    r = √((R - R₀)^2 + z^2)
    r > a && throw(OverflowError("out of vacuum vessel"))
    θ = atan(z, R - R₀)
    Bζ = Bζ0 * R₀ / R
    Bθ = r * Bζ / R₀ / q_profile(r / a)
    ζ = atan(y, x)

    Bx = -Bζ * sin(ζ) - Bθ * sin(θ) * cos(ζ)
    By = Bζ * cos(ζ) - Bθ * sin(θ) * sin(ζ)
    Bz = Bθ * cos(θ)

    return SVector(Bx, By, Bz)
end

"""
    getB_zpinch(x, y, z, I, a) -> SVector{3}
    getB_zpinch(r, I, a) -> SVector{3}

Get magnetic field from a Z-pinch configuration.
Reference: [Z-pinch](https://en.wikipedia.org/wiki/Z-pinch)

# Arguments

  - `r`: location, vector of length 3 [m]
  - `x,y,z`: location [m]
  - `I::Float`: current in the wire [A].
  - `a::Float`: radius of the wire [m].
"""
function getB_zpinch(x, y, z, I, a)
    return getB_zpinch(SVector(x, y, z), I, a)
end

function getB_zpinch(pos, I, a)
    x, y = pos[1], pos[2]
    r = hypot(x, y)
    if r < a
        factor = μ₀ * I / (2π * a^2)
        Bx = -factor * y
        By = factor * x
    else
        factor = μ₀ * I / (2π * r^2)
        Bx = -factor * y
        By = factor * x
    end

    return SVector(Bx, By, 0.0)
end
