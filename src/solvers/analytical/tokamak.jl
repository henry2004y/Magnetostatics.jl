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
@inline function getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma)
    return getB_tokamak_coil(SVector(x, y, z), a, b, ICoils, IPlasma)
end

@inline function getB_tokamak_coil(r, a, b, ICoils, IPlasma)
    x, y, z = r
    # Original author used 2a for radius in their loop configuration
    R_coil = a * 2

    Bx, By, Bz = 0.0, 0.0, 0.0

    # magnetic field of the coils
    for i in 0:15
        θ = π / 16 + i * π / 8 # angle between the i-th coil and the x-axis

        # Coil center and normal
        R_major = b + R_coil
        center = SVector(R_major * cos(θ), R_major * sin(θ), 0.0)

        # Normal is toroidal direction (perpendicular to poloidal plane)
        # Poloidal plane is at angle θ. Normal is (-sin(θ), cos(θ), 0)
        normal = SVector(-sin(θ), cos(θ), 0.0)

        cl = CurrentLoop(R_coil, ICoils, center, normal)
        B_coil = getB_loop(r, cl)

        Bx += B_coil[1]
        By += B_coil[2]
        Bz += B_coil[3]
    end

    # magnetic field of the plasma current
    σ = R_coil / 3 # parameter of the Gauss curve
    ϕ = atan(y, x)
    # distance to center of plasma ring
    R_p = R_coil + b # Plasma ring radius
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
@inline function getB_tokamak_profile(x, y, z, q_profile, a, R₀, Bζ0)
    return getB_tokamak_profile(SVector(x, y, z), q_profile, a, R₀, Bζ0)
end

@inline function getB_tokamak_profile(pos, q_profile, a, R₀, Bζ0)
    x, y, z = pos
    R = √(x^2 + y^2)
    r = √((R - R₀)^2 + z^2)
    r > a && throw(OverflowError("out of vacuum vessel"))
    θ = atan(z, R - R₀)
    Bζ = Bζ0 * R₀ / R
    Bθ = r * Bζ / R₀ / q_profile(r / a)
    ζ = atan(y, x)
    sinζ, cosζ = sincos(ζ)
    sinθ, cosθ = sincos(θ)

    Bx = -Bζ * sinζ - Bθ * sinθ * cosζ
    By = Bζ * cosζ - Bθ * sinθ * sinζ
    Bz = Bθ * cosθ

    return SVector(Bx, By, Bz)
end
