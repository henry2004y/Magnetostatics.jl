"""
    ParkerSpiral{T} <: AbstractMagneticField

Parker spiral model of the heliospheric magnetic field (IMF) in the solar wind:
- Br(r, θ) = B0 * (r0 / r)^2
- Bϕ(r, θ) = -B0 * (r0 / r)^2 * (Ω * (r - r0) / v_sw) * sin(θ)
- Bθ = 0

# Fields
- `B0::T`: Reference radial magnetic field strength at radius `r0` [T].
- `r0::T`: Reference radius [m] (e.g. 1 AU = 1.496e11 m or 1 R_sun = 6.96e8 m).
- `v_sw::T`: Solar wind radial velocity [m/s] (default: 400 km/s).
- `Omega::T`: Solar angular rotation speed [rad/s] (default: 2.86e-6 rad/s).
- `center::SVector{3, T}`: Origin coordinates of the central star [m].
"""
struct ParkerSpiral{T} <: AbstractMagneticField
    B0::T
    r0::T
    v_sw::T
    Omega::T
    center::SVector{3, T}

    function ParkerSpiral(
            B0, r0;
            v_sw = 4.0e5,
            Omega = 2.865e-6,
            center = SA[0.0, 0.0, 0.0]
        )
        T = promote_type(typeof(B0), typeof(r0), typeof(v_sw), typeof(Omega), eltype(center))
        T_float = T <: AbstractFloat ? T : Float64
        return new{T_float}(
            T_float(B0), T_float(r0), T_float(v_sw), T_float(Omega),
            SVector{3, T_float}(center)
        )
    end
end

@inline function (field::ParkerSpiral{T})(pos) where {T}
    rx = pos[1] - field.center[1]
    ry = pos[2] - field.center[2]
    rz = pos[3] - field.center[3]

    r = sqrt(rx^2 + ry^2 + rz^2)
    rho = hypot(rx, ry)

    if r < 1.0e-10
        return SVector{3, T}(zero(T), zero(T), zero(T))
    end

    (; B0, r0, v_sw, Omega) = field
    r_ratio = r0 / r
    Br = B0 * r_ratio^2

    # Azimuthal magnetic field component
    sin_theta = rho / r
    Bphi = -Br * (Omega * (r - r0) / v_sw) * sin_theta

    inv_r = inv(r)
    r_hat = SVector{3, T}(rx * inv_r, ry * inv_r, rz * inv_r)

    if rho < 1.0e-15
        return Br * r_hat
    end

    inv_rho = inv(rho)
    phi_hat = SVector{3, T}(-ry * inv_rho, rx * inv_rho, zero(T))

    return Br * r_hat + Bphi * phi_hat
end

"""
    GuideFieldHarrisSheet{T} <: AbstractMagneticField

Harris current sheet equilibrium with an out-of-plane guide field `Bg` and optional
normal field `Bn`:
- Bx(z) = B0 * tanh(z / L)
- By(z) = Bg
- Bz(z) = Bn

# Fields
- `B0::T`: Reconnecting asymptotic magnetic field [T].
- `L::T`: Half-thickness of the current sheet [m].
- `Bg::T`: Uniform guide field [T].
- `Bn::T`: Normal magnetic field [T].
"""
struct GuideFieldHarrisSheet{T} <: AbstractMagneticField
    B0::T
    L::T
    Bg::T
    Bn::T

    function GuideFieldHarrisSheet(B0, L; Bg = 0.0, Bn = 0.0)
        T = promote_type(typeof(B0), typeof(L), typeof(Bg), typeof(Bn))
        T_float = T <: AbstractFloat ? T : Float64
        return new{T_float}(T_float(B0), T_float(L), T_float(Bg), T_float(Bn))
    end
end

@inline function (field::GuideFieldHarrisSheet{T})(pos) where {T}
    z = pos[3]
    (; B0, L, Bg, Bn) = field
    bx = B0 * tanh(z / L)
    return SVector{3, T}(bx, Bg, Bn)
end

function current_density(field::GuideFieldHarrisSheet{T}, pos) where {T}
    z = pos[3]
    (; B0, L) = field
    jy = (B0 / (T(μ₀) * L)) * sech(z / L)^2
    return SVector{3, T}(zero(T), jy, zero(T))
end

function vector_potential(field::GuideFieldHarrisSheet{T}, pos) where {T}
    x, z = pos[1], pos[3]
    (; B0, L, Bg, Bn) = field
    z_norm = abs(z / L)
    val = z_norm + log1p(exp(-2 * z_norm)) - log(T(2))
    ay = -B0 * L * val + Bn * x
    ax = -Bg * z
    return SVector{3, T}(ax, ay, zero(T))
end

"""
    GoldHoyleFluxTube{T} <: AbstractMagneticField

Gold-Hoyle uniformly twisted force-free cylindrical magnetic flux tube (coronal loop):
- Bz(r) = B0 / (1 + τ² r²)
- Bθ(r) = B0 * τ * r / (1 + τ² r²)
- Br = 0

# Fields
- `B0::T`: Axial magnetic field strength on axis [T].
- `tau::T`: Twist parameter τ = 2π / L_pitch [1/m].
- `center::SVector{2, T}`: Coordinates (x0, y0) of the flux tube axis [m].
"""
struct GoldHoyleFluxTube{T} <: AbstractMagneticField
    B0::T
    tau::T
    center::SVector{2, T}

    function GoldHoyleFluxTube(B0, tau; center = SA[0.0, 0.0])
        T = promote_type(typeof(B0), typeof(tau), eltype(center))
        T_float = T <: AbstractFloat ? T : Float64
        return new{T_float}(T_float(B0), T_float(tau), SVector{2, T_float}(center))
    end
end

@inline function (field::GoldHoyleFluxTube{T})(pos) where {T}
    x = pos[1] - field.center[1]
    y = pos[2] - field.center[2]
    r = hypot(x, y)
    (; B0, tau) = field

    denom = one(T) + (tau * r)^2
    bz = B0 / denom

    if r < 1.0e-15
        return SVector{3, T}(zero(T), zero(T), bz)
    end

    b_theta = (B0 * tau * r) / denom
    inv_r = inv(r)
    bx = -b_theta * (y * inv_r)
    by = b_theta * (x * inv_r)

    return SVector{3, T}(bx, by, bz)
end
