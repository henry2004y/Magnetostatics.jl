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
        T_float = float(T)
        return new{T_float}(
            T_float(B0), T_float(r0), T_float(v_sw), T_float(Omega),
            SVector{3, T_float}(center)
        )
    end
end

@inline function (field::ParkerSpiral{T})(pos) where {T}
    FT = _promote_float(T, pos)
    rx = FT(pos[1]) - field.center[1]
    ry = FT(pos[2]) - field.center[2]
    rz = FT(pos[3]) - field.center[3]
    r2 = rx^2 + ry^2 + rz^2
    r2 < 1.0e-20 && return zero(SVector{3, FT})

    r = sqrt(r2)
    inv_r = inv(r)
    (; B0, r0, v_sw, Omega) = field
    Br = B0 * (r0 * inv_r)^2
    bphi_over_rho = -Br * (Omega * (r - r0) / v_sw) * inv_r
    return SVector{3, FT}(
        Br * rx * inv_r - bphi_over_rho * ry,
        Br * ry * inv_r + bphi_over_rho * rx,
        Br * rz * inv_r
    )
end

function current_density(field::ParkerSpiral{T}, pos) where {T}
    FT = _promote_float(T, pos)
    rx = FT(pos[1]) - field.center[1]
    ry = FT(pos[2]) - field.center[2]
    rz = FT(pos[3]) - field.center[3]
    rho2 = rx^2 + ry^2
    r2 = rho2 + rz^2
    r2 < 1.0e-20 && return zero(SVector{3, FT})

    r = sqrt(r2)
    inv_r = inv(r)
    (; B0, r0, v_sw, Omega) = field
    cos_theta = rz * inv_r
    pref = B0 * r0^2 * Omega / (FT(μ₀) * v_sw * r2 * r)
    jr = -2 * pref * (r - r0) * cos_theta
    jt_over_rho = pref * r0 * inv_r
    j_perp_fac = jr * inv_r + jt_over_rho * cos_theta
    return SVector{3, FT}(
        j_perp_fac * rx,
        j_perp_fac * ry,
        jr * cos_theta - pref * r0 * rho2 / r2
    )
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
        T_float = float(T)
        return new{T_float}(T_float(B0), T_float(L), T_float(Bg), T_float(Bn))
    end
end

@inline function (field::GuideFieldHarrisSheet{T})(pos) where {T}
    FT = _promote_float(T, pos)
    z = FT(pos[3])
    (; B0, L, Bg, Bn) = field
    bx = B0 * tanh(z / L)
    return SVector{3, FT}(bx, Bg, Bn)
end

function current_density(field::GuideFieldHarrisSheet{T}, pos) where {T}
    FT = _promote_float(T, pos)
    z = FT(pos[3])
    (; B0, L) = field
    jy = (B0 / (FT(μ₀) * L)) * sech(z / L)^2
    return SVector{3, FT}(zero(FT), jy, zero(FT))
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
        T_float = float(T)
        return new{T_float}(T_float(B0), T_float(tau), SVector{2, T_float}(center))
    end
end

@inline function (field::GoldHoyleFluxTube{T})(pos) where {T}
    FT = _promote_float(T, pos)
    x = FT(pos[1]) - field.center[1]
    y = FT(pos[2]) - field.center[2]
    r2 = x * x + y * y
    (; B0, tau) = field

    denom = one(FT) + tau^2 * r2
    fac = (B0 * tau) / denom
    return SVector{3, FT}(-fac * y, fac * x, B0 / denom)
end

function current_density(field::GoldHoyleFluxTube{T}, pos) where {T}
    FT = _promote_float(T, pos)
    x = FT(pos[1]) - field.center[1]
    y = FT(pos[2]) - field.center[2]
    denom = one(FT) + field.tau^2 * (x * x + y * y)
    return (2 * field.tau / (FT(μ₀) * denom)) * field(pos)
end
