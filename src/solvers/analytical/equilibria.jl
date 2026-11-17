"""
    SolovevEquilibrium{T} <: AbstractMagneticField

Analytical Solov'ev equilibrium for an axisymmetric tokamak plasma solving the
Grad-Shafranov equation Δ*ψ = -μ₀ R² p'(ψ) - F F'(ψ).

# Fields
- `B0::T`: Toroidal magnetic field on axis at `R0` [T].
- `R0::T`: Major radius of the magnetic axis [m].
- `kappa::T`: Plasma elongation κ.
- `q0::T`: On-axis safety factor.
- `psi0::T`: Characteristic poloidal flux scale `psi0 = (R0 * B0 * kappa) / (2 * q0)`.
"""
struct SolovevEquilibrium{T} <: AbstractMagneticField
    B0::T
    R0::T
    kappa::T
    q0::T
    psi0::T

    function SolovevEquilibrium(B0, R0, kappa = 1.0, q0 = 1.0)
        T = promote_type(typeof(B0), typeof(R0), typeof(kappa), typeof(q0))
        T_float = float(T)
        b0 = T_float(B0)
        r0 = T_float(R0)
        k = T_float(kappa)
        q = T_float(q0)
        psi0 = (r0 * b0 * k) / (2 * q)
        return new{T_float}(b0, r0, k, q, psi0)
    end
end

"""
    poloidal_flux(eq::SolovevEquilibrium, pos) -> Real

Compute the poloidal flux function ψ(R, Z) at Cartesian location `pos = [x, y, z]`.
"""
function poloidal_flux(eq::SolovevEquilibrium{T}, pos) where {T}
    FT = _promote_float(T, pos)
    x, y, z = FT(pos[1]), FT(pos[2]), FT(pos[3])
    R2 = x * x + y * y
    (; R0, kappa, psi0) = eq
    inv_R04 = inv(R0^4)
    term_z = (R2 * z^2) / (kappa^2)
    term_r = FT(0.25) * (R2 - R0^2)^2
    return psi0 * inv_R04 * (term_z + term_r)
end

@inline function (eq::SolovevEquilibrium{T})(pos) where {T}
    FT = _promote_float(T, pos)
    x, y, z = FT(pos[1]), FT(pos[2]), FT(pos[3])
    R2 = x * x + y * y
    R2 < 1.0e-20 && return zero(SVector{3, FT})

    (; B0, R0, kappa, psi0) = eq
    inv_R04 = inv(R0^4)

    br_over_R = -2 * psi0 * inv_R04 * z / kappa^2
    bphi_over_R = (R0 * B0) / R2
    Bx = br_over_R * x - bphi_over_R * y
    By = br_over_R * y + bphi_over_R * x
    Bz = psi0 * inv_R04 * ((2 * z^2) / kappa^2 + (R2 - R0^2))

    return SVector{3, FT}(Bx, By, Bz)
end

function current_density(eq::SolovevEquilibrium{T}, pos) where {T}
    FT = _promote_float(T, pos)
    x, y = FT(pos[1]), FT(pos[2])
    x * x + y * y < 1.0e-20 && return zero(SVector{3, FT})
    (; R0, kappa, psi0) = eq
    j_over_R = -(2 * psi0 * (one(FT) + inv(kappa^2))) / (FT(μ₀) * R0^4)
    return SVector{3, FT}(-j_over_R * y, j_over_R * x, zero(FT))
end

"""
    TaylorState{T} <: AbstractMagneticField

Cylindrical relaxed force-free magnetic equilibrium satisfying ∇ × B = α B
(Woltjer-Taylor state / Lundquist field) in terms of Bessel functions J₀ and J₁.

# Fields
- `B0::T`: Magnetic field strength on axis [T].
- `alpha::T`: Force-free parameter α [1/m].
- `center::SVector{2, T}`: Transverse center coordinates (x0, y0) [m].
"""
struct TaylorState{T} <: AbstractMagneticField
    B0::T
    alpha::T
    center::SVector{2, T}

    function TaylorState(B0, alpha; center = SA[0.0, 0.0])
        T = promote_type(typeof(B0), typeof(alpha), eltype(center))
        T_float = float(T)
        return new{T_float}(T_float(B0), T_float(alpha), SVector{2, T_float}(center))
    end
end

@inline function (field::TaylorState{T})(pos) where {T}
    FT = _promote_float(T, pos)
    x = FT(pos[1]) - field.center[1]
    y = FT(pos[2]) - field.center[2]
    r2 = x * x + y * y
    (; B0, alpha) = field

    ar2 = alpha^2 * r2
    if ar2 < 1.0e-16
        b_over_r = B0 * (alpha / 2) * (one(FT) - ar2 / 8)
        bz = B0 * (one(FT) - ar2 / 4)
        return SVector{3, FT}(-b_over_r * y, b_over_r * x, bz)
    end

    r = sqrt(r2)
    ar = alpha * r
    b_over_r = (B0 * besselj1(ar)) / r
    bz = B0 * besselj0(ar)
    return SVector{3, FT}(-b_over_r * y, b_over_r * x, bz)
end

function current_density(field::TaylorState{T}, pos) where {T}
    FT = _promote_float(T, pos)
    return (field.alpha / FT(μ₀)) * field(pos)
end

"""
    ScrewPinch{T} <: AbstractMagneticField

Cylindrical screw pinch combining a uniform axial field `Bz0` with the azimuthal field
from an axial current `I` inside radius `a`.

# Fields
- `I::T`: Axial current [A].
- `a::T`: Pinch radius [m].
- `Bz0::T`: Axial magnetic field [T].
- `center::SVector{2, T}`: Transverse center coordinates [m].
"""
struct ScrewPinch{T} <: AbstractMagneticField
    I::T
    a::T
    Bz0::T
    center::SVector{2, T}

    function ScrewPinch(I, a, Bz0; center = SA[0.0, 0.0])
        T = promote_type(typeof(I), typeof(a), typeof(Bz0), eltype(center))
        T_float = float(T)
        return new{T_float}(
            T_float(I), T_float(a), T_float(Bz0), SVector{2, T_float}(center)
        )
    end
end

@inline function (field::ScrewPinch{T})(pos) where {T}
    FT = _promote_float(T, pos)
    x = FT(pos[1]) - field.center[1]
    y = FT(pos[2]) - field.center[2]
    r2 = x * x + y * y
    (; I, a, Bz0) = field

    if a > zero(T) && r2 <= a^2
        factor = (FT(μ₀) * I) / (2 * FT(π) * a^2)
        return SVector{3, FT}(-factor * y, factor * x, FT(Bz0))
    end
    r2 < 1.0e-30 && return SVector{3, FT}(zero(FT), zero(FT), FT(Bz0))
    factor = (FT(μ₀) * I) / (2 * FT(π) * r2)
    return SVector{3, FT}(-factor * y, factor * x, FT(Bz0))
end

function current_density(field::ScrewPinch{T}, pos) where {T}
    FT = _promote_float(T, pos)
    x = FT(pos[1]) - field.center[1]
    y = FT(pos[2]) - field.center[2]
    r2 = x * x + y * y
    jz = (r2 <= field.a^2 && field.a > zero(T)) ? FT(field.I / (FT(π) * field.a^2)) :
        zero(FT)
    return SVector{3, FT}(zero(FT), zero(FT), jz)
end
