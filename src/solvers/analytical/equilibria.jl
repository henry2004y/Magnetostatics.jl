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
        T_float = T <: AbstractFloat ? T : Float64
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
    x, y, z = pos[1], pos[2], pos[3]
    R = hypot(x, y)
    (; R0, kappa, psi0) = eq
    inv_R04 = inv(R0^4)
    term_z = (R^2 * z^2) / (kappa^2)
    term_r = 0.25 * (R^2 - R0^2)^2
    return psi0 * inv_R04 * (term_z + term_r)
end

@inline function (eq::SolovevEquilibrium{T})(pos) where {T}
    x, y, z = pos[1], pos[2], pos[3]
    R = hypot(x, y)
    if R < 1.0e-10
        return SVector{3, T}(zero(T), zero(T), zero(T))
    end

    (; B0, R0, kappa, psi0) = eq
    inv_R04 = inv(R0^4)

    # Poloidal field from derivatives of psi:
    # B_R = -1/R * dpsi/dZ
    # B_Z =  1/R * dpsi/dR
    # B_phi = R0 * B0 / R
    BR = -psi0 * inv_R04 * (2 * z) / (kappa^2)
    BZ = psi0 * inv_R04 * ((2 * z^2) / (kappa^2) + (R^2 - R0^2))
    Bphi = (R0 * B0) / R

    inv_R = inv(R)
    cos_phi = x * inv_R
    sin_phi = y * inv_R

    Bx = BR * cos_phi - Bphi * sin_phi
    By = BR * sin_phi + Bphi * cos_phi
    Bz = BZ

    return SVector{3, T}(Bx, By, Bz)
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
        T_float = T <: AbstractFloat ? T : Float64
        return new{T_float}(T_float(B0), T_float(alpha), SVector{2, T_float}(center))
    end
end

@inline function (field::TaylorState{T})(pos) where {T}
    x = pos[1] - field.center[1]
    y = pos[2] - field.center[2]
    r = hypot(x, y)
    (; B0, alpha) = field

    if r < 1.0e-15
        return SVector{3, T}(zero(T), zero(T), B0)
    end

    inv_r = inv(r)
    cos_theta = x * inv_r
    sin_theta = y * inv_r

    ar = alpha * r
    b_theta = B0 * besselj1(ar)
    bz = B0 * besselj0(ar)

    bx = -b_theta * sin_theta
    by = b_theta * cos_theta

    return SVector{3, T}(bx, by, bz)
end

function current_density(field::TaylorState{T}, pos) where {T}
    B_vec = field(pos)
    return (field.alpha / T(μ₀)) * B_vec
end

function vector_potential(field::TaylorState{T}, pos) where {T}
    B_vec = field(pos)
    return inv(field.alpha) * B_vec
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
        T_float = T <: AbstractFloat ? T : Float64
        return new{T_float}(
            T_float(I), T_float(a), T_float(Bz0), SVector{2, T_float}(center)
        )
    end
end

@inline function (field::ScrewPinch{T})(pos) where {T}
    x = pos[1] - field.center[1]
    y = pos[2] - field.center[2]
    r = hypot(x, y)
    (; I, a, Bz0) = field

    if r < 1.0e-15
        return SVector{3, T}(zero(T), zero(T), Bz0)
    end

    factor = r <= a ? (T(μ₀) * I) / (2π * a^2) : (T(μ₀) * I) / (2π * r^2)
    bx = -factor * y
    by = factor * x
    bz = Bz0

    return SVector{3, T}(bx, by, bz)
end

function current_density(field::ScrewPinch{T}, pos) where {T}
    x = pos[1] - field.center[1]
    y = pos[2] - field.center[2]
    r = hypot(x, y)
    jz = r <= field.a ? field.I / (π * field.a^2) : zero(T)
    return SVector{3, T}(zero(T), zero(T), jz)
end
