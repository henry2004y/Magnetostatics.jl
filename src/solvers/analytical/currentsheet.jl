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

    function HarrisSheet(B0::Real, L::Real)
        T = float(promote_type(typeof(B0), typeof(L)))
        return new{T}(T(B0), T(L))
    end
    HarrisSheet{T}(B0::Real, L::Real) where {T} = new{T}(T(B0), T(L))
end

function (field::HarrisSheet{T})(r) where {T}
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    (; B0, L) = field
    return SVector{3, FT}(B0 * tanh(r[3] / L), zero(FT), zero(FT))
end

function current_density(field::HarrisSheet{T}, r) where {T}
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    z = r[3]
    (; B0, L) = field
    j_y = B0 / (FT(μ₀) * L) * sech(z / L)^2
    return SVector{3, FT}(zero(FT), j_y, zero(FT))
end

"""
    AsymmetricHarrisSheet{T} <: AbstractMagneticField

Asymmetric Harris current sheet:
Bx(z) = (B1 + B2)/2 + (B1 - B2)/2 * tanh(z/L)
"""
struct AsymmetricHarrisSheet{T} <: AbstractMagneticField
    B1::T
    B2::T
    L::T

    function AsymmetricHarrisSheet(B1::Real, B2::Real, L::Real)
        T = float(promote_type(typeof(B1), typeof(B2), typeof(L)))
        return new{T}(T(B1), T(B2), T(L))
    end
    AsymmetricHarrisSheet{T}(B1::Real, B2::Real, L::Real) where {T} =
        new{T}(T(B1), T(B2), T(L))
end

function (field::AsymmetricHarrisSheet{T})(r) where {T}
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    z = r[3]
    (; B1, B2, L) = field
    B_x = (B1 + B2) / 2 + (B1 - B2) / 2 * tanh(z / L)
    return SVector{3, FT}(B_x, zero(FT), zero(FT))
end

function current_density(field::AsymmetricHarrisSheet{T}, r) where {T}
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    z = r[3]
    (; B1, B2, L) = field
    j_y = (B1 - B2) / (2 * FT(μ₀) * L) * sech(z / L)^2
    return SVector{3, FT}(zero(FT), j_y, zero(FT))
end

"""
    ForceFreeHarrisSheet{T} <: AbstractMagneticField

Force-free Harris current sheet:
Bx(z) = B0 * tanh(z/L)
By(z) = B0 * sech(z/L)
"""
struct ForceFreeHarrisSheet{T} <: AbstractMagneticField
    B0::T
    L::T

    function ForceFreeHarrisSheet(B0::Real, L::Real)
        T = float(promote_type(typeof(B0), typeof(L)))
        return new{T}(T(B0), T(L))
    end
    ForceFreeHarrisSheet{T}(B0::Real, L::Real) where {T} = new{T}(T(B0), T(L))
end

function (field::ForceFreeHarrisSheet{T})(r) where {T}
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    z = r[3]
    (; B0, L) = field
    B_x = B0 * tanh(z / L)
    B_y = B0 * sech(z / L)
    return SVector{3, FT}(B_x, B_y, zero(FT))
end

function current_density(field::ForceFreeHarrisSheet{T}, r) where {T}
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    z = r[3]
    (; B0, L) = field
    sech_zL = sech(z / L)
    j_x = B0 / (FT(μ₀) * L) * sech_zL * tanh(z / L)
    j_y = B0 / (FT(μ₀) * L) * sech_zL^2
    return SVector{3, FT}(j_x, j_y, zero(FT))
end

"""
    BifurcatedHarrisSheet{T} <: AbstractMagneticField

Bifurcated Harris current sheet:
Bx(z) = B0/2 * [tanh((z - d)/L) + tanh((z + d)/L)]
"""
struct BifurcatedHarrisSheet{T} <: AbstractMagneticField
    B0::T
    L::T
    d::T

    function BifurcatedHarrisSheet(B0::Real, L::Real, d::Real)
        T = float(promote_type(typeof(B0), typeof(L), typeof(d)))
        return new{T}(T(B0), T(L), T(d))
    end
    BifurcatedHarrisSheet{T}(B0::Real, L::Real, d::Real) where {T} =
        new{T}(T(B0), T(L), T(d))
end

function (field::BifurcatedHarrisSheet{T})(r) where {T}
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    z = r[3]
    (; B0, L, d) = field
    B_x = B0 / 2 * (tanh((z - d) / L) + tanh((z + d) / L))
    return SVector{3, FT}(B_x, zero(FT), zero(FT))
end

function current_density(field::BifurcatedHarrisSheet{T}, r) where {T}
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    z = r[3]
    (; B0, L, d) = field
    j_y = B0 / (2 * FT(μ₀) * L) * (sech((z - d) / L)^2 + sech((z + d) / L)^2)
    return SVector{3, FT}(zero(FT), j_y, zero(FT))
end

"""
    FadeevIsland{T} <: AbstractMagneticField

Fadeev magnetic island equilibrium:
Ay(x, z) = -B0 * L * ln[cosh(z/L) + ε * cos(x/Lx)]
Bx = ∂Ay/∂z = -B0 * sinh(z/L) / [cosh(z/L) + ε * cos(x/Lx)]
Bz = -∂Ay/∂x = -B0 * L * ε * sin(x/Lx) / Lx / [cosh(z/L) + ε * cos(x/Lx)]
"""
struct FadeevIsland{T} <: AbstractMagneticField
    B0::T
    L::T
    Lx::T
    ε::T

    function FadeevIsland(B0::Real, L::Real, Lx::Real, ε::Real)
        T = float(promote_type(typeof(B0), typeof(L), typeof(Lx), typeof(ε)))
        return new{T}(T(B0), T(L), T(Lx), T(ε))
    end
    FadeevIsland{T}(B0::Real, L::Real, Lx::Real, ε::Real) where {T} =
        new{T}(T(B0), T(L), T(Lx), T(ε))
end

function (field::FadeevIsland{T})(r) where {T}
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    x, z = r[1], r[3]
    (; B0, L, Lx, ε) = field

    denom = cosh(z / L) + ε * cos(x / Lx)
    B_x = -B0 * sinh(z / L) / denom
    B_z = -B0 * L * ε * sin(x / Lx) / Lx / denom

    return SVector{3, FT}(B_x, zero(FT), B_z)
end

function current_density(field::FadeevIsland{T}, r) where {T}
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    x, z = r[1], r[3]
    (; B0, L, Lx, ε) = field

    cosh_zL = cosh(z / L)
    cos_xLx = cos(x / Lx)
    f = cosh_zL + ε * cos_xLx
    prod = cosh_zL * cos_xLx

    # jy = -1/μ₀ * ∇²Ay
    # Following the derivation for general L, Lx:
    term1 = (1 / L) * (1 + ε * prod)
    term2 = (L * ε / Lx^2) * (ε + prod)
    j_y = B0 / (FT(μ₀) * f^2) * (term1 - term2)

    return SVector{3, FT}(zero(FT), j_y, zero(FT))
end
