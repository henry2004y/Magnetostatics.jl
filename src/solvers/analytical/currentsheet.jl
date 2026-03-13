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
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    T = eltype(r)
    return SVector(field.B0 * tanh(r[3] / field.L), zero(T), zero(T))
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
end

function (field::AsymmetricHarrisSheet)(r)
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    T = eltype(r)
    z = r[3]
    L = field.L
    B_x = (field.B1 + field.B2) / 2 + (field.B1 - field.B2) / 2 * tanh(z / L)
    return SVector(B_x, zero(T), zero(T))
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
end

function (field::ForceFreeHarrisSheet)(r)
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    T = eltype(r)
    z = r[3]
    L = field.L
    B_x = field.B0 * tanh(z / L)
    B_y = field.B0 * sech(z / L)
    return SVector(B_x, B_y, zero(T))
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
end

function (field::BifurcatedHarrisSheet)(r)
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    T = eltype(r)
    z = r[3]
    L = field.L
    d = field.d
    B_x = field.B0 / 2 * (tanh((z - d) / L) + tanh((z + d) / L))
    return SVector(B_x, zero(T), zero(T))
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
end

function (field::FadeevIsland)(r)
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    T = eltype(r)
    x, z = r[1], r[3]
    B0, L, Lx, ε = field.B0, field.L, field.Lx, field.ε

    denom = cosh(z / L) + ε * cos(x / Lx)
    B_x = -B0 * sinh(z / L) / denom
    B_z = -B0 * L * ε * sin(x / Lx) / Lx / denom

    return SVector(B_x, zero(T), B_z)
end
