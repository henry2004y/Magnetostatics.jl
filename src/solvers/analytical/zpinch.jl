"""
    ZPinch{T} <: AbstractMagneticField

Magnetic field of a Z-pinch configuration (cylindrical wire of radius `a` carrying current `I`).

# Fields
- `I::T`: Total axial current [A].
- `a::T`: Pinch radius [m].
- `center::SVector{2, T}`: Center coordinates (x0, y0) in the transverse plane [m].
"""
struct ZPinch{T} <: AbstractMagneticField
    I::T
    a::T
    center::SVector{2, T}

    function ZPinch(I, a; center = SA[0.0, 0.0])
        T = promote_type(typeof(I), typeof(a), eltype(center))
        T_float = float(T)
        return new{T_float}(T_float(I), T_float(a), SVector{2, T_float}(center))
    end
end

@inline function (field::ZPinch{T})(pos) where {T}
    FT = _promote_float(T, pos)
    x = FT(pos[1]) - field.center[1]
    y = FT(pos[2]) - field.center[2]
    r2 = x * x + y * y
    (; I, a) = field

    if a > zero(T) && r2 <= a^2
        factor = (FT(μ₀) * I) / (2 * FT(π) * a^2)
        return SVector{3, FT}(-factor * y, factor * x, zero(FT))
    end
    r2 < 1.0e-30 && return zero(SVector{3, FT})
    factor = (FT(μ₀) * I) / (2 * FT(π) * r2)
    return SVector{3, FT}(-factor * y, factor * x, zero(FT))
end

function current_density(field::ZPinch{T}, pos) where {T}
    FT = _promote_float(T, pos)
    x = FT(pos[1]) - field.center[1]
    y = FT(pos[2]) - field.center[2]
    r2 = x * x + y * y
    jz = (r2 <= field.a^2 && field.a > zero(T)) ? FT(field.I / (FT(π) * field.a^2)) :
        zero(FT)
    return SVector{3, FT}(zero(FT), zero(FT), jz)
end

"""
    getB_zpinch(x, y, z, I, a) -> SVector{3}
    getB_zpinch(r, I, a) -> SVector{3}

Get magnetic field from a Z-pinch configuration.
"""
@inline function getB_zpinch(x, y, z, I, a)
    return getB_zpinch(SVector(x, y, z), I, a)
end

@inline function getB_zpinch(pos, I, a)
    return ZPinch(I, a)(pos)
end
