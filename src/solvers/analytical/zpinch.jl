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
        T_float = T <: AbstractFloat ? T : Float64
        return new{T_float}(T_float(I), T_float(a), SVector{2, T_float}(center))
    end
end

@inline function (field::ZPinch{T})(pos) where {T}
    x = pos[1] - field.center[1]
    y = pos[2] - field.center[2]
    r = hypot(x, y)
    (; I, a) = field

    if r < 1.0e-15
        return SVector{3, T}(zero(T), zero(T), zero(T))
    end

    factor = r <= a ? (T(μ₀) * I) / (2π * a^2) : (T(μ₀) * I) / (2π * r^2)
    return SVector{3, T}(-factor * y, factor * x, zero(T))
end

function current_density(field::ZPinch{T}, pos) where {T}
    x = pos[1] - field.center[1]
    y = pos[2] - field.center[2]
    r = hypot(x, y)
    jz = r <= field.a ? field.I / (π * field.a^2) : zero(T)
    return SVector{3, T}(zero(T), zero(T), jz)
end

function vector_potential(field::ZPinch{T}, pos) where {T}
    x = pos[1] - field.center[1]
    y = pos[2] - field.center[2]
    r = hypot(x, y)
    (; I, a) = field
    mu0_4pi_I = T(μ0_4π) * I

    az = if r <= a
        -mu0_4pi_I * (r / a)^2
    else
        -mu0_4pi_I * (one(T) + 2 * log(r / a))
    end

    return SVector{3, T}(zero(T), zero(T), az)
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
