"""
    UniformField{T} <: AbstractMagneticField

Uniform magnetic field with constant value `B0`.

# Fields
- `B0::SVector{3, T}`: Magnetic field vector.
"""
struct UniformField{T} <: AbstractMagneticField
    B0::SVector{3, T}
end

@inline function (field::UniformField{T})(x, y, z) where {T}
    FT = float(promote_type(T, typeof(x), typeof(y), typeof(z)))
    return SVector{3, FT}(field.B0)
end

@inline function (field::UniformField{T})(r) where {T}
    return SVector{3, _promote_float(T, r)}(field.B0)
end

@inline function current_density(field::UniformField{T}, r) where {T}
    return zero(SVector{3, _promote_float(T, r)})
end
