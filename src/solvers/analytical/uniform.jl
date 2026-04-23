"""
    UniformField{T} <: AbstractMagneticField

Uniform magnetic field with constant value `B0`.

# Fields
- `B0::SVector{3, T}`: Magnetic field vector.
"""
struct UniformField{T} <: AbstractMagneticField
    B0::SVector{3, T}
end

@inline function (field::UniformField)(x, y, z)
    return field.B0
end

@inline function (field::UniformField)(r)
    return field.B0
end
