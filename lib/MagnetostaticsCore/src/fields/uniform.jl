"""
    UniformField{T} <: AbstractMagneticField

Uniform magnetic field with constant value `B0`.

# Fields
- `B0::SVector{3, T}`: Magnetic field vector.
"""
struct UniformField{T} <: AbstractMagneticField
    B0::SVector{3, T}
end

function UniformField(B0::AbstractVector)
    T = float(eltype(B0))
    return UniformField{T}(SVector{3, T}(B0[1], B0[2], B0[3]))
end

function UniformField(bx::Real, by::Real, bz::Real)
    T = float(promote_type(typeof(bx), typeof(by), typeof(bz)))
    return UniformField{T}(SVector{3, T}(T(bx), T(by), T(bz)))
end

@inline function (field::UniformField{T})(x, y, z) where {T}
    FT = float(promote_type(T, typeof(x), typeof(y), typeof(z)))
    return SVector{3, FT}(field.B0)
end

@inline function (field::UniformField{T})(r) where {T}
    return SVector{3, _promote_float(T, r)}(field.B0)
end

@inline function current_density(::UniformField{T}, r) where {T}
    return zero(SVector{3, _promote_float(T, r)})
end

@inline function vector_potential(field::UniformField{T}, r) where {T}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    @inbounds r_vec = SVector{3, FT}(r[1], r[2], r[3])
    return FT(0.5) * cross(field.B0, r_vec)
end

@inline function field_and_potential(field::UniformField, r)
    return (field(r), vector_potential(field, r))
end
