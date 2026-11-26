"""
    AbstractMagneticField

Abstract supertype for all magnetic field representations.
"""
abstract type AbstractMagneticField end

@inline (field::AbstractMagneticField)(x, y, z) = field(SVector(x, y, z))

@inline function (field::AbstractMagneticField)(dr::AbstractVector, r)
    B = field(r)
    dr[1] = B[1]
    dr[2] = B[2]
    dr[3] = B[3]
    return dr
end

"""
    AbstractCurrentSource

Abstract supertype for all current sources.
"""
abstract type AbstractCurrentSource end

"""
    AbstractSolver

Abstract supertype for all magnetic field solvers.
"""
abstract type AbstractSolver end

"""
    current_density(field::AbstractMagneticField, r)

Return the current density vector [A/m²] at position `r` for the given magnetic field model.
"""
function current_density end

"""
    vector_potential(field::AbstractMagneticField, r)

Return the magnetic vector potential **A** at position `r` for the given field model.
"""
function vector_potential end

"""
    field_and_potential(field::AbstractMagneticField, r) -> (B, A)

Return `(B, A)` at position `r`, where `B = field(r)` and `A = vector_potential(field, r)`.
"""
@inline function field_and_potential(field::AbstractMagneticField, r)
    return (field(r), vector_potential(field, r))
end

@inline function field_and_potential(field::AbstractMagneticField, x, y, z)
    return field_and_potential(field, SVector(x, y, z))
end

"""
    evaluate!(ws, field::AbstractMagneticField, r)

Evaluate `field` at `r` using preallocated scratch workspace `ws`.
"""
function evaluate! end

function poloidal_flux end

@inline _pos_eltype(::AbstractVector{T}) where {T} = T
@inline function _pos_eltype(r::Tuple)
    n = length(r)
    n == 2 && return promote_type(typeof(r[1]), typeof(r[2]))
    return promote_type(typeof(r[1]), typeof(r[2]), typeof(r[3]))
end
@inline _pos_eltype(r) = eltype(r)

@inline _promote_float(::Type{T}, r) where {T} = float(promote_type(T, _pos_eltype(r)))
@inline _promote_float(r) = float(_pos_eltype(r))

"""
    NullField <: AbstractMagneticField

A magnetic field that is zero everywhere.
"""
struct NullField <: AbstractMagneticField end

@inline (::NullField)(r) = zero(SVector{3, _promote_float(r)})
@inline function (::NullField)(x, y, z)
    return zero(SVector{3, float(promote_type(typeof(x), typeof(y), typeof(z)))})
end
@inline current_density(::NullField, r) = zero(SVector{3, _promote_float(r)})
@inline vector_potential(::NullField, r) = zero(SVector{3, _promote_float(r)})
@inline function field_and_potential(::NullField, r)
    z = zero(SVector{3, _promote_float(r)})
    return (z, z)
end
