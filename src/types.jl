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

# issue #79: float(...) preserves ForwardDiff.Dual<:Real while mapping Integer to Float64
@inline _pos_eltype(::AbstractVector{T}) where {T} = T
@inline function _pos_eltype(r::Tuple)
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

"""
    vector_potential(field::AbstractMagneticField, r)

Return the magnetic vector potential vector **A** at position `r` for the given field model.
"""
function vector_potential end
