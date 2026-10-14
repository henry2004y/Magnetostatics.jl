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
    NullField <: AbstractMagneticField

A magnetic field that is zero everywhere.
"""
struct NullField <: AbstractMagneticField end

(::NullField)(r::AbstractVector{T}) where {T} = zero(SVector{3, T})
(::NullField)(x, y, z) = zero(SVector{3, promote_type(typeof(x), typeof(y), typeof(z))})

"""
    vector_potential(field::AbstractMagneticField, r)

Return the magnetic vector potential vector **A** at position `r` for the given field model.
"""
function vector_potential end
