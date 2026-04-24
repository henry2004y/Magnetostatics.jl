"""
    AbstractMagneticField

Abstract supertype for all magnetic field representations.
"""
abstract type AbstractMagneticField end

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
