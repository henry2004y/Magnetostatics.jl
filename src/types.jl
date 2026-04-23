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
    vector_potential(field::AbstractMagneticField, r)

Return the magnetic vector potential vector **A** at position `r` for the given field model.
"""
function vector_potential end
