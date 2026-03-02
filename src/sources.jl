"""
    Wire{T} <: AbstractCurrentSource

A unified structure for a current carrying wire. 
It can represent a straight wire or a loop depending on the geometry of points.

# Fields
- `points::Vector{SVector{3, T}}`: The points defining the wire geometry.
- `current::T`: The current flowing through the wire in Amperes.
"""
struct Wire{T} <: AbstractCurrentSource
    points::Vector{SVector{3, T}}
    current::T
end

"""
    CurrentLoop{T}

A circular current loop.

# Fields

  - `radius::T`: Radius of the loop [m].
  - `current::T`: Current in the loop [A].
  - `center::SVector{3, T}`: Position of the loop center [m].
  - `normal::SVector{3, T}`: Unit normal vector of the loop.
"""
struct CurrentLoop{T} <: AbstractCurrentSource
    radius::T
    current::T
    center::SVector{3, T}
    normal::SVector{3, T}

    function CurrentLoop(radius, current, center::AbstractVector, normal::AbstractVector)
        T = promote_type(typeof(radius), typeof(current), eltype(center), eltype(normal))
        n_hat = normalize(SVector{3, T}(normal))
        return new{T}(T(radius), T(current), SVector{3, T}(center), n_hat)
    end
end

"""
    SurfaceCurrentMesh{T} <: AbstractCurrentSource

Discretized surface-current distribution on a closed surface.  Each patch `i`
carries a surface current density **K** [A/m] expressed in the local tangent
frame `(tangents1[i], tangents2[i])`.

Created by [`compute_surface_current`](@ref) for a [`ConductingSphere`](@ref) boundary.
"""
struct SurfaceCurrentMesh{T} <: AbstractCurrentSource
    centers::Vector{SVector{3, T}}
    normals::Vector{SVector{3, T}}
    tangents1::Vector{SVector{3, T}}
    tangents2::Vector{SVector{3, T}}
    areas::Vector{T}
    K::Vector{SVector{2, T}}
end
