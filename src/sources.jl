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

"""
    TorusKnot{T} <: AbstractCurrentSource

A torus knot (p, q), where p and q are coprime integers.

# Fields
- `R::T`: Major radius of the torus.
- `r::T`: Minor radius of the torus.
- `p::Int`: Number of times the knot winds around the symmetry axis.
- `q::Int`: Number of times the knot winds through the "hole" of the torus.
- `current::T`: Current in the wire [A].
- `center::SVector{3, T}`: Center of the torus.
- `normal::SVector{3, T}`: Unit normal vector of the torus axis.
- `up::SVector{3, T}`: Unit vector defining the zero-angle (t=0) direction.
"""
struct TorusKnot{T} <: AbstractCurrentSource
    R::T
    r::T
    p::Int
    q::Int
    current::T
    center::SVector{3, T}
    normal::SVector{3, T}
    up::SVector{3, T}

    function TorusKnot(R, r, p, q, current, center = SA[0, 0, 0], normal = SA[0, 0, 1], up = SA[1, 0, 0])
        T = promote_type(
            typeof(R), typeof(r), typeof(current), eltype(center),
            eltype(normal), eltype(up)
        )
        n_hat = normalize(SVector{3, T}(normal))
        # Ensure 'up' is perpendicular to 'normal'
        u_hat = SVector{3, T}(up)
        u_hat = normalize(u_hat - dot(u_hat, n_hat) * n_hat)
        return new{T}(
            T(R), T(r), Int(p), Int(q), T(current), SVector{3, T}(center),
            n_hat, u_hat
        )
    end
end

"""
    TrefoilKnot(R, r, current, center=[0,0,0], normal=[0,0,1], up=[1,0,0])

A convenience constructor for a (2,3) torus knot, also known as a trefoil knot.
"""
function TrefoilKnot(R, r, current, center = [0, 0, 0], normal = [0, 0, 1], up = [1, 0, 0])
    return TorusKnot(R, r, 2, 3, current, center, normal, up)
end

"""
    MultiSource{S} <: AbstractCurrentSource

Container grouping multiple `AbstractCurrentSource` objects (as a heterogeneous `Tuple` or
`AbstractVector`) so solvers such as `BiotSavart`, `VectorPotential`, and
`compute_surface_current` evaluate their combined field cleanly without concatenating wires.
"""
struct MultiSource{S} <: AbstractCurrentSource
    sources::S

    function MultiSource{S}(sources::S) where {S}
        !isempty(sources) ||
            throw(ArgumentError("MultiSource requires at least one source"))
        return new{S}(sources)
    end
end

function MultiSource(sources::Tuple{Vararg{AbstractCurrentSource}})
    return MultiSource{typeof(sources)}(sources)
end
MultiSource(sources::AbstractCurrentSource...) = MultiSource(sources)
function MultiSource(sources::AbstractVector{<:AbstractCurrentSource})
    return MultiSource{typeof(sources)}(sources)
end
