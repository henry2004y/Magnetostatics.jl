"""
    AbstractBoundary

Abstract supertype for all magnetic field boundary conditions.
"""
abstract type AbstractBoundary end

"""
    OpenBoundary <: AbstractBoundary

Free-space boundary condition: the field decays naturally to zero at infinity.
This is the default for point-source solvers such as `BiotSavart`.
"""
struct OpenBoundary <: AbstractBoundary end

"""
    PeriodicBoundary <: AbstractBoundary

Periodic boundary condition: sources are replicated at integer multiples of the given
period in each Cartesian direction before field summation.

# Fields
- `period::SVector{3,Float64}`: box side lengths (Lx, Ly, Lz) [m].
"""
struct PeriodicBoundary <: AbstractBoundary
    period::SVector{3, Float64}
end

"""
    ConductingWall <: AbstractBoundary

Planar perfectly-conducting wall boundary condition, implemented via the
method of images.  The wall is an infinite plane whose normal points along one
of the Cartesian axes.

# Fields
- `axis::Int`: axis perpendicular to the wall (1 = x, 2 = y, 3 = z).
- `position::Float64`: coordinate of the wall along that axis [m].
"""
struct ConductingWall <: AbstractBoundary
    axis::Int
    position::Float64
end

"""
    ConductingSphere <: AbstractBoundary

Perfectly-conducting spherical boundary condition solved numerically via a
boundary-element method (BEM).  The sphere surface is discretized into patches;
a linear system is solved for the induced surface current **K** that cancels the
normal component of the incident field at every patch centre.

# Fields
- `center::SVector{3,Float64}`: centre of the sphere [m].
- `radius::Float64`: radius of the sphere [m].
- `n_theta::Int`: number of patches in the polar (latitude) direction.
- `n_phi::Int`: number of patches in the azimuthal (longitude) direction.
"""
struct ConductingSphere <: AbstractBoundary
    center::SVector{3, Float64}
    radius::Float64
    n_theta::Int
    n_phi::Int
end

"""
    PrescribedBoundary{F} <: AbstractBoundary

Dirichlet-type background boundary condition: a prescribed magnetic field
`field_func(r)` is added to the source-computed field everywhere in space.

# Fields
- `field_func::F`: callable with signature `(r::SVector{3}) -> SVector{3}`.
"""
struct PrescribedBoundary{F} <: AbstractBoundary
    field_func::F
end
