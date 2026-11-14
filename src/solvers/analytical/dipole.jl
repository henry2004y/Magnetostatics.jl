"""
    Dipole{T} <: AbstractMagneticField

Magnetic field of a point magnetic dipole moment `M` located at `center`
(defaults to the origin).

# Fields
- `M::SVector{3, T}`: Magnetic dipole moment [A·m²].
- `center::SVector{3, T}`: Position of the dipole [m].
"""
struct Dipole{T} <: AbstractMagneticField
    M::SVector{3, T}
    center::SVector{3, T}
end

function Dipole(M::AbstractVector)
    T = eltype(M) <: AbstractFloat ? eltype(M) : Float64
    return Dipole{T}(SVector{3, T}(M), zero(SVector{3, T}))
end

function Dipole(M::AbstractVector, center::AbstractVector)
    T_raw = promote_type(eltype(M), eltype(center))
    T = T_raw <: AbstractFloat ? T_raw : Float64
    return Dipole{T}(SVector{3, T}(M), SVector{3, T}(center))
end

@inline function (field::Dipole{T})(x, y, z) where {T}
    R = promote_type(T, typeof(x), typeof(y), typeof(z))
    FT = R <: AbstractFloat ? R : Float64
    r = SVector{3, FT}(
        x - field.center[1], y - field.center[2], z - field.center[3]
    )
    r_mag = norm(r)
    if r_mag < 1.0e-10
        return zero(SVector{3, FT})
    end

    n = r / r_mag
    return (FT(μ0_4π) / r_mag^3) * (3 * dot(field.M, n) * n - field.M)
end

@inline function (field::Dipole)(r)
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    @inbounds return field(r[1], r[2], r[3])
end

@inline function current_density(field::Dipole{T}, r) where {T}
    R = promote_type(T, eltype(r))
    return zero(SVector{3, R <: AbstractFloat ? R : Float64})
end
