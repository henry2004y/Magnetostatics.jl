"""
    Dipole{T} <: AbstractMagneticField

Magnetic field of a dipole moment M at the origin.

# Fields
- `M::SVector{3, T}`: Magnetic dipole moment.
"""
struct Dipole{T} <: AbstractMagneticField
    M::SVector{3, T}
end

@inline function (field::Dipole)(x::T, y::T, z::T) where {T}
    r = SVector(x, y, z)
    r_mag = norm(r)
    if r_mag < 1.0e-10
        return @SVector zeros(T, 3) # Singularity at origin
    end

    n = r / r_mag

    return (μ0_4π / r_mag^3) * (3 * dot(field.M, n) * n - field.M)
end

@inline function (field::Dipole)(r)
    @boundscheck length(r) >= 3 || throw(ArgumentError("Input must have at least 3 elements."))
    @inbounds return field(r[1], r[2], r[3])
end
