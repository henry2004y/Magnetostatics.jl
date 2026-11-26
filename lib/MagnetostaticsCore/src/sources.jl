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
        T = float(promote_type(typeof(radius), typeof(current), eltype(center), eltype(normal)))
        n_hat = normalize(SVector{3, T}(normal))
        return new{T}(T(radius), T(current), SVector{3, T}(center), n_hat)
    end

    function CurrentLoop{T}(
            radius::T, current::T, center::SVector{3, T},
            normal::SVector{3, T}, ::Val{false}
        ) where {T}
        return new{T}(radius, current, center, normal)
    end
end

@inline function _current_loop_unchecked(
        radius::T, current::T, center::SVector{3, T}, n_hat::SVector{3, T}
    ) where {T}
    return CurrentLoop{T}(radius, current, center, n_hat, Val(false))
end
