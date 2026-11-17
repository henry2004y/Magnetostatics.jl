"""
    SphereMagnet{T} <: AbstractMagneticField

Uniformly magnetized solid sphere with magnetization `M` [A/m], `radius` [m], and
`center` [m]. Produces uniform interior field `2/3 * μ₀ * M` and exact dipole exterior
field for `r > radius`.
"""
struct SphereMagnet{T} <: AbstractMagneticField
    M::SVector{3, T}
    radius::T
    center::SVector{3, T}

    function SphereMagnet(
            M::AbstractVector, radius, center::AbstractVector = SVector(0.0, 0.0, 0.0)
        )
        T_raw = promote_type(eltype(M), typeof(radius), eltype(center))
        T = float(T_raw)
        R = T(radius)
        R > zero(T) || throw(ArgumentError("radius must be positive"))
        return new{T}(
            SVector{3, T}(M[1], M[2], M[3]),
            R,
            SVector{3, T}(center[1], center[2], center[3])
        )
    end
end

@inline function (field::SphereMagnet{T})(r) where {T}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    @inbounds r_rel = SVector{3, FT}(
        r[1] - field.center[1], r[2] - field.center[2], r[3] - field.center[3]
    )
    r_sq = dot(r_rel, r_rel)
    (; M, radius) = field
    if r_sq <= radius^2
        return SVector{3, FT}((FT(2) / FT(3)) * FT(μ₀) * M)
    end
    rmag = sqrt(r_sq)
    n = r_rel / rmag
    m = (FT(4) / FT(3)) * FT(π) * radius^3 * M
    return (FT(μ0_4π) / (r_sq * rmag)) * (3 * dot(m, n) * n - m)
end

@inline function current_density(::SphereMagnet{T}, r) where {T}
    return zero(SVector{3, _promote_float(T, r)})
end

"""
    CylinderMagnet{T} <: AbstractMagneticField

Uniformly axially magnetized solid circular cylinder with `radius` [m], `length` [m],
axial magnetization `M` [A/m] along `normal`, and center `center` [m]. Evaluated via the
exact bound surface current equivalence `K_b = M`.
"""
struct CylinderMagnet{T} <: AbstractMagneticField
    radius::T
    length::T
    M::T
    center::SVector{3, T}
    normal::SVector{3, T}
    _solenoid::Solenoid{T}

    function CylinderMagnet(
            radius, length, M;
            center = SVector(0.0, 0.0, 0.0), normal = SVector(0.0, 0.0, 1.0)
        )
        T_raw = promote_type(
            typeof(radius), typeof(length), typeof(M), eltype(center), eltype(normal)
        )
        T = float(T_raw)
        sol = Solenoid(
            T(radius), T(length), T(M) * T(length);
            center = SVector{3, T}(center[1], center[2], center[3]),
            normal = SVector{3, T}(normal[1], normal[2], normal[3])
        )
        return new{T}(sol.radius, sol.length, T(M), sol.center, sol.normal, sol)
    end
end

@inline function (field::CylinderMagnet)(r)
    return field._solenoid(r)
end

@inline function poloidal_flux(field::CylinderMagnet, r)
    return poloidal_flux(field._solenoid, r)
end

@inline function current_density(::CylinderMagnet{T}, r) where {T}
    return zero(SVector{3, _promote_float(T, r)})
end

"""
    CuboidMagnet{T} <: AbstractMagneticField

Uniformly magnetized rectangular prism (bar magnet) with side lengths
`dimensions = (lx, ly, lz)` [m], magnetization vector `M` [A/m] in local coordinates,
`center` [m], and orientation matrix `rotation`. Evaluated via the exact cancellation-free
8-corner magnetic surface charge summation.
"""
struct CuboidMagnet{T} <: AbstractMagneticField
    dimensions::SVector{3, T}
    M::SVector{3, T}
    center::SVector{3, T}
    rotation::SMatrix{3, 3, T, 9}

    function CuboidMagnet(
            dimensions::AbstractVector,
            M::AbstractVector;
            center::AbstractVector = SVector(0.0, 0.0, 0.0),
            rotation::AbstractMatrix = SMatrix{3, 3, Float64, 9}(I)
        )
        T_raw = promote_type(
            eltype(dimensions), eltype(M), eltype(center), eltype(rotation)
        )
        T = float(T_raw)
        dims = SVector{3, T}(dimensions[1], dimensions[2], dimensions[3])
        all(>(zero(T)), dims) || throw(ArgumentError("dimensions must be positive"))
        R = SMatrix{3, 3, T, 9}(rotation)
        isapprox(R' * R, I; atol = 1.0e-6) ||
            throw(ArgumentError("rotation matrix must be orthogonal"))
        det(R) > zero(T) ||
            throw(ArgumentError("rotation matrix must be proper (det ≈ +1)"))
        return new{T}(
            dims,
            SVector{3, T}(M[1], M[2], M[3]),
            SVector{3, T}(center[1], center[2], center[3]),
            R
        )
    end
end

@inline function _stable_log_minus(R::T, d::T, rp2::T) where {T <: Real}
    fmin = floatmin(float(T))
    return if d > zero(T)
        log(max(rp2 / (R + d), fmin))
    else
        log(max(R - d, fmin))
    end
end

@inline function _cuboid_atan(num::T, den::T) where {T <: Real}
    if abs(den) >= abs(num)
        !(abs(den) > zero(T)) && return zero(T)
        return atan(num / den)
    end
    half_pi = T(π) / 2
    if den > zero(T)
        return (num > zero(T) ? half_pi : -half_pi) - atan(den / num)
    elseif den < zero(T)
        return (num > zero(T) ? -half_pi : half_pi) - atan(den / num)
    end
    return -atan(den / num)
end

function (field::CuboidMagnet{T})(r) where {T}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    @inbounds dr = SVector{3, FT}(
        r[1] - field.center[1], r[2] - field.center[2], r[3] - field.center[3]
    )
    r_loc = field.rotation' * dr
    half = T(0.5) * field.dimensions
    (; M) = field
    corner_tol = T(1.0e-14) * maximum(field.dimensions)

    Hx = zero(FT)
    Hy = zero(FT)
    Hz = zero(FT)
    for sz in (-one(FT), one(FT)), sy in (-one(FT), one(FT)), sx in (-one(FT), one(FT))
        s = sx * sy * sz
        dx = r_loc[1] - sx * half[1]
        dy = r_loc[2] - sy * half[2]
        dz = r_loc[3] - sz * half[3]
        R2 = dx^2 + dy^2 + dz^2
        R2 <= corner_tol^2 && continue
        R = sqrt(R2)
        lx = _stable_log_minus(R, dx, dy^2 + dz^2)
        ly = _stable_log_minus(R, dy, dz^2 + dx^2)
        lz = _stable_log_minus(R, dz, dx^2 + dy^2)
        tx = _cuboid_atan(dy * dz, dx * R)
        ty = _cuboid_atan(dz * dx, dy * R)
        tz = _cuboid_atan(dx * dy, dz * R)
        Hx += s * (M[1] * tx + M[2] * lz + M[3] * ly)
        Hy += s * (M[1] * lz + M[2] * ty + M[3] * lx)
        Hz += s * (M[1] * ly + M[2] * lx + M[3] * tz)
    end
    H = SVector{3, FT}(Hx, Hy, Hz) / (4 * FT(π))
    inside = abs(r_loc[1]) <= half[1] && abs(r_loc[2]) <= half[2] &&
        abs(r_loc[3]) <= half[3]
    B_loc = FT(μ₀) * (inside ? H + M : H)
    return field.rotation * B_loc
end

@inline function current_density(::CuboidMagnet{T}, r) where {T}
    return zero(SVector{3, _promote_float(T, r)})
end
