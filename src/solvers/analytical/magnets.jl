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
        T = T_raw <: AbstractFloat ? T_raw : Float64
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
    @inbounds r_rel = SVector{3, T}(
        r[1] - field.center[1], r[2] - field.center[2], r[3] - field.center[3]
    )
    rmag = norm(r_rel)
    (; M, radius) = field
    if rmag <= radius
        return (T(2) / T(3)) * T(μ₀) * M
    end
    n = r_rel / rmag
    m = (T(4) / T(3)) * T(π) * radius^3 * M
    return (T(μ0_4π) / rmag^3) * (3 * dot(m, n) * n - m)
end

@inline function current_density(::SphereMagnet{T}, r) where {T}
    R = promote_type(T, eltype(r))
    return zero(SVector{3, R <: AbstractFloat ? R : Float64})
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
        T = T_raw <: AbstractFloat ? T_raw : Float64
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
    R = promote_type(T, eltype(r))
    return zero(SVector{3, R <: AbstractFloat ? R : Float64})
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
        T = T_raw <: AbstractFloat ? T_raw : Float64
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

@inline function _stable_log_minus(R::T, d::T, rp2::T) where {T <: AbstractFloat}
    return if d > zero(T)
        log(max(rp2 / (R + d), floatmin(T)))
    else
        log(max(R - d, floatmin(T)))
    end
end

function (field::CuboidMagnet{T})(r) where {T}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    @inbounds dr = SVector{3, T}(
        r[1] - field.center[1], r[2] - field.center[2], r[3] - field.center[3]
    )
    r_loc = field.rotation' * dr
    half = T(0.5) * field.dimensions
    (; M) = field
    corner_tol = T(1.0e-14) * maximum(field.dimensions)

    Hx = zero(T)
    Hy = zero(T)
    Hz = zero(T)
    for sz in (-one(T), one(T)), sy in (-one(T), one(T)), sx in (-one(T), one(T))
        s = sx * sy * sz
        dx = r_loc[1] - sx * half[1]
        dy = r_loc[2] - sy * half[2]
        dz = r_loc[3] - sz * half[3]
        R = sqrt(dx^2 + dy^2 + dz^2)
        R <= corner_tol && continue
        lx = _stable_log_minus(R, dx, dy^2 + dz^2)
        ly = _stable_log_minus(R, dy, dz^2 + dx^2)
        lz = _stable_log_minus(R, dz, dx^2 + dy^2)
        tx = dx == zero(T) ? zero(T) : atan((dy * dz) / (dx * R))
        ty = dy == zero(T) ? zero(T) : atan((dz * dx) / (dy * R))
        tz = dz == zero(T) ? zero(T) : atan((dx * dy) / (dz * R))
        Hx += s * (M[1] * tx + M[2] * lz + M[3] * ly)
        Hy += s * (M[1] * lz + M[2] * ty + M[3] * lx)
        Hz += s * (M[1] * ly + M[2] * lx + M[3] * tz)
    end
    H = SVector{3, T}(Hx, Hy, Hz) / (4 * T(π))
    inside = abs(r_loc[1]) <= half[1] && abs(r_loc[2]) <= half[2] &&
        abs(r_loc[3]) <= half[3]
    B_loc = T(μ₀) * (inside ? H + M : H)
    return field.rotation * B_loc
end

@inline function current_density(::CuboidMagnet{T}, r) where {T}
    R = promote_type(T, eltype(r))
    return zero(SVector{3, R <: AbstractFloat ? R : Float64})
end
