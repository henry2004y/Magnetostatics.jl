const _GL8_NODES = SVector{8, Float64}(
    -0.9602898564975363, -0.7966664774136267, -0.525532409916329,
    -0.1834346424956498, 0.1834346424956498, 0.525532409916329,
    0.7966664774136267, 0.9602898564975363
)
const _GL8_WEIGHTS = SVector{8, Float64}(
    0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
    0.362683783378362, 0.362683783378362, 0.3137066458778873,
    0.2223810344533745, 0.1012285362903763
)

"""
    Solenoid{T} <: AbstractMagneticField

Exact magnetic field and vector potential of a finite-length cylindrical solenoid carrying
total ampere-turns `current = turns * I_turn`. Supports zero-thickness sheets and radially
thick windings (`thickness > 0`) via 8-point Gauss-Legendre radial integration.

# Fields
- `radius::T`: Mean winding radius [m].
- `length::T`: Axial length [m].
- `current::T`: Total ampere-turns `N * I` [A].
- `thickness::T`: Radial winding thickness `ΔR >= 0` [m].
- `center::SVector{3, T}`: Geometric center [m].
- `normal::SVector{3, T}`: Unit axis vector in the positive field direction.
"""
struct Solenoid{T} <: AbstractMagneticField
    radius::T
    length::T
    current::T
    thickness::T
    center::SVector{3, T}
    normal::SVector{3, T}

    function Solenoid(
            radius, length, current;
            turns = 1,
            thickness = 0.0,
            center = SVector(0.0, 0.0, 0.0),
            normal = SVector(0.0, 0.0, 1.0)
        )
        T_raw = promote_type(
            typeof(radius), typeof(length), typeof(current),
            typeof(turns), typeof(thickness), eltype(center), eltype(normal)
        )
        T = T_raw <: AbstractFloat ? T_raw : Float64
        R = T(radius)
        L = T(length)
        th = T(thickness)
        R > zero(T) || throw(ArgumentError("radius must be positive"))
        L > zero(T) || throw(ArgumentError("length must be positive"))
        (zero(T) <= th && th < 2 * R) ||
            throw(ArgumentError("thickness must satisfy 0 <= thickness < 2 * radius"))
        n_hat = normalize(SVector{3, T}(normal))
        c = SVector{3, T}(center[1], center[2], center[3])
        return new{T}(R, L, T(turns) * T(current), th, c, n_hat)
    end
end

@inline function _solenoid_sheet_b_rz(
        R::T, L::T, I_tot::T, rho::T, z::T
    ) where {T <: AbstractFloat}
    b0 = T(μ₀) * I_tot / L
    half_l = T(0.5) * L
    zp = z + half_l
    zm = z - half_l
    if rho <= T(1.0e-12) * R
        bz = T(0.5) * b0 * (zp / hypot(zp, R) - zm / hypot(zm, R))
        return zero(T), bz
    end
    rp = R + rho
    rm = R - rho
    dp = hypot(zp, rp)
    dm = hypot(zm, rp)
    kp = hypot(zp, rm) / dp
    km = hypot(zm, rm) / dm
    if kp < T(1.0e-14) || km < T(1.0e-14)
        return zero(T), zero(T)
    end
    ap = R / dp
    am = R / dm
    bp = zp / dp
    bm = zm / dm
    gam = rm / rp
    p = gam * gam
    b_rho = (b0 / T(π)) *
        (ap * _cel(kp, one(T), one(T), -one(T)) - am * _cel(km, one(T), one(T), -one(T)))
    b_z = (b0 * R / (T(π) * rp)) *
        (bp * _cel(kp, p, one(T), gam) - bm * _cel(km, p, one(T), gam))
    return b_rho, b_z
end

@inline function _solenoid_a_kernel(kc::T, p::T) where {T <: AbstractFloat}
    c1 = _cel(kc, one(T), zero(T), one(T))
    iszero(p) && return c1
    return c1 - _cel(kc, p, zero(T), p)
end

@inline function _solenoid_sheet_a_phi(
        R::T, L::T, I_tot::T, rho::T, z::T
    ) where {T <: AbstractFloat}
    rho <= T(1.0e-12) * R && return zero(T)
    b0 = T(μ₀) * I_tot / L
    half_l = T(0.5) * L
    zp = z + half_l
    zm = z - half_l
    rp = R + rho
    rm = R - rho
    dp = hypot(zp, rp)
    dm = hypot(zm, rp)
    kp = hypot(zp, rm) / dp
    km = hypot(zm, rm) / dm
    if kp < T(1.0e-14) || km < T(1.0e-14)
        return zero(T)
    end
    bp = zp / dp
    bm = zm / dm
    gam = rm / rp
    p = gam * gam
    term = bp * _solenoid_a_kernel(kp, p) - bm * _solenoid_a_kernel(km, p)
    return (b0 * R / T(π)) * term
end

@inline function _solenoid_local_cyl(field::Solenoid{T}, pos) where {T}
    @boundscheck length(pos) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    @inbounds r_rel = SVector{3, T}(
        pos[1] - field.center[1], pos[2] - field.center[2], pos[3] - field.center[3]
    )
    z_loc = dot(r_rel, field.normal)
    rho_vec = r_rel - z_loc * field.normal
    rho = norm(rho_vec)
    return rho_vec, rho, z_loc
end

@inline function _gl8_b_rz(
        a::T, b::T, L::T, I_tot::T, rho::T, z::T
    ) where {T <: AbstractFloat}
    mid = T(0.5) * (a + b)
    half = T(0.5) * (b - a)
    br = zero(T)
    bz = zero(T)
    @inbounds for q in 1:8
        r_q = mid + half * T(_GL8_NODES[q])
        w_q = T(0.5) * T(_GL8_WEIGHTS[q])
        br_q, bz_q = _solenoid_sheet_b_rz(r_q, L, I_tot, rho, z)
        br += w_q * br_q
        bz += w_q * bz_q
    end
    return br, bz
end

@inline function _solenoid_thick_b_rz(
        R::T, L::T, I_tot::T, th::T, rho::T, z::T
    ) where {T <: AbstractFloat}
    r_in = R - T(0.5) * th
    r_out = R + T(0.5) * th
    if r_in < rho && rho < r_out
        br1, bz1 = _gl8_b_rz(r_in, rho, L, I_tot, rho, z)
        br2, bz2 = _gl8_b_rz(rho, r_out, L, I_tot, rho, z)
        inv_th = inv(th)
        w1 = (rho - r_in) * inv_th
        w2 = (r_out - rho) * inv_th
        return w1 * br1 + w2 * br2, w1 * bz1 + w2 * bz2
    end
    return _gl8_b_rz(r_in, r_out, L, I_tot, rho, z)
end

@inline function _gl8_a_phi(
        a::T, b::T, L::T, I_tot::T, rho::T, z::T
    ) where {T <: AbstractFloat}
    mid = T(0.5) * (a + b)
    half = T(0.5) * (b - a)
    a_phi = zero(T)
    @inbounds for q in 1:8
        r_q = mid + half * T(_GL8_NODES[q])
        w_q = T(0.5) * T(_GL8_WEIGHTS[q])
        a_phi += w_q * _solenoid_sheet_a_phi(r_q, L, I_tot, rho, z)
    end
    return a_phi
end

@inline function _solenoid_thick_a_phi(
        R::T, L::T, I_tot::T, th::T, rho::T, z::T
    ) where {T <: AbstractFloat}
    r_in = R - T(0.5) * th
    r_out = R + T(0.5) * th
    if r_in < rho && rho < r_out
        a1 = _gl8_a_phi(r_in, rho, L, I_tot, rho, z)
        a2 = _gl8_a_phi(rho, r_out, L, I_tot, rho, z)
        inv_th = inv(th)
        return (rho - r_in) * inv_th * a1 + (r_out - rho) * inv_th * a2
    end
    return _gl8_a_phi(r_in, r_out, L, I_tot, rho, z)
end

function (field::Solenoid{T})(pos) where {T}
    rho_vec, rho, z_loc = _solenoid_local_cyl(field, pos)
    (; radius, length, current, thickness, normal) = field
    if thickness <= zero(T)
        b_rho, b_z = _solenoid_sheet_b_rz(radius, length, current, rho, z_loc)
    else
        b_rho, b_z = _solenoid_thick_b_rz(radius, length, current, thickness, rho, z_loc)
    end
    if rho <= T(1.0e-12) * radius
        return b_z * normal
    end
    rho_hat = rho_vec / rho
    return b_rho * rho_hat + b_z * normal
end

function poloidal_flux(field::Solenoid{T}, pos) where {T}
    _, rho, z_loc = _solenoid_local_cyl(field, pos)
    rho <= T(1.0e-12) * field.radius && return zero(T)
    if field.thickness <= zero(T)
        a_phi = _solenoid_sheet_a_phi(
            field.radius, field.length, field.current, rho, z_loc
        )
        return rho * a_phi
    end
    a_phi = _solenoid_thick_a_phi(
        field.radius, field.length, field.current, field.thickness, rho, z_loc
    )
    return rho * a_phi
end

function current_density(field::Solenoid{T}, pos) where {T}
    (; radius, length, current, thickness, normal) = field
    thickness <= zero(T) && return zero(SVector{3, T})
    rho_vec, rho, z_loc = _solenoid_local_cyl(field, pos)
    half_l = T(0.5) * length
    half_th = T(0.5) * thickness
    if abs(z_loc) <= half_l && abs(rho - radius) <= half_th && rho > T(1.0e-12) * radius
        j_mag = current / (length * thickness)
        phi_hat = cross(normal, rho_vec / rho)
        return j_mag * phi_hat
    end
    return zero(SVector{3, T})
end

"""
    HelmholtzCoil(radius, current; center = SA[0.0, 0.0, 0.0], normal = SA[0.0, 0.0, 1.0])

Construct a [`MagneticMirror`](@ref) pair separated by `distance = radius` to produce a
highly uniform magnetic field near `center`.
"""
function HelmholtzCoil(
        radius, current;
        center = SVector(0.0, 0.0, 0.0), normal = SVector(0.0, 0.0, 1.0)
    )
    return MagneticMirror(radius, radius, current; center = center, normal = normal)
end

"""
    AntiHelmholtzCoil{T} <: AbstractMagneticField

Two coaxial circular coils carrying equal and opposite currents, creating a 3D quadrupole
magnetic trap with a zero-field null at `center`. Defaults to optimal gradient-linearity
spacing `distance = sqrt(3) * radius`.
"""
struct AntiHelmholtzCoil{T} <: AbstractMagneticField
    loop1::CurrentLoop{T}
    loop2::CurrentLoop{T}

    function AntiHelmholtzCoil(
            radius, current;
            distance = sqrt(3) * radius,
            center = SVector(0.0, 0.0, 0.0),
            normal = SVector(0.0, 0.0, 1.0)
        )
        T_raw = promote_type(
            typeof(radius), typeof(current), typeof(distance),
            eltype(center), eltype(normal)
        )
        T = T_raw <: AbstractFloat ? T_raw : Float64
        R = T(radius)
        I_c = T(current)
        half_d = T(0.5) * T(distance)
        c = SVector{3, T}(center[1], center[2], center[3])
        n_hat = normalize(SVector{3, T}(normal))
        l1 = _current_loop_unchecked(R, I_c, c - half_d * n_hat, n_hat)
        l2 = _current_loop_unchecked(R, -I_c, c + half_d * n_hat, n_hat)
        return new{T}(l1, l2)
    end
end

@inline function (field::AntiHelmholtzCoil)(r)
    return getB_loop(r, field.loop1) + getB_loop(r, field.loop2)
end

@inline function current_density(::AntiHelmholtzCoil{T}, r) where {T}
    R = promote_type(T, eltype(r))
    return zero(SVector{3, R <: AbstractFloat ? R : Float64})
end
