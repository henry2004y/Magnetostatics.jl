"""
    InfiniteWire{T} <: AbstractMagneticField

Magnetic field of an infinite straight current-carrying wire. When `core_radius > 0`,
uses Plummer regularization `ρ² + core_radius²` preserving `∇ × A = B`; otherwise models
a finite solid cylinder of radius `r` with uniform internal current density.

# Fields
- `I::T`: Current in the wire [A].
- `r::T`: Radius of the wire [m].
- `center::SVector{3, T}`: A point on the wire [m].
- `direction::SVector{3, T}`: Unit vector pointing in the direction of the current.
- `core_radius::T`: Optional Plummer softening radius (0 uses solid-wire radius `r`).
- `r_ref::T`: Reference radius for logarithmic vector potential gauge under Plummer core.
- `k::T`: Prefactor (`μ₀ / 4π` in SI, so azimuthal prefactor is `2 * k * I`).
"""
struct InfiniteWire{T} <: AbstractMagneticField
    I::T
    r::T
    center::SVector{3, T}
    direction::SVector{3, T}
    core_radius::T
    r_ref::T
    k::T

    function InfiniteWire(
            I, r,
            center = SVector(0, 0, 0), direction = SVector(0, 0, 1);
            core_radius::Real = 0,
            r_ref::Real = 1,
            k::Union{Real, Nothing} = nothing
        )
        T = k === nothing ?
            float(
                promote_type(
                    typeof(I), typeof(r), eltype(center), eltype(direction),
                    typeof(core_radius), typeof(r_ref)
                )
            ) :
            float(
                promote_type(
                    typeof(I), typeof(r), eltype(center), eltype(direction),
                    typeof(core_radius), typeof(r_ref), typeof(k)
                )
            )
        core_radius >= 0 || throw(ArgumentError("core_radius must be non-negative"))
        r_ref > 0 || throw(ArgumentError("r_ref must be positive"))
        dir_hat = normalize(SVector{3, T}(direction[1], direction[2], direction[3]))
        k_val = k === nothing ? T(μ0_4π) : T(k)
        return new{T}(
            T(I), T(r),
            SVector{3, T}(center[1], center[2], center[3]),
            dir_hat,
            T(core_radius),
            T(r_ref),
            k_val
        )
    end
end

@inline function (field::InfiniteWire{T})(pos) where {T}
    @boundscheck length(pos) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, pos)
    (; I, r, center, direction, core_radius, k) = field
    @inbounds r_pos = SVector{3, FT}(pos[1], pos[2], pos[3]) - center

    z_local = dot(r_pos, direction)
    rho_vec = r_pos - z_local * direction
    rho_sq = dot(rho_vec, rho_vec)
    azimuthal_vec = cross(direction, rho_vec)

    if core_radius > zero(T)
        rho2_eff = rho_sq + FT(core_radius)^2
        return (2 * FT(k) * I / rho2_eff) * azimuthal_vec
    end
    if r > zero(T) && rho_sq <= r^2
        return (2 * FT(k) * I / (r * r)) * azimuthal_vec
    end
    rho_sq < 1.0e-30 && return zero(SVector{3, FT})
    return (2 * FT(k) * I / rho_sq) * azimuthal_vec
end

@inline (field::InfiniteWire)(x, y, z) = field(SVector(x, y, z))

function current_density(field::InfiniteWire{T}, pos) where {T}
    @boundscheck length(pos) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, pos)
    (; I, r, center, direction, core_radius) = field
    core_radius > zero(T) && return zero(SVector{3, FT})
    @inbounds r_pos = SVector{3, FT}(pos[1], pos[2], pos[3]) - center
    z_local = dot(r_pos, direction)
    rho_vec = r_pos - z_local * direction
    rho_sq = dot(rho_vec, rho_vec)
    if rho_sq <= r^2 && r > zero(T)
        return SVector{3, FT}((I / (FT(π) * r^2)) * direction)
    end
    return zero(SVector{3, FT})
end

function vector_potential(field::InfiniteWire{T}, r) where {T}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    (; I, center, direction, core_radius, r_ref, k) = field
    a = field.r
    @inbounds r_pos = SVector{3, FT}(r[1], r[2], r[3]) - center
    z_local = dot(r_pos, direction)
    rho_vec = r_pos - z_local * direction
    rho_sq = dot(rho_vec, rho_vec)
    pref = FT(k) * I
    if core_radius > zero(T)
        rho2_eff = rho_sq + FT(core_radius)^2
        return (-pref * log(rho2_eff / FT(r_ref)^2)) * direction
    end
    if a > zero(T) && rho_sq <= a^2
        return (-pref * (rho_sq / a^2)) * direction
    end
    rho_sq < 1.0e-30 && return zero(SVector{3, FT})
    return (-pref * (one(FT) + log(rho_sq / a^2))) * direction
end

"""
    WireSegment{T} <: AbstractMagneticField

Finite straight current-carrying wire segment from `x1` to `x2` with optional Plummer
core regularization `core_radius`. Satisfies `∇ × A = B` identically for closed loops
even when `core_radius > 0`.

# Fields
- `x1::SVector{3, T}`: Segment start point.
- `x2::SVector{3, T}`: Segment end point.
- `I::T`: Current [A].
- `core_radius::T`: Plummer softening radius (`ρ_eff² = ρ² + core_radius²`).
- `k::T`: Prefactor (defaults to `μ₀ / 4π`; pass `k = 1.0` for normalized code units).
"""
struct WireSegment{T} <: AbstractMagneticField
    x1::SVector{3, T}
    x2::SVector{3, T}
    I::T
    core_radius::T
    k::T
    dir::SVector{3, T}
    length::T

    function WireSegment(
            x1::AbstractVector,
            x2::AbstractVector,
            I;
            core_radius::Real = 0,
            k::Union{Real, Nothing} = nothing
        )
        T = k === nothing ?
            float(promote_type(eltype(x1), eltype(x2), typeof(I), typeof(core_radius))) :
            float(
                promote_type(
                    eltype(x1), eltype(x2), typeof(I), typeof(core_radius), typeof(k)
                )
            )
        core_radius >= 0 || throw(ArgumentError("core_radius must be non-negative"))
        p1 = SVector{3, T}(x1[1], x1[2], x1[3])
        p2 = SVector{3, T}(x2[1], x2[2], x2[3])
        diff = p2 - p1
        L = sqrt(dot(diff, diff))
        L > zero(T) || throw(ArgumentError("WireSegment length must be positive"))
        k_val = k === nothing ? T(μ0_4π) : T(k)
        return new{T}(p1, p2, T(I), T(core_radius), k_val, diff / L, L)
    end
end

@inline function _segment_geom(w::WireSegment{T}, r) where {T}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    @inbounds r_vec = SVector{3, FT}(r[1], r[2], r[3]) - w.x1
    s = dot(r_vec, w.dir)
    p = r_vec - s * w.dir
    rho2_eff = dot(p, p) + FT(w.core_radius)^2
    s1 = -s
    s2 = FT(w.length) - s
    r1 = sqrt(rho2_eff + s1 * s1)
    r2 = sqrt(rho2_eff + s2 * s2)
    return (FT, p, rho2_eff, s1, s2, r1, r2)
end

@inline function (w::WireSegment)(r)
    FT, p, rho2_eff, s1, s2, r1, r2 = _segment_geom(w, r)
    rho2_eff < 1.0e-30 && return zero(SVector{3, FT})
    fac = (FT(w.k) * FT(w.I) / rho2_eff) * (s2 / r2 - s1 / r1)
    return fac * cross(SVector{3, FT}(w.dir), p)
end

@inline function vector_potential(w::WireSegment, r)
    FT, _, rho2_eff, s1, s2, r1, r2 = _segment_geom(w, r)
    num = s2 >= zero(FT) ? (s2 + r2) : (rho2_eff / max(r2 - s2, FT(1.0e-30)))
    den = s1 >= zero(FT) ? (s1 + r1) : (rho2_eff / max(r1 - s1, FT(1.0e-30)))
    a_mag = FT(w.k) * FT(w.I) * log(max(num / max(den, FT(1.0e-30)), FT(1.0e-30)))
    return a_mag * SVector{3, FT}(w.dir)
end

@inline function field_and_potential(w::WireSegment, r)
    FT, p, rho2_eff, s1, s2, r1, r2 = _segment_geom(w, r)
    B = if rho2_eff < 1.0e-30
        zero(SVector{3, FT})
    else
        fac = (FT(w.k) * FT(w.I) / rho2_eff) * (s2 / r2 - s1 / r1)
        fac * cross(SVector{3, FT}(w.dir), p)
    end
    num = s2 >= zero(FT) ? (s2 + r2) : (rho2_eff / max(r2 - s2, FT(1.0e-30)))
    den = s1 >= zero(FT) ? (s1 + r1) : (rho2_eff / max(r1 - s1, FT(1.0e-30)))
    a_mag = FT(w.k) * FT(w.I) * log(max(num / max(den, FT(1.0e-30)), FT(1.0e-30)))
    return (B, a_mag * SVector{3, FT}(w.dir))
end

@inline current_density(::WireSegment{T}, r) where {T} = zero(SVector{3, _promote_float(T, r)})

"""
    CurrentLoopAnalytic{T} <: AbstractMagneticField

Analytical magnetic field of a circular current loop with optional Plummer axial softening
`core_radius` and scaling prefactor `k` (defaults to `μ₀ / 4π`). Preserves `∇ × A = B`
identically even when `core_radius > 0`.
"""
struct CurrentLoopAnalytic{T} <: AbstractMagneticField
    loop::CurrentLoop{T}
    core_radius::T
    k::T

    function CurrentLoopAnalytic(
            loop::CurrentLoop{T};
            core_radius::Real = 0,
            k::Union{Real, Nothing} = nothing
        ) where {T}
        FT = k === nothing ?
            float(promote_type(T, typeof(core_radius))) :
            float(promote_type(T, typeof(core_radius), typeof(k)))
        core_radius >= 0 || throw(ArgumentError("core_radius must be non-negative"))
        lp = FT === T ? loop :
            _current_loop_unchecked(
                FT(loop.radius), FT(loop.current),
                SVector{3, FT}(loop.center), SVector{3, FT}(loop.normal)
            )
        k_val = k === nothing ? FT(μ0_4π) : FT(k)
        return new{FT}(lp, FT(core_radius), k_val)
    end
end

function _loop_BE(
        x, y, z, loop::CurrentLoop{T}, core_radius, k_pref;
        need_B::Bool = true, need_A::Bool = true
    ) where {T}
    FT = float(promote_type(T, typeof(x), typeof(y), typeof(z), typeof(core_radius)))
    (; radius, current, center, normal) = loop
    R = FT(radius)
    I_cur = FT(current)
    kf = k_pref === nothing ? FT(μ0_4π) : FT(k_pref)
    n_hat = SVector{3, FT}(normal)

    r_rel = SVector{3, FT}(x - center[1], y - center[2], z - center[3])
    z_local = dot(r_rel, n_hat)
    rho_vec = r_rel - z_local * n_hat
    rho_sq = dot(rho_vec, rho_vec)
    z_eff2 = z_local * z_local + FT(core_radius)^2

    # Exact on-axis polynomial branch avoids sqrt(rho_sq) so ForwardDiff at rho == 0 works.
    if rho_sq < (FT(1.0e-8) * R)^2
        d2 = R * R + z_eff2
        d2 < FT(1.0e-30) && return (zero(SVector{3, FT}), zero(SVector{3, FT}))
        d = sqrt(d2)
        b_z = (2 * FT(π) * kf * I_cur * R * R) / (d2 * d)
        b_over_rho = (3 * z_local / (2 * d2)) * b_z
        B = need_B ? (b_over_rho * rho_vec + b_z * n_hat) : zero(SVector{3, FT})
        A = need_A ? (FT(0.5) * b_z * cross(n_hat, rho_vec)) : zero(SVector{3, FT})
        return (B, A)
    end

    rho = sqrt(rho_sq)
    beta2 = (R + rho)^2 + z_eff2
    alpha2 = (R - rho)^2 + z_eff2
    alpha2 < (FT(1.0e-12) * R)^2 && return (zero(SVector{3, FT}), zero(SVector{3, FT}))
    beta = sqrt(beta2)
    k_sq = (4 * R * rho) / beta2
    K_val, E_val = _ellipke(k_sq)

    B = if need_B
        fac = (2 * kf * I_cur) / beta
        b_z = fac * (K_val + ((R * R - rho_sq - z_eff2) / alpha2) * E_val)
        b_over_rho = if k_sq < FT(1.0e-4)
            beta5 = beta2 * beta2 * beta
            ((3 * FT(π) * kf * I_cur * R * R * z_local) / beta5) *
                (one(FT) + FT(1.25) * k_sq + FT(175 / 128) * k_sq * k_sq)
        else
            (fac * z_local / rho_sq) * (((R * R + rho_sq + z_eff2) / alpha2) * E_val - K_val)
        end
        b_over_rho * rho_vec + b_z * n_hat
    else
        zero(SVector{3, FT})
    end

    A = if need_A
        a_over_rho = if k_sq < FT(1.0e-4)
            ((FT(π) * kf * I_cur * R * R) / (beta2 * beta)) *
                (one(FT) + FT(0.75) * k_sq + FT(75 / 128) * k_sq * k_sq)
        else
            ((2 * kf * I_cur) / (rho_sq * beta)) *
                ((R * R + rho_sq + z_eff2) * K_val - beta2 * E_val)
        end
        a_over_rho * cross(n_hat, rho_vec)
    else
        zero(SVector{3, FT})
    end

    return (B, A)
end

"""
    getB_loop(r, loop::CurrentLoop; core_radius=0, k=nothing)
    getB_loop(x, y, z, loop::CurrentLoop; core_radius=0, k=nothing)

Calculate the magnetic field `B` [T] at 3D point `r` from a `CurrentLoop`.
"""
@inline function getB_loop(
        r, loop::CurrentLoop; core_radius::Real = 0, k::Union{Real, Nothing} = nothing
    )
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    @inbounds return getB_loop(r[1], r[2], r[3], loop; core_radius = core_radius, k = k)
end

@inline function getB_loop(
        x, y, z, loop::CurrentLoop{T};
        core_radius::Real = 0, k::Union{Real, Nothing} = nothing
    ) where {T}
    B, _ = _loop_BE(x, y, z, loop, core_radius, k; need_B = true, need_A = false)
    return B
end

@inline function (field::CurrentLoopAnalytic)(r)
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    B, _ = _loop_BE(
        r[1], r[2], r[3], field.loop, field.core_radius, field.k;
        need_B = true, need_A = false
    )
    return B
end

@inline function vector_potential(field::CurrentLoopAnalytic, r)
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    _, A = _loop_BE(
        r[1], r[2], r[3], field.loop, field.core_radius, field.k;
        need_B = false, need_A = true
    )
    return A
end

@inline function field_and_potential(field::CurrentLoopAnalytic, r)
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    return _loop_BE(
        r[1], r[2], r[3], field.loop, field.core_radius, field.k;
        need_B = true, need_A = true
    )
end

@inline function current_density(::CurrentLoopAnalytic{T}, pos) where {T}
    return zero(SVector{3, _promote_float(T, pos)})
end

"""
    ToroidalField{T} <: AbstractMagneticField

Axisymmetric toroidal magnetic field `B_φ = k_tor * ρ / (ρ² + core_radius²)` around
`axis` passing through `center`, with vector potential
`A = -0.5 * k_tor * log((ρ² + core_radius²) / r_ref²) * axis` (`∇ × A = B` identically).
"""
struct ToroidalField{T} <: AbstractMagneticField
    k_tor::T
    center::SVector{3, T}
    axis::SVector{3, T}
    core_radius::T
    r_ref::T

    function ToroidalField(
            k_tor;
            center::AbstractVector = SVector(0, 0, 0),
            axis::AbstractVector = SVector(0, 0, 1),
            core_radius::Real = 0,
            r_ref::Real = 1
        )
        T = float(
            promote_type(
                typeof(k_tor), eltype(center), eltype(axis),
                typeof(core_radius), typeof(r_ref)
            )
        )
        core_radius >= 0 || throw(ArgumentError("core_radius must be non-negative"))
        r_ref > 0 || throw(ArgumentError("r_ref must be positive"))
        ax = normalize(SVector{3, T}(axis[1], axis[2], axis[3]))
        return new{T}(
            T(k_tor),
            SVector{3, T}(center[1], center[2], center[3]),
            ax,
            T(core_radius),
            T(r_ref)
        )
    end
end

@inline function _toroidal_r2_eff(tf::ToroidalField{T}, r) where {T}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    @inbounds dr = SVector{3, FT}(
        r[1] - tf.center[1], r[2] - tf.center[2], r[3] - tf.center[3]
    )
    ax = SVector{3, FT}(tf.axis)
    zeta = dot(dr, ax)
    r_perp = dr - zeta * ax
    r2_eff = dot(r_perp, r_perp) + FT(tf.core_radius)^2
    return (FT, ax, r_perp, r2_eff)
end

@inline function (tf::ToroidalField)(r)
    FT, ax, r_perp, r2_eff = _toroidal_r2_eff(tf, r)
    r2_eff < 1.0e-30 && return zero(SVector{3, FT})
    return (FT(tf.k_tor) / r2_eff) * cross(ax, r_perp)
end

@inline function vector_potential(tf::ToroidalField, r)
    FT, ax, _, r2_eff = _toroidal_r2_eff(tf, r)
    a_ax = -FT(0.5) * FT(tf.k_tor) * log(max(r2_eff / FT(tf.r_ref)^2, FT(1.0e-30)))
    return a_ax * ax
end

@inline function field_and_potential(tf::ToroidalField, r)
    FT, ax, r_perp, r2_eff = _toroidal_r2_eff(tf, r)
    B = r2_eff < 1.0e-30 ? zero(SVector{3, FT}) :
        (FT(tf.k_tor) / r2_eff) * cross(ax, r_perp)
    a_ax = -FT(0.5) * FT(tf.k_tor) * log(max(r2_eff / FT(tf.r_ref)^2, FT(1.0e-30)))
    return (B, a_ax * ax)
end

@inline current_density(::ToroidalField{T}, r) where {T} = zero(SVector{3, _promote_float(T, r)})
