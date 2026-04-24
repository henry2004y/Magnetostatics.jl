"""
    AnalyticalMagnetosphere <: AbstractMagneticField

A generalized analytical model for planetary magnetospheres.
It superposes internal fields, shielding fields, and magnetotail structures, 
optionally including a bow shock and magnetosheath.

# Fields
- `dipole_intrinsic`: Internal planetary dipole field.
- `dipole_mp`: Image dipole field for magnetopause shielding.
- `tail`: Field model for the magnetotail (e.g., `HarrisSheet`).
- `imf`: Interplanetary magnetic field (IMF), usually a `UniformField`.
- `tail_z`: Background field contribution inside the magnetosphere.
- `image_pos`: Position of the image dipole.
- `r_mp`: Magnetopause standoff distance.
- `d_mp`: Magnetopause transition layer thickness.
- `r_bs`: Bow shock standoff distance.
- `a_bs`: Bow shock flaring parameter.
- `b_bs`: Bow shock hyperboloid parameter.
- `d_bs`: Bow shock transition layer thickness.
- `r_c`: Magnetosheath compression ratio.
- `has_shock::Bool`: Whether to include the bow shock and magnetosheath regions.
"""
struct AnalyticalMagnetosphere{D1, D2, T, U1, U2, S} <: AbstractMagneticField
    dipole_intrinsic::D1
    dipole_mp::D2
    tail::T
    imf::U1
    tail_z::U2
    image_pos::SVector{3, S}
    r_mp::S
    d_mp::S
    r_bs::S
    a_bs::S
    b_bs::S
    d_bs::S
    r_c::S
    has_shock::Bool
end

function AnalyticalMagnetosphere(;
        dipole_intrinsic = NullField(),
        dipole_mp = NullField(),
        tail = NullField(),
        imf = NullField(),
        tail_z = NullField(),
        image_pos = SA[20.0, 0.0, 0.0],
        r_mp = 10.0,
        d_mp = 1.0,
        r_bs = 13.0,
        a_bs = 0.04,
        b_bs = 2.0,
        d_bs = 0.5,
        r_c = 3.0,
        has_shock = true
    )
    S = promote_type(
        typeof(r_mp), typeof(d_mp), typeof(r_bs), typeof(a_bs),
        typeof(b_bs), typeof(d_bs), typeof(r_c), eltype(image_pos)
    )
    return AnalyticalMagnetosphere(
        dipole_intrinsic, dipole_mp, tail, imf, tail_z,
        SVector{3, S}(image_pos), S(r_mp), S(d_mp), S(r_bs), S(a_bs), S(b_bs), S(d_bs), S(r_c),
        has_shock
    )
end

@inline function tanh_weight(dist, d)
    t = tanh(dist / d)
    w = 0.5 * (1 - t)
    # dw/ddist = -0.5 * sech^2(dist/d) / d = -0.5 * (1 - t^2) / d
    dw = -0.5 * (1 - t^2) / d
    return w, dw
end

# Helper functions for the magnetosphere components
@inline inner_field(f::AnalyticalMagnetosphere, r) =
    f.dipole_intrinsic(r) + f.dipole_mp(r .- f.image_pos) + f.tail(r) + f.tail_z(r)

@inline inner_potential(f::AnalyticalMagnetosphere, r) =
    vector_potential(f.dipole_intrinsic, r) + vector_potential(f.dipole_mp, r .- f.image_pos) +
    vector_potential(f.tail, r) + vector_potential(f.tail_z, r)

"""
    draped_imf_field(imf, r, r_mp)

Return the IMF field at position `r` after being "draped" around a magnetopause
with standoff distance `r_mp`. For a `UniformField`, this uses an analytical
potential flow transformation that ensures ∇⋅B = 0 in the magnetosheath.
For other field types, it falls back to the original field.
"""
@inline function draped_imf_field(imf::UniformField, r, r_mp)
    B_imf = imf.B0
    bx = B_imf[1] - 0.5 * (B_imf[2] * r[2] + B_imf[3] * r[3]) / r_mp
    return SVector(bx, B_imf[2], B_imf[3])
end

@inline draped_imf_field(imf::AbstractMagneticField, r, r_mp) = imf(r)

@inline function outer_potential(f::AnalyticalMagnetosphere, r, dist_mp)
    A_sw = vector_potential(f.imf, r)
    if !f.has_shock
        return A_sw
    end

    r_prime = SVector(dist_mp, r[2], r[3])
    A_sh = f.r_c * vector_potential(f.imf, r_prime)

    y2z2 = r[2]^2 + r[3]^2
    rad_term = sqrt(f.a_bs * y2z2 + f.b_bs^2)
    dist_bs = r[1] - f.r_bs + rad_term - f.b_bs
    w_bs = 0.5 * (1 - tanh(dist_bs / f.d_bs))

    return w_bs * A_sh + (1 - w_bs) * A_sw
end

@inline function outer_region_eval(f::AnalyticalMagnetosphere, r, dist_mp)
    # pristine IMF potential and field
    A_sw = vector_potential(f.imf, r)
    B_sw = f.imf(r)

    if !f.has_shock
        return A_sw, B_sw
    end

    # Magnetosheath (potential flow transformation)
    r_prime = SVector(dist_mp, r[2], r[3])
    A_sh = f.r_c * vector_potential(f.imf, r_prime)
    B_sh = f.r_c * draped_imf_field(f.imf, r, f.r_mp)

    # Bow shock geometry
    y2z2 = r[2]^2 + r[3]^2
    rad_term = sqrt(f.a_bs * y2z2 + f.b_bs^2)
    dist_bs = r[1] - f.r_bs + rad_term - f.b_bs
    w_bs, dw_bs = tanh_weight(dist_bs, f.d_bs)

    # Weighting gradient
    grad_dist_bs = SVector(one(eltype(r)), f.a_bs * r[2] / rad_term, f.a_bs * r[3] / rad_term)
    grad_w_bs = dw_bs * grad_dist_bs

    A_out = w_bs * A_sh + (1 - w_bs) * A_sw
    B_out = w_bs * B_sh + (1 - w_bs) * B_sw + cross(grad_w_bs, A_sh - A_sw)

    return A_out, B_out
end

function (f::AnalyticalMagnetosphere)(r)
    # Geometry and weighting function
    x, y, z = r[1], r[2], r[3]
    dist_mp = x - f.r_mp + (y^2 + z^2) / (2 * f.r_mp)
    w_mp, dw_mp = tanh_weight(dist_mp, f.d_mp)

    # Gradient of weighting function
    grad_dist = SVector(one(eltype(r)), y / f.r_mp, z / f.r_mp)
    grad_w_mp = dw_mp * grad_dist

    # Field and Potential evaluation
    B_in, A_in = inner_field(f, r), inner_potential(f, r)
    A_out, B_out = outer_region_eval(f, r, dist_mp)

    # Total B
    return w_mp * B_in + (1 - w_mp) * B_out + cross(grad_w_mp, A_in - A_out)
end
