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

# Helper functions for the magnetosphere components
@inline inner_field(f::AnalyticalMagnetosphere, r) =
    f.dipole_intrinsic(r) + f.dipole_mp(r .- f.image_pos) + f.tail(r) + f.tail_z(r)

@inline inner_potential(f::AnalyticalMagnetosphere, r) =
    vector_potential(f.dipole_intrinsic, r) + vector_potential(f.dipole_mp, r .- f.image_pos) +
    vector_potential(f.tail, r) + vector_potential(f.tail_z, r)

@inline function sheath_potential(f::AnalyticalMagnetosphere, r)
    dist_mp = r[1] - f.r_mp + (r[2]^2 + r[3]^2) / (2 * f.r_mp)
    r_prime = SVector(dist_mp, r[2], r[3])
    return f.r_c * vector_potential(f.imf, r_prime)
end

@inline function sheath_field(f::AnalyticalMagnetosphere, r)
    # Check if IMF has a constant B0 field for analytical draping
    if hasfield(typeof(f.imf), :B0)
        B_imf = f.imf.B0
        bx = B_imf[1] - 0.5 * (B_imf[2] * r[2] + B_imf[3] * r[3]) / f.r_mp
        return f.r_c * SVector(bx, B_imf[2], B_imf[3])
    else
        return f.r_c * f.imf(r)
    end
end

@inline function outer_potential(f::AnalyticalMagnetosphere, r)
    if !f.has_shock
        return vector_potential(f.imf, r)
    end
    rad_term = sqrt(f.a_bs * (r[2]^2 + r[3]^2) + f.b_bs^2)
    dist_bs = r[1] - f.r_bs + rad_term - f.b_bs
    w_bs = 0.5 * (1 - tanh(dist_bs / f.d_bs))
    return w_bs * sheath_potential(f, r) + (1 - w_bs) * vector_potential(f.imf, r)
end

@inline function outer_field(f::AnalyticalMagnetosphere, r)
    if !f.has_shock
        return f.imf(r)
    end
    rad_term = sqrt(f.a_bs * (r[2]^2 + r[3]^2) + f.b_bs^2)
    dist_bs = r[1] - f.r_bs + rad_term - f.b_bs
    w_bs = 0.5 * (1 - tanh(dist_bs / f.d_bs))

    dw_ddist = -0.5 * sech(dist_bs / f.d_bs)^2 / f.d_bs
    grad_dist = SVector(one(eltype(r)), f.a_bs * r[2] / rad_term, f.a_bs * r[3] / rad_term)
    grad_w_bs = dw_ddist * grad_dist

    A_sh, A_sw = sheath_potential(f, r), vector_potential(f.imf, r)
    B_sh, B_sw = sheath_field(f, r), f.imf(r)

    return w_bs * B_sh + (1 - w_bs) * B_sw + cross(grad_w_bs, A_sh - A_sw)
end

function (f::AnalyticalMagnetosphere)(r)
    # Geometry and weighting function
    x, y, z = r[1], r[2], r[3]
    dist_mp = x - f.r_mp + (y^2 + z^2) / (2 * f.r_mp)
    w_mp = 0.5 * (1 - tanh(dist_mp / f.d_mp))

    # Gradient of weighting function
    dw_ddist = -0.5 * sech(dist_mp / f.d_mp)^2 / f.d_mp
    grad_dist = SVector(one(eltype(r)), y / f.r_mp, z / f.r_mp)
    grad_w_mp = dw_ddist * grad_dist

    # Field and Potential evaluation
    B_in, B_out = inner_field(f, r), outer_field(f, r)
    A_in, A_out = inner_potential(f, r), outer_potential(f, r)

    # Total B
    return w_mp * B_in + (1 - w_mp) * B_out + cross(grad_w_mp, A_in - A_out)
end
