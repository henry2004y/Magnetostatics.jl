"""
    AnalyticalMagnetosphere <: AbstractMagneticField

A generalized analytical model for planetary magnetospheres.
It superposes internal fields, shielding fields, and magnetotail structures, 
optionally including a bow shock and magnetosheath.

# Fields
- `intrinsic_dipole`: Internal planetary dipole field.
- `shielding_dipole`: Image dipole field for magnetopause shielding.
- `tail_field`: Field model for the magnetotail (e.g., `HarrisSheet`).
- `imf_field`: Interplanetary magnetic field (IMF), usually a `UniformField`.
- `magnetosphere_background`: Background field contribution inside the magnetosphere.
- `shielding_dipole_pos`: Position of the image dipole.
- `mp_standoff`: Magnetopause standoff distance.
- `mp_thickness`: Magnetopause transition layer thickness.
- `bs_standoff`: Bow shock standoff distance.
- `bs_flaring`: Bow shock flaring parameter.
- `bs_hyperboloid`: Bow shock hyperboloid parameter.
- `bs_thickness`: Bow shock transition layer thickness.
- `ms_compression`: Magnetosheath compression ratio.
- `has_shock::Bool`: Whether to include the bow shock and magnetosheath regions.
"""
struct AnalyticalMagnetosphere{D1, D2, T, U1, U2, S} <: AbstractMagneticField
    intrinsic_dipole::D1
    shielding_dipole::D2
    tail_field::T
    imf_field::U1
    magnetosphere_background::U2
    shielding_dipole_pos::SVector{3, S}
    mp_standoff::S
    mp_thickness::S
    bs_standoff::S
    bs_flaring::S
    bs_hyperboloid::S
    bs_thickness::S
    ms_compression::S
    has_shock::Bool
end

function AnalyticalMagnetosphere(;
        intrinsic_dipole = NullField(),
        shielding_dipole = NullField(),
        tail_field = NullField(),
        imf_field = NullField(),
        magnetosphere_background = NullField(),
        shielding_dipole_pos = SA[20.0, 0.0, 0.0],
        mp_standoff = 10.0,
        mp_thickness = 1.0,
        bs_standoff = 13.0,
        bs_flaring = 0.04,
        bs_hyperboloid = 2.0,
        bs_thickness = 0.5,
        ms_compression = 3.0,
        has_shock = true
    )
    S = promote_type(
        typeof(mp_standoff), typeof(mp_thickness),
        typeof(bs_standoff), typeof(bs_flaring),
        typeof(bs_hyperboloid), typeof(bs_thickness),
        typeof(ms_compression), eltype(shielding_dipole_pos)
    )
    return AnalyticalMagnetosphere(
        intrinsic_dipole, shielding_dipole, tail_field, imf_field, magnetosphere_background,
        SVector{3, S}(shielding_dipole_pos), S(mp_standoff),
        S(mp_thickness), S(bs_standoff), S(bs_flaring),
        S(bs_hyperboloid), S(bs_thickness), S(ms_compression),
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
    f.intrinsic_dipole(r) + f.shielding_dipole(r .- f.shielding_dipole_pos) +
    f.tail_field(r) + f.magnetosphere_background(r)

@inline inner_potential(f::AnalyticalMagnetosphere, r) =
    vector_potential(f.intrinsic_dipole, r) +
    vector_potential(f.shielding_dipole, r .- f.shielding_dipole_pos) +
    vector_potential(f.tail_field, r) + vector_potential(f.magnetosphere_background, r)

"""
    draped_imf_field(imf, r, mp_standoff)

Return the IMF field at position `r` after being "draped" around a magnetopause
with standoff distance `mp_standoff`. For a `UniformField`, this uses an
analytical potential flow transformation that ensures ∇⋅B = 0 in the magnetosheath.
For other field types, it falls back to the original field.
"""
@inline function draped_imf_field(imf::UniformField, r, mp_standoff)
    B_imf = imf.B0
    bx = B_imf[1] - 0.5 * (B_imf[2] * r[2] + B_imf[3] * r[3]) / mp_standoff
    return SVector(bx, B_imf[2], B_imf[3])
end

@inline draped_imf_field(imf::AbstractMagneticField, r, mp_standoff) = imf(r)

@inline function outer_potential(f::AnalyticalMagnetosphere, r, dist_mp)
    A_sw = vector_potential(f.imf_field, r)
    if !f.has_shock
        return A_sw
    end

    r_prime = SVector(dist_mp, r[2], r[3])
    A_sh = f.ms_compression * vector_potential(f.imf_field, r_prime)

    y2z2 = r[2]^2 + r[3]^2
    rad_term = sqrt(f.bs_flaring * y2z2 + f.bs_hyperboloid^2)
    dist_bs = r[1] - f.bs_standoff + rad_term - f.bs_hyperboloid
    w_bs = 0.5 * (1 - tanh(dist_bs / f.bs_thickness))

    return w_bs * A_sh + (1 - w_bs) * A_sw
end

@inline function outer_region_eval(f::AnalyticalMagnetosphere, r, dist_mp)
    # pristine IMF potential and field
    A_sw = vector_potential(f.imf_field, r)
    B_sw = f.imf_field(r)

    if !f.has_shock
        return A_sw, B_sw
    end

    # Magnetosheath (potential flow transformation)
    r_prime = SVector(dist_mp, r[2], r[3])
    A_sh = f.ms_compression * vector_potential(f.imf_field, r_prime)
    B_sh = f.ms_compression * draped_imf_field(f.imf_field, r, f.mp_standoff)

    # Bow shock geometry
    y2z2 = r[2]^2 + r[3]^2
    rad_term = sqrt(f.bs_flaring * y2z2 + f.bs_hyperboloid^2)
    dist_bs = r[1] - f.bs_standoff + rad_term - f.bs_hyperboloid
    w_bs, dw_bs = tanh_weight(dist_bs, f.bs_thickness)

    # Weighting gradient
    grad_dist_bs = SVector(
        one(eltype(r)), f.bs_flaring * r[2] / rad_term, f.bs_flaring * r[3] / rad_term
    )
    grad_w_bs = dw_bs * grad_dist_bs

    A_out = w_bs * A_sh + (1 - w_bs) * A_sw
    B_out = w_bs * B_sh + (1 - w_bs) * B_sw + cross(grad_w_bs, A_sh - A_sw)

    return A_out, B_out
end

function (f::AnalyticalMagnetosphere)(r)
    # Geometry and weighting function
    x, y, z = r[1], r[2], r[3]
    dist_mp = x - f.mp_standoff + (y^2 + z^2) / (2 * f.mp_standoff)
    w_mp, dw_mp = tanh_weight(dist_mp, f.mp_thickness)

    # Gradient of weighting function
    grad_dist = SVector(one(eltype(r)), y / f.mp_standoff, z / f.mp_standoff)
    grad_w_mp = dw_mp * grad_dist

    # Field and Potential evaluation
    B_in, A_in = inner_field(f, r), inner_potential(f, r)
    A_out, B_out = outer_region_eval(f, r, dist_mp)

    # Total B
    return w_mp * B_in + (1 - w_mp) * B_out + cross(grad_w_mp, A_in - A_out)
end
