"""
    FLARECoils{T} <: AbstractMagneticField

Magnetic field configuration of the Facility for Laboratory Reconnection
Experiments (FLARE) at Princeton University / PPPL.

FLARE scales up the cylindrical flux-core reconnection concept pioneered by MRX:
- Two internal flux cores at radius `R_pf` and locations `z = ±z_pf`, housing
  Poloidal Field (PF) coils.
- External Equilibrium Field (EF) coils at radius `R_ef` and `z = ±z_ef` to
  control the radial position of the current sheet.
- Optional midplane plasma current `I_p` at radius `R_p` representing the inductive
  reconnection layer.
- Optional out-of-plane (toroidal) guide field `B_guide`.

# Fields
- `R_pf::T`: Flux core major radius [m] (nominal: 0.75 m).
- `z_pf::T`: Flux core axial half-separation [m] (nominal: 0.40 m).
- `I_pf::T`: Current in each PF coil [A].
- `R_ef::T`: Radius of external EF coils [m] (nominal: 1.30 m).
- `z_ef::T`: Axial half-separation of external EF coils [m] (nominal: 0.70 m).
- `I_ef::T`: Current in each EF coil [A].
- `R_p::T`: Reconnection layer plasma major radius [m] (nominal: 0.75 m).
- `I_p::T`: Induced plasma current in the current sheet [A].
- `B_guide::T`: Toroidal guide field at `R_p` [T].
"""
struct FLARECoils{T} <: AbstractMagneticField
    R_pf::T
    z_pf::T
    I_pf::T
    R_ef::T
    z_ef::T
    I_ef::T
    R_p::T
    I_p::T
    B_guide::T

    function FLARECoils(;
            R_pf = 0.75,
            z_pf = 0.4,
            I_pf = 50.0e3,
            R_ef = 1.3,
            z_ef = 0.7,
            I_ef = 15.0e3,
            R_p = 0.75,
            I_p = 0.0,
            B_guide = 0.0
        )
        T = promote_type(
            typeof(R_pf), typeof(z_pf), typeof(I_pf),
            typeof(R_ef), typeof(z_ef), typeof(I_ef),
            typeof(R_p), typeof(I_p), typeof(B_guide)
        )
        T_float = T <: AbstractFloat ? T : Float64
        return new{T_float}(
            T_float(R_pf), T_float(z_pf), T_float(I_pf),
            T_float(R_ef), T_float(z_ef), T_float(I_ef),
            T_float(R_p), T_float(I_p), T_float(B_guide)
        )
    end
end

@inline function (field::FLARECoils{T})(pos) where {T}
    x, y, z = pos[1], pos[2], pos[3]
    R = hypot(x, y)
    (; R_pf, z_pf, I_pf, R_ef, z_ef, I_ef, R_p, I_p, B_guide) = field

    # PF coils in upper and lower flux cores
    cl_pf_up = CurrentLoop(
        R_pf, I_pf, SVector(zero(T), zero(T), z_pf), SVector(zero(T), zero(T), one(T))
    )
    cl_pf_dn = CurrentLoop(
        R_pf, I_pf, SVector(zero(T), zero(T), -z_pf), SVector(zero(T), zero(T), one(T))
    )

    B = getB_loop(pos, cl_pf_up) + getB_loop(pos, cl_pf_dn)

    # EF coils for equilibrium shaping
    if abs(I_ef) > 0.0
        cl_ef_up = CurrentLoop(
            R_ef, I_ef, SVector(zero(T), zero(T), z_ef), SVector(zero(T), zero(T), one(T))
        )
        cl_ef_dn = CurrentLoop(
            R_ef, I_ef, SVector(zero(T), zero(T), -z_ef), SVector(zero(T), zero(T), one(T))
        )
        B += getB_loop(pos, cl_ef_up) + getB_loop(pos, cl_ef_dn)
    end

    # Plasma current sheet at midplane
    if abs(I_p) > 0.0
        cl_p = CurrentLoop(
            R_p, I_p, SVector(zero(T), zero(T), zero(T)), SVector(zero(T), zero(T), one(T))
        )
        B += getB_loop(pos, cl_p)
    end

    # Toroidal guide field B_phi = B_guide * (R_p / R)
    if abs(B_guide) > 0.0 && R > 1.0e-10
        Bphi = B_guide * (R_p / R)
        inv_R = inv(R)
        Bx_guide = -Bphi * (y * inv_R)
        By_guide = Bphi * (x * inv_R)
        B += SVector{3, T}(Bx_guide, By_guide, zero(T))
    end

    return B
end

function vector_potential(field::FLARECoils{T}, pos) where {T}
    (; R_pf, z_pf, I_pf, R_ef, z_ef, I_ef, R_p, I_p) = field
    solver = VectorPotential()

    cl_pf_up = CurrentLoop(
        R_pf, I_pf, SVector(zero(T), zero(T), z_pf), SVector(zero(T), zero(T), one(T))
    )
    cl_pf_dn = CurrentLoop(
        R_pf, I_pf, SVector(zero(T), zero(T), -z_pf), SVector(zero(T), zero(T), one(T))
    )

    A = solve(solver, cl_pf_up, pos) + solve(solver, cl_pf_dn, pos)

    if abs(I_ef) > 0.0
        cl_ef_up = CurrentLoop(
            R_ef, I_ef, SVector(zero(T), zero(T), z_ef), SVector(zero(T), zero(T), one(T))
        )
        cl_ef_dn = CurrentLoop(
            R_ef, I_ef, SVector(zero(T), zero(T), -z_ef), SVector(zero(T), zero(T), one(T))
        )
        A += solve(solver, cl_ef_up, pos) + solve(solver, cl_ef_dn, pos)
    end

    if abs(I_p) > 0.0
        cl_p = CurrentLoop(
            R_p, I_p, SVector(zero(T), zero(T), zero(T)), SVector(zero(T), zero(T), one(T))
        )
        A += solve(solver, cl_p, pos)
    end

    return A
end

"""
    poloidal_flux(field::FLARECoils, pos) -> Real

Compute the poloidal flux function ψ(R, z) = R · A_ϕ(R, z) for the FLARE device.
Contour levels of ψ correspond to magnetic field lines in the poloidal plane.
"""
function poloidal_flux(field::FLARECoils, pos)
    R = hypot(pos[1], pos[2])
    R < 1.0e-10 && return 0.0
    A = vector_potential(field, pos)
    # A_phi = cross(z_hat, R_hat) · A = (-y * Ax + x * Ay) / R
    A_phi = (-pos[2] * A[1] + pos[1] * A[2]) / R
    return R * A_phi
end
