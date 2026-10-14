"""
    b_unit(B::AbstractVector) -> SVector{3}
    b_unit(field::AbstractMagneticField, r) -> SVector{3}

Compute the unit vector `b = B / |B|` along the magnetic field.
"""
@inline function b_unit(B::AbstractVector{T}) where {T}
    B_mag = norm(B)
    if B_mag < 1.0e-15
        return SVector{3, T}(zero(T), zero(T), zero(T))
    end
    return SVector{3, T}(B[1] / B_mag, B[2] / B_mag, B[3] / B_mag)
end

@inline b_unit(field::AbstractMagneticField, r) = b_unit(field(r))

"""
    field_magnitude(field::AbstractMagneticField, r) -> Real

Compute the magnitude |B| of the magnetic field at position `r`.
"""
@inline field_magnitude(field::AbstractMagneticField, r) = norm(field(r))

"""
    grad_b(field::AbstractMagneticField, r; h=1.0e-5) -> SVector{3}

Compute the gradient ∇|B| of the magnetic field magnitude at position `r`
using central differences with step size `h`.
"""
function grad_b(field::AbstractMagneticField, r; h = 1.0e-5)
    r_vec = SVector{3, Float64}(r[1], r[2], r[3])
    inv_2h = inv(2 * h)

    dx = SVector(h, 0.0, 0.0)
    dy = SVector(0.0, h, 0.0)
    dz = SVector(0.0, 0.0, h)

    gx = (field_magnitude(field, r_vec + dx) - field_magnitude(field, r_vec - dx)) * inv_2h
    gy = (field_magnitude(field, r_vec + dy) - field_magnitude(field, r_vec - dy)) * inv_2h
    gz = (field_magnitude(field, r_vec + dz) - field_magnitude(field, r_vec - dz)) * inv_2h

    return SVector(gx, gy, gz)
end

"""
    magnetic_curvature(field::AbstractMagneticField, r; h=1.0e-5) -> SVector{3}

Compute the field line curvature vector κ = (b · ∇)b at position `r`
using directional central differences along the local unit vector `b`.
"""
function magnetic_curvature(field::AbstractMagneticField, r; h = 1.0e-5)
    r_vec = SVector{3, Float64}(r[1], r[2], r[3])
    b0 = b_unit(field, r_vec)
    norm(b0) < 1.0e-15 && return SVector(0.0, 0.0, 0.0)

    dr = h * b0
    b_plus = b_unit(field, r_vec + dr)
    b_minus = b_unit(field, r_vec - dr)

    return (b_plus - b_minus) / (2 * h)
end

"""
    drift_velocities(field::AbstractMagneticField, r; q=1.0, m=1.0, v_perp=1.0, v_par=1.0, h=1.0e-5)

Compute the guiding center grad-B and curvature drift velocities:
- `v_gradB = (m * v_perp^2) / (2 * q * B^2) * (b × ∇|B|)`
- `v_curv  = (m * v_par^2)  / (q * B)       * (b × κ)`

Returns a NamedTuple `(; v_gradB, v_curv, v_drift)`.
"""
function drift_velocities(
        field::AbstractMagneticField, r;
        q = 1.0, m = 1.0, v_perp = 1.0, v_par = 1.0, h = 1.0e-5
    )
    B_vec = SVector{3, Float64}(field(r))
    B_mag = norm(B_vec)
    B_mag < 1.0e-15 && return (;
        v_gradB = SVector(0.0, 0.0, 0.0),
        v_curv = SVector(0.0, 0.0, 0.0),
        v_drift = SVector(0.0, 0.0, 0.0)
    )

    b = B_vec / B_mag
    ∇B = grad_b(field, r; h = h)
    κ = magnetic_curvature(field, r; h = h)

    coeff_gradB = (m * v_perp^2) / (2 * q * B_mag^2)
    v_gradB = coeff_gradB * cross(b, ∇B)

    coeff_curv = (m * v_par^2) / (q * B_mag)
    v_curv = coeff_curv * cross(b, κ)

    v_drift = v_gradB + v_curv

    return (; v_gradB, v_curv, v_drift)
end

"""
    mirror_ratio(field::AbstractMagneticField, r_min, r_max) -> Real

Compute the magnetic mirror ratio `Rm = |B(r_max)| / |B(r_min)|`.
"""
function mirror_ratio(field::AbstractMagneticField, r_min, r_max)
    B_min = field_magnitude(field, r_min)
    B_max = field_magnitude(field, r_max)
    B_min < 1.0e-15 && throw(ArgumentError("B_min is zero at $r_min"))
    return B_max / B_min
end

"""
    loss_cone_angle(field::AbstractMagneticField, r_min, r_max) -> Real

Compute the loss cone half-angle in radians `θ_loss = asin(1 / sqrt(Rm))`
where `Rm = mirror_ratio(field, r_min, r_max) >= 1`.
"""
function loss_cone_angle(field::AbstractMagneticField, r_min, r_max)
    Rm = mirror_ratio(field, r_min, r_max)
    Rm < 1.0 && throw(ArgumentError("Mirror ratio Rm must be >= 1, got $Rm"))
    return asin(1 / sqrt(Rm))
end

"""
    magnetic_helicity(A::AbstractVector, B::AbstractVector) -> Real

Compute the local magnetic helicity density `h = A · B`.
"""
@inline magnetic_helicity(A::AbstractVector, B::AbstractVector) = dot(A, B)
