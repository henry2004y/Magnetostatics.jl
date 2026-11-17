"""
    field_jacobian(field::AbstractMagneticField, r; h::Real = 1.0e-5) -> SMatrix{3, 3}

Compute the spatial Jacobian `J[i, j] = ∂B_i / ∂x_j` at 3D Cartesian position `r` using
second-order central finite differences with step size `h`.
"""
function field_jacobian(field::AbstractMagneticField, r; h::Real = 1.0e-5)
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    isfinite(h) && h > 0 || throw(ArgumentError("h must be finite and positive"))
    T = _promote_float(r)
    h_T = T(h)
    inv2h = inv(T(2) * h_T)
    @inbounds p = SVector{3, T}(T(r[1]), T(r[2]), T(r[3]))
    ex = SVector{3, T}(h_T, zero(T), zero(T))
    ey = SVector{3, T}(zero(T), h_T, zero(T))
    ez = SVector{3, T}(zero(T), zero(T), h_T)
    c1 = (SVector{3}(field(p + ex)) - SVector{3}(field(p - ex))) * inv2h
    c2 = (SVector{3}(field(p + ey)) - SVector{3}(field(p - ey))) * inv2h
    c3 = (SVector{3}(field(p + ez)) - SVector{3}(field(p - ez))) * inv2h
    return hcat(c1, c2, c3)
end

@inline function _eigvals_3x3(M::SMatrix{3, 3, Float64, 9})
    if !all(isfinite, M)
        nan_c = ComplexF64(NaN, NaN)
        return SVector{3, ComplexF64}(nan_c, nan_c, nan_c)
    end
    if norm(M - M') <= 1.0e-14 * norm(M)
        ev = eigvals(Symmetric(M))
        return SVector{3, ComplexF64}(
            ComplexF64(ev[1], 0.0), ComplexF64(ev[2], 0.0), ComplexF64(ev[3], 0.0)
        )
    end
    mu = tr(M) / 3.0
    B = M - mu * one(SMatrix{3, 3, Float64, 9})
    s = norm(B)
    s == 0.0 && return SVector{3, ComplexF64}(mu, mu, mu)
    C = B / s
    c2 = tr(C)
    c1 = (C[1, 1] * C[2, 2] - C[1, 2] * C[2, 1]) +
        (C[1, 1] * C[3, 3] - C[1, 3] * C[3, 1]) +
        (C[2, 2] * C[3, 3] - C[2, 3] * C[3, 2])
    c0 = det(C)
    p = c1 - (c2 * c2) / 3.0
    q = -c0 + (c1 * c2) / 3.0 - (2.0 / 27.0) * (c2 * c2 * c2)
    half_q = 0.5 * q
    third_p = p / 3.0
    disc = half_q * half_q + third_p * third_p * third_p
    t1 = if disc > 0.0
        u = cbrt(-half_q - copysign(sqrt(disc), half_q))
        v = u == 0.0 ? 0.0 : -third_p / u
        u + v
    else
        r = sqrt(max(0.0, -third_p))
        theta = atan(sqrt(max(0.0, -disc)), -half_q) / 3.0
        half_q <= 0.0 ? 2.0 * r * cos(theta) : 2.0 * r * cos(theta + 2.0 * pi / 3.0)
    end
    r1 = t1 + c2 / 3.0
    for _ in 1:2
        fval = ((r1 - c2) * r1 + c1) * r1 - c0
        df = (3.0 * r1 - 2.0 * c2) * r1 + c1
        abs(df) > 1.0e-12 && (r1 -= fval / df)
    end
    b1 = r1 - c2
    b0 = r1 * r1 > abs(c0) ? c0 / r1 : c1 + r1 * b1
    half_b1 = -0.5 * b1
    disc_q = half_b1 * half_b1 - b0
    e1 = ComplexF64(mu + s * r1, 0.0)
    if disc_q >= 0.0
        sq = sqrt(disc_q)
        x_far = half_b1 + (half_b1 >= 0.0 ? sq : -sq)
        x_near = x_far == 0.0 ? 0.0 : b0 / x_far
        e2 = ComplexF64(mu + s * x_near, 0.0)
        e3 = ComplexF64(mu + s * x_far, 0.0)
    else
        im_part = s * sqrt(-disc_q)
        re_part = mu + s * half_b1
        e2 = ComplexF64(re_part, -im_part)
        e3 = ComplexF64(re_part, im_part)
    end
    less_eig(a, b) = real(a) < real(b) || (real(a) == real(b) && imag(a) < imag(b))
    less_eig(e2, e1) && ((e1, e2) = (e2, e1))
    less_eig(e3, e2) && ((e2, e3) = (e3, e2))
    less_eig(e2, e1) && ((e1, e2) = (e2, e1))
    return SVector{3, ComplexF64}(e1, e2, e3)
end

@inline function _classify_null_3d(
        M::SMatrix{3, 3, Float64, 9}, sv::SVector{3, Float64},
        vals::SVector{3, ComplexF64}, tol::Float64
    )
    all(isfinite, M) || return :degenerate
    normM = norm(M)
    normM == 0.0 && return :degenerate
    abs(tr(M)) <= 0.25 * normM || return :degenerate
    if sv[3] <= tol * sv[1]
        sv[2] <= tol * sv[1] && return :degenerate
        disc2 = tr(M * M)
        thresh = (tol * normM)^2
        disc2 > thresh && return :X
        disc2 < -thresh && return :O
        return :degenerate
    end
    detM = det(M)
    max_im = max(abs(imag(vals[1])), abs(imag(vals[2])), abs(imag(vals[3])))
    is_spiral = max_im > tol * normM
    if detM > 0.0
        return is_spiral ? :As : :A
    end
    return is_spiral ? :Bs : :B
end

"""
    classify_magnetic_null(J::AbstractMatrix; tol::Real = 1.0e-6) -> Symbol

Classify a magnetic null from its `3 × 3` field Jacobian `J[i, j] = ∂B_i / ∂x_j` according
to the Parnell et al. (1996) eigenvalue taxonomy:
- `:A`, `:B`: 3D radial nulls (`det(J) > 0` for fan inflow / spine outflow Type A;
  `det(J) < 0` for spine inflow / fan outflow Type B).
- `:As`, `:Bs`: 3D spiral nulls with complex-conjugate fan eigenvalues.
- `:X`, `:O`: 2D hyperbolic X-point or elliptic O-point null line (`rank(J) == 2`).
- `:degenerate`: vanishing, rank-1, parabolic, or strongly non-solenoidal Jacobian.
"""
function classify_magnetic_null(J::AbstractMatrix; tol::Real = 1.0e-6)
    size(J) == (3, 3) || throw(ArgumentError("Jacobian must be a 3 × 3 matrix"))
    isfinite(tol) && tol >= 0 || throw(ArgumentError("tol must be finite and nonnegative"))
    M = SMatrix{3, 3, Float64, 9}(J)
    all(isfinite, M) || return :degenerate
    return _classify_null_3d(M, svdvals(M), _eigvals_3x3(M), Float64(tol))
end

@inline function _in_bounds_3d(p::SVector{3, Float64}, bounds)
    bounds === nothing && return true
    return bounds[1] <= p[1] <= bounds[2] &&
        bounds[3] <= p[2] <= bounds[4] &&
        bounds[5] <= p[3] <= bounds[6]
end

@inline function _in_bounds_2d(R::Float64, Z::Float64, bounds)
    bounds === nothing && return true
    return bounds[1] <= R <= bounds[2] && bounds[3] <= Z <= bounds[4]
end

function _degen_null_result(x, res, J::SMatrix{3, 3, Float64, 9})
    eigs = _eigvals_3x3(J)
    return (;
        point = x, residual = res, status = :degenerate, kind = :degenerate,
        jacobian = J, eigenvalues = eigs,
    )
end

@inline function _newton_step_3d(
        J::SMatrix{3, 3, Float64, 9}, sv::SVector{3, Float64}, f::SVector{3, Float64}
    )
    sv[3] > 1.0e-10 * sv[1] && return -(J \ f)
    F = svd(J)
    cutoff = 1.0e-8 * F.S[1]
    s_inv = SVector{3, Float64}(
        F.S[1] > cutoff ? inv(F.S[1]) : 0.0,
        F.S[2] > cutoff ? inv(F.S[2]) : 0.0,
        F.S[3] > cutoff ? inv(F.S[3]) : 0.0,
    )
    return -(F.V * (s_inv .* (F.U' * f)))
end

@inline function _try_eval_3d(field::F, cand::SVector{3, Float64}) where {F}
    try
        return SVector{3, Float64}(field(cand)), true
    catch err
        err isa DomainError && return zero(SVector{3, Float64}), false
        rethrow()
    end
end

"""
    locate_magnetic_null(field::AbstractMagneticField, start;
                         bounds = nothing, h = 1.0e-5, tol = 1.0e-11, maxiters = 30)

Locate a 3D magnetic null (`B = 0`) starting from `start` via damped Newton-Raphson with
Armijo backtracking, optional box `bounds = (xmin, xmax, ymin, ymax, zmin, zmax)`, and
rank-aware pseudoinverse fallback for 2D null lines.

Returns a named tuple `(; point, residual, status, kind, jacobian, eigenvalues)` where
`status` is `:converged`, `:degenerate`, or `:maxiters`, and `kind` is the topological class
from [`classify_magnetic_null`](@ref) (or `:unconverged`).
"""
function locate_magnetic_null(
        field::F, start;
        bounds = nothing, h::Real = 1.0e-5, tol::Real = 1.0e-11, maxiters::Integer = 30
    ) where {F <: AbstractMagneticField}
    length(start) >= 3 || throw(ArgumentError("start must have 3 coordinates"))
    isfinite(h) && h > 0 || throw(ArgumentError("h must be finite and positive"))
    isfinite(tol) && tol > 0 || throw(ArgumentError("tol must be finite and positive"))
    maxiters >= 1 || throw(ArgumentError("maxiters must be >= 1"))
    b_tup = if bounds === nothing
        nothing
    else
        length(bounds) == 6 && all(isfinite, bounds) &&
            all(bounds[2k] > bounds[2k - 1] for k in 1:3) ||
            throw(ArgumentError("bounds must be 6 increasing finite values"))
        NTuple{6, Float64}(Float64.(Tuple(bounds)))
    end
    x = SVector{3, Float64}(Float64(start[1]), Float64(start[2]), Float64(start[3]))
    all(isfinite, x) || throw(ArgumentError("start must be finite"))
    _in_bounds_3d(x, b_tup) || throw(ArgumentError("start must lie inside bounds"))
    h_f = Float64(h)
    tol_f = Float64(tol)

    f = SVector{3, Float64}(field(x))
    res = norm(f)
    for _ in 1:maxiters
        J = field_jacobian(field, x; h = h_f)
        all(isfinite, J) || return _degen_null_result(x, res, J)
        sv = svdvals(J)
        sv[1] <= 1.0e-14 && return _degen_null_result(x, res, J)
        step_vec = _newton_step_3d(J, sv, f)
        x_scale = max(1.0, norm(x))
        if res <= tol_f && norm(step_vec) <= tol_f * x_scale
            eigs = _eigvals_3x3(J)
            return (;
                point = x, residual = res, status = :converged,
                kind = _classify_null_3d(J, sv, eigs, 1.0e-6),
                jacobian = J, eigenvalues = eigs,
            )
        end
        accepted = false
        alpha = 1.0
        for _ in 1:16
            cand = x + alpha * step_vec
            if _in_bounds_3d(cand, b_tup)
                f_cand, ok = _try_eval_3d(field, cand)
                if ok
                    res_cand = norm(f_cand)
                    if res_cand < res
                        x = cand
                        f = f_cand
                        res = res_cand
                        accepted = true
                        break
                    end
                end
            end
            alpha *= 0.5
        end
        accepted || break
    end

    J_end = field_jacobian(field, x; h = h_f)
    all(isfinite, J_end) || return _degen_null_result(x, res, J_end)
    sv_end = svdvals(J_end)
    sv_end[1] <= 1.0e-14 && return _degen_null_result(x, res, J_end)
    dx_end = _newton_step_3d(J_end, sv_end, f)
    eigs_end = _eigvals_3x3(J_end)
    converged = res <= tol_f && norm(dx_end) <= 100.0 * tol_f * max(1.0, norm(x))
    kind_end = converged ?
        _classify_null_3d(J_end, sv_end, eigs_end, 1.0e-6) : :unconverged
    return (;
        point = x, residual = res,
        status = converged ? :converged : :maxiters,
        kind = kind_end,
        jacobian = J_end, eigenvalues = eigs_end,
    )
end

@inline function _poloidal_br_bz(
        field::F, R::Float64, Z::Float64, cphi::Float64, sphi::Float64
    ) where {F <: AbstractMagneticField}
    B = SVector{3, Float64}(field(SVector(R * cphi, R * sphi, Z)))
    return SVector{2, Float64}(B[1] * cphi + B[2] * sphi, B[3])
end

@inline function _try_poloidal_br_bz(
        field::F, R::Float64, Z::Float64, cphi::Float64, sphi::Float64
    ) where {F <: AbstractMagneticField}
    try
        return _poloidal_br_bz(field, R, Z, cphi, sphi), true
    catch err
        err isa DomainError && return zero(SVector{2, Float64}), false
        rethrow()
    end
end

function _poloidal_jacobian(
        field::F, R::Float64, Z::Float64,
        cphi::Float64, sphi::Float64, h::Float64
    ) where {F <: AbstractMagneticField}
    inv2h = inv(2.0 * h)
    cR = (
        _poloidal_br_bz(field, R + h, Z, cphi, sphi) -
            _poloidal_br_bz(field, R - h, Z, cphi, sphi)
    ) * inv2h
    cZ = (
        _poloidal_br_bz(field, R, Z + h, cphi, sphi) -
            _poloidal_br_bz(field, R, Z - h, cphi, sphi)
    ) * inv2h
    return SMatrix{2, 2, Float64, 4}(cR[1], cR[2], cZ[1], cZ[2])
end

function _classify_poloidal(J2::SMatrix{2, 2, Float64, 4})
    return det(J2) < 0.0 ? :X_point : :O_point
end

"""
    locate_poloidal_critical_point(field::AbstractMagneticField, start;
                                   phi = nothing, bounds = nothing,
                                   h = 1.0e-5, tol = 1.0e-11, maxiters = 30)

Locate a poloidal critical point (`B_R = 0, B_Z = 0`, i.e., a magnetic axis O-point or
poloidal X-point) in the cylindrical plane at azimuth `phi`. `start` may be `(R, Z)` or
3D `(x, y, z)`; when `phi === nothing`, a 3D `start` uses `atan(y, x)` and a 2D `start`
uses `phi = 0.0`. `bounds`, if given, is `(Rmin, Rmax, Zmin, Zmax)` with `Rmin > 0`.

Returns `(; R, Z, point, residual, status, kind, jacobian)` where `kind` is `:O_point`,
`:X_point`, `:degenerate`, or `:unconverged`.
"""
function locate_poloidal_critical_point(
        field::F, start;
        phi::Union{Nothing, Real} = nothing, bounds = nothing,
        h::Real = 1.0e-5, tol::Real = 1.0e-11, maxiters::Integer = 30
    ) where {F <: AbstractMagneticField}
    len = length(start)
    len == 2 || len == 3 ||
        throw(ArgumentError("start must have 2 coordinates (R, Z) or 3 coordinates"))
    isfinite(h) && h > 0 || throw(ArgumentError("h must be finite and positive"))
    isfinite(tol) && tol > 0 || throw(ArgumentError("tol must be finite and positive"))
    maxiters >= 1 || throw(ArgumentError("maxiters must be >= 1"))
    b_tup = if bounds === nothing
        nothing
    else
        length(bounds) == 4 && all(isfinite, bounds) &&
            bounds[1] > 0 && bounds[2] > bounds[1] && bounds[4] > bounds[3] ||
            throw(ArgumentError("bounds must be (Rmin, Rmax, Zmin, Zmax) with Rmin > 0"))
        NTuple{4, Float64}(Float64.(Tuple(bounds)))
    end
    R, Z, phi_val = if len == 2
        Float64(start[1]), Float64(start[2]), phi === nothing ? 0.0 : Float64(phi)
    else
        ang = phi === nothing ? atan(Float64(start[2]), Float64(start[1])) : Float64(phi)
        hypot(Float64(start[1]), Float64(start[2])), Float64(start[3]), ang
    end
    isfinite(R) && isfinite(Z) && isfinite(phi_val) && R > Float64(h) ||
        throw(ArgumentError("start R must be finite and > h"))
    _in_bounds_2d(R, Z, b_tup) || throw(ArgumentError("start must lie inside bounds"))
    sphi, cphi = sincos(phi_val)
    h_f = Float64(h)
    tol_f = Float64(tol)

    f = _poloidal_br_bz(field, R, Z, cphi, sphi)
    res = norm(f)
    for _ in 1:maxiters
        J2 = _poloidal_jacobian(field, R, Z, cphi, sphi, h_f)
        if !all(isfinite, J2)
            return (;
                R, Z, point = SVector(R * cphi, R * sphi, Z), residual = res,
                status = :degenerate, kind = :degenerate, jacobian = J2,
            )
        end
        sv = svdvals(J2)
        if sv[2] <= 1.0e-12 * max(sv[1], 1.0e-15)
            return (;
                R, Z, point = SVector(R * cphi, R * sphi, Z), residual = res,
                status = :degenerate, kind = :degenerate, jacobian = J2,
            )
        end
        step = -(J2 \ f)
        r_scale = max(1.0, hypot(R, Z))
        if res <= tol_f && norm(step) <= tol_f * r_scale
            return (;
                R, Z, point = SVector(R * cphi, R * sphi, Z), residual = res,
                status = :converged, kind = _classify_poloidal(J2), jacobian = J2,
            )
        end
        accepted = false
        alpha = 1.0
        for _ in 1:16
            R_c = R + alpha * step[1]
            Z_c = Z + alpha * step[2]
            if R_c > h_f && _in_bounds_2d(R_c, Z_c, b_tup)
                f_c, ok = _try_poloidal_br_bz(field, R_c, Z_c, cphi, sphi)
                if ok
                    res_c = norm(f_c)
                    if res_c < res
                        R, Z = R_c, Z_c
                        f = f_c
                        res = res_c
                        accepted = true
                        break
                    end
                end
            end
            alpha *= 0.5
        end
        accepted || break
    end

    J_end = _poloidal_jacobian(field, R, Z, cphi, sphi, h_f)
    if !all(isfinite, J_end)
        return (;
            R, Z, point = SVector(R * cphi, R * sphi, Z), residual = res,
            status = :degenerate, kind = :degenerate, jacobian = J_end,
        )
    end
    sv_end = svdvals(J_end)
    if sv_end[2] <= 1.0e-12 * max(sv_end[1], 1.0e-15)
        return (;
            R, Z, point = SVector(R * cphi, R * sphi, Z), residual = res,
            status = :degenerate, kind = :degenerate, jacobian = J_end,
        )
    end
    step_end = -(J_end \ f)
    converged = res <= tol_f && norm(step_end) <= 100.0 * tol_f * max(1.0, hypot(R, Z))
    return (;
        R, Z, point = SVector(R * cphi, R * sphi, Z), residual = res,
        status = converged ? :converged : :maxiters,
        kind = converged ? _classify_poloidal(J_end) : :unconverged,
        jacobian = J_end,
    )
end

function _basis_eval_matrix_3d(
        field::FieldReconstruction3D{T}, point::SVector{3, T}
    ) where {T}
    ncenters = size(field.centers, 1)
    n_rbf = 3 + 3 * ncenters
    n_vac = _harmonic_count(field.vacuum_degree)
    n_dip = size(field.dipole_positions, 1)
    n_src = length(field.source_fields)
    ncoeff = n_rbf + n_vac + 3 * n_dip + n_src
    E = zeros(T, 3, ncoeff)
    for k in 1:3
        E[k, k] = one(T)
    end
    scale_T = field.length_scale
    if ncenters > 0
        scratch = _KernelScratch3D{T}()
        _collect_kernel_3d!(
            scratch, point, field.centers, scale_T, field.kernel,
            field.domain, field.grid_bins
        )
        @inbounds for idx in eachindex(scratch.indices)
            g = scratch.grads[idx]
            off = 3 + 3 * (scratch.indices[idx] - 1)
            E[2, off + 1] = scale_T * g[3]
            E[3, off + 1] = -scale_T * g[2]
            E[1, off + 2] = -scale_T * g[3]
            E[3, off + 2] = scale_T * g[1]
            E[1, off + 3] = scale_T * g[2]
            E[2, off + 3] = -scale_T * g[1]
        end
    end
    if field.vacuum_degree >= 1
        _foreach_harmonic_grad_3d(
            point, field.vacuum_center, field.vacuum_scale, field.vacuum_degree
        ) do idx, grad
            @inbounds for m in 1:3
                E[m, n_rbf + idx] = grad[m]
            end
        end
    end
    offset_dip = n_rbf + n_vac
    if n_dip > 0
        @inbounds for j in 1:n_dip
            rd = SVector{3, T}(
                field.dipole_positions[j, 1],
                field.dipole_positions[j, 2],
                field.dipole_positions[j, 3]
            )
            M_j = _dipole_basis_matrix(point, rd, _dipole_scale(field.domain, rd))
            col_off = offset_dip + 3 * (j - 1)
            for k in 1:3, m in 1:3
                E[m, col_off + k] = M_j[m, k]
            end
        end
    end
    offset_src = offset_dip + 3 * n_dip
    for (idx, sf) in enumerate(field.source_fields)
        b_sf = _eval_source_term(Val(false), sf, point)
        @inbounds for m in 1:3
            E[m, offset_src + idx] = b_sf[m]
        end
    end
    return E
end

function _require_covariance(field)
    hasproperty(field.diagnostics, :coefficient_covariance) &&
        field.diagnostics.coefficient_covariance !== nothing &&
        return field.diagnostics.coefficient_covariance
    throw(ArgumentError("reconstruct with store_covariance = true to compute covariance"))
end

"""
    field_covariance(field, pos) -> SMatrix{3, 3}

Return the `3 × 3` linear-Gaussian prediction covariance matrix `Cov(B(pos))` at Cartesian
position `pos` for a [`FieldReconstruction3D`](@ref) or
[`AxisymmetricFieldReconstruction`](@ref) reconstructed with `store_covariance = true`.
"""
function field_covariance(field::FieldReconstruction3D{T}, pos) where {T}
    Σ_c = _require_covariance(field)
    pt = _check_3d_position(field, pos)
    point = SVector{3, T}(T(pt[1]), T(pt[2]), T(pt[3]))
    E = _basis_eval_matrix_3d(field, point)
    Σ_B = E * Σ_c * transpose(E)
    return SMatrix{3, 3, T, 9}(T(0.5) * (Σ_B + transpose(Σ_B)))
end

function field_covariance(field::AxisymmetricFieldReconstruction{T}, pos) where {T}
    cov_pair = _require_covariance(field)
    len = length(pos)
    len == 2 || len == 3 ||
        throw(ArgumentError("position must have two (R, Z) or three (x, y, z) coordinates"))
    R = len == 2 ? T(pos[1]) : T(hypot(pos[1], pos[2]))
    Z = len == 2 ? T(pos[2]) : T(pos[3])
    Rmin, Rmax, Zmin, Zmax = field.domain
    Rmin <= R <= Rmax && Zmin <= Z <= Zmax ||
        throw(DomainError(pos, "position is outside the reconstruction domain"))
    cphi, sphi = len == 2 ? (one(T), zero(T)) : (T(pos[1]) / R, T(pos[2]) / R)

    Σ_psi = cov_pair.psi
    Σ_tor = cov_pair.toroidal
    n_psi = size(Σ_psi, 1)
    ncenters = size(field.centers, 1)
    e_R = zeros(T, n_psi)
    e_Z = zeros(T, n_psi)
    e_phi = zeros(T, ncenters)
    inv_R = inv(R)
    _foreach_basis_term(field, R, Z) do k, value, dR, dZ, _
        @inbounds begin
            e_R[k] = -dZ * inv_R
            e_Z[k] = dR * inv_R
            k <= ncenters && (e_phi[k] = value)
        end
    end
    v_RR = dot(e_R, Σ_psi * e_R)
    v_RZ = dot(e_R, Σ_psi * e_Z)
    v_ZZ = dot(e_Z, Σ_psi * e_Z)
    v_pp = dot(e_phi, Σ_tor * e_phi)
    Σ_cyl = @SMatrix [
        v_RR zero(T) v_RZ;
        zero(T) v_pp zero(T);
        v_RZ zero(T) v_ZZ
    ]
    Q = @SMatrix [
        cphi -sphi zero(T);
        sphi cphi zero(T);
        zero(T) zero(T) one(T)
    ]
    Σ_B = Q * Σ_cyl * transpose(Q)
    return SMatrix{3, 3, T, 9}(T(0.5) * (Σ_B + transpose(Σ_B)))
end
