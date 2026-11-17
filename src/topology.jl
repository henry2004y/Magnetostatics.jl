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
    normM = norm(M)
    normM == 0 && return :degenerate
    abs(tr(M)) <= 0.25 * normM || return :degenerate
    M_dense = Matrix(M)
    sv = svdvals(M_dense)
    if sv[3] <= tol * sv[1]
        sv[2] <= tol * sv[1] && return :degenerate
        disc2 = tr(M * M)
        thresh = (tol * normM)^2
        disc2 > thresh && return :X
        disc2 < -thresh && return :O
        return :degenerate
    end
    detM = det(M)
    vals = eigvals(M_dense)
    is_spiral = maximum(abs ∘ imag, vals) > tol * normM
    if detM > 0
        return is_spiral ? :As : :A
    end
    return is_spiral ? :Bs : :B
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

function _degen_null_result(x, res, J, J_mat)
    eigs = all(isfinite, J) ? eigvals(J_mat) : ComplexF64[NaN, NaN, NaN]
    return (;
        point = x, residual = res, status = :degenerate, kind = :degenerate,
        jacobian = J, eigenvalues = SVector{3, ComplexF64}(eigs[1], eigs[2], eigs[3]),
    )
end

function _newton_step_3d(
        J_mat::Matrix{Float64}, sv::Vector{Float64}, f::SVector{3, Float64}
    )
    dx = if sv[3] > 1.0e-10 * sv[1]
        -(J_mat \ Vector(f))
    else
        -pinv(J_mat; rtol = 1.0e-8) * Vector(f)
    end
    return SVector{3, Float64}(dx[1], dx[2], dx[3])
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
        field::AbstractMagneticField, start;
        bounds = nothing, h::Real = 1.0e-5, tol::Real = 1.0e-11, maxiters::Integer = 30
    )
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

    for _ in 1:maxiters
        f = SVector{3, Float64}(field(x))
        res = norm(f)
        J = field_jacobian(field, x; h = h_f)
        J_mat = Matrix(J)
        if !all(isfinite, J)
            return _degen_null_result(x, res, J, J_mat)
        end
        sv = svdvals(J_mat)
        if sv[1] <= 1.0e-14
            return _degen_null_result(x, res, J, J_mat)
        end
        step_vec = _newton_step_3d(J_mat, sv, f)
        x_scale = max(1.0, norm(x))
        if res <= tol && norm(step_vec) <= tol * x_scale
            eigs = eigvals(J_mat)
            return (;
                point = x, residual = res, status = :converged,
                kind = classify_magnetic_null(J), jacobian = J,
                eigenvalues = SVector{3, ComplexF64}(eigs[1], eigs[2], eigs[3]),
            )
        end
        accepted = false
        alpha = 1.0
        for _ in 1:16
            cand = x + alpha * step_vec
            if _in_bounds_3d(cand, b_tup)
                try
                    f_cand = SVector{3, Float64}(field(cand))
                    if norm(f_cand) < res
                        x = cand
                        accepted = true
                        break
                    end
                catch err
                    err isa DomainError || rethrow()
                end
            end
            alpha *= 0.5
        end
        accepted || break
    end

    f_end = SVector{3, Float64}(field(x))
    res_end = norm(f_end)
    J_end = field_jacobian(field, x; h = h_f)
    J_end_mat = Matrix(J_end)
    all(isfinite, J_end) || return _degen_null_result(x, res_end, J_end, J_end_mat)
    sv_end = svdvals(J_end_mat)
    sv_end[1] <= 1.0e-14 && return _degen_null_result(x, res_end, J_end, J_end_mat)
    dx_end = _newton_step_3d(J_end_mat, sv_end, f_end)
    eigs_end = eigvals(J_end_mat)
    converged = res_end <= tol && norm(dx_end) <= 100.0 * tol * max(1.0, norm(x))
    return (;
        point = x, residual = res_end,
        status = converged ? :converged : :maxiters,
        kind = converged ? classify_magnetic_null(J_end) : :unconverged,
        jacobian = J_end,
        eigenvalues = SVector{3, ComplexF64}(eigs_end[1], eigs_end[2], eigs_end[3]),
    )
end

@inline function _poloidal_br_bz(
        field::AbstractMagneticField, R::Float64, Z::Float64, cphi::Float64, sphi::Float64
    )
    B = SVector{3, Float64}(field(SVector(R * cphi, R * sphi, Z)))
    return SVector{2, Float64}(B[1] * cphi + B[2] * sphi, B[3])
end

function _poloidal_jacobian(
        field::AbstractMagneticField, R::Float64, Z::Float64,
        cphi::Float64, sphi::Float64, h::Float64
    )
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
        field::AbstractMagneticField, start;
        phi::Union{Nothing, Real} = nothing, bounds = nothing,
        h::Real = 1.0e-5, tol::Real = 1.0e-11, maxiters::Integer = 30
    )
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

    for _ in 1:maxiters
        f = _poloidal_br_bz(field, R, Z, cphi, sphi)
        res = norm(f)
        J2 = _poloidal_jacobian(field, R, Z, cphi, sphi, h_f)
        if !all(isfinite, J2) ||
                svdvals(Matrix(J2))[2] <= 1.0e-12 * max(svdvals(Matrix(J2))[1], 1.0e-15)
            return (;
                R, Z, point = SVector(R * cphi, R * sphi, Z), residual = res,
                status = :degenerate, kind = :degenerate, jacobian = J2,
            )
        end
        step = -(J2 \ f)
        r_scale = max(1.0, hypot(R, Z))
        if res <= tol && norm(step) <= tol * r_scale
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
                try
                    f_c = _poloidal_br_bz(field, R_c, Z_c, cphi, sphi)
                    if norm(f_c) < res
                        R, Z = R_c, Z_c
                        accepted = true
                        break
                    end
                catch err
                    err isa DomainError || rethrow()
                end
            end
            alpha *= 0.5
        end
        accepted || break
    end

    f_end = _poloidal_br_bz(field, R, Z, cphi, sphi)
    res_end = norm(f_end)
    J_end = _poloidal_jacobian(field, R, Z, cphi, sphi, h_f)
    sv_end = all(isfinite, J_end) ? svdvals(Matrix(J_end)) : zeros(2)
    if !all(isfinite, J_end) || sv_end[2] <= 1.0e-12 * max(sv_end[1], 1.0e-15)
        return (;
            R, Z, point = SVector(R * cphi, R * sphi, Z), residual = res_end,
            status = :degenerate, kind = :degenerate, jacobian = J_end,
        )
    end
    step_end = -(J_end \ f_end)
    converged = res_end <= tol && norm(step_end) <= 100.0 * tol * max(1.0, hypot(R, Z))
    return (;
        R, Z, point = SVector(R * cphi, R * sphi, Z), residual = res_end,
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
