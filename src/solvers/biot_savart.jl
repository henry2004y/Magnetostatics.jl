"""
    BiotSavart{BC <: AbstractBoundary} <: AbstractSolver

Solver that uses the Biot-Savart law to compute the magnetic field.
The boundary condition `BC` controls how sources outside the physical domain
are treated (e.g., image sources for a conducting wall or periodic lattice sums),
and optional `core_radius > 0` applies Plummer filament softening to avoid
on-wire singularities.

The default constructor `BiotSavart()` uses `OpenBoundary()` with `core_radius = 0.0`.
"""
struct BiotSavart{BC <: AbstractBoundary} <: AbstractSolver
    bc::BC
    core_radius::Float64

    function BiotSavart(
            bc::BC = OpenBoundary(); core_radius::Real = 0.0
        ) where {BC <: AbstractBoundary}
        core_radius < 0 && throw(ArgumentError("core_radius must be >= 0"))
        return new{BC}(bc, Float64(core_radius))
    end
end

@inline function (solver::BiotSavart)(source, r)
    return solve(solver, source, r)
end

@inline _source_eltype(::Wire{T}) where {T} = T
@inline _source_eltype(::CurrentLoop{T}) where {T} = T
@inline _source_eltype(::SurfaceCurrentMesh{T}) where {T} = T
@inline _source_eltype(::MultiSource) = Float64
@inline _source_eltype(ms::MultiSource{<:Tuple}) = _tuple_source_eltype(ms.sources)
@inline _tuple_source_eltype(s::Tuple{A}) where {A} = _source_eltype(s[1])
@inline function _tuple_source_eltype(s::Tuple)
    return promote_type(_source_eltype(s[1]), _tuple_source_eltype(Base.tail(s)))
end
@inline function _source_eltype(sources::AbstractVector{<:AbstractCurrentSource})
    isempty(sources) && return Float64
    T = _source_eltype(first(sources))
    for s in sources
        T = promote_type(T, _source_eltype(s))
    end
    return T
end

const _BiotSavartSource = Union{
    AbstractCurrentSource, AbstractVector{<:AbstractCurrentSource},
}

# ---------------------------------------------------------------------------
# Core kernels (open / free-space)
# ---------------------------------------------------------------------------

function _solve_bs(
        ::OpenBoundary, core_radius::Float64, source::Wire{T},
        r_in::SVector{3}
    ) where {T}
    FT = float(promote_type(T, eltype(r_in)))
    r = SVector{3, FT}(r_in[1], r_in[2], r_in[3])
    B = zero(SVector{3, FT})
    points = source.points
    I = source.current
    rc2 = FT(core_radius * core_radius)

    @inbounds for i in 1:(length(points) - 1)
        r1 = points[i]
        r2 = points[i + 1]
        a = r - r1
        b = r - r2
        dl = r2 - r1
        cross_prod = cross(dl, a)
        cross_sq = dot(cross_prod, cross_prod)
        if core_radius > 0.0
            dl_sq = dot(dl, dl)
            dl_sq < FT(1.0e-30) && continue
            na = sqrt(dot(a, a) + rc2)
            nb = sqrt(dot(b, b) + rc2)
            term = dot(dl, a) / na - dot(dl, b) / nb
            B += (FT(μ0_4π) * I) * (cross_prod / (cross_sq + rc2 * dl_sq)) * term
        else
            cross_sq < FT(1.0e-20) && continue
            na = sqrt(dot(a, a))
            nb = sqrt(dot(b, b))
            term = dot(dl, a) / na - dot(dl, b) / nb
            B += (FT(μ0_4π) * I) * (cross_prod / cross_sq) * term
        end
    end

    return B
end

function _solve_bs(
        ::OpenBoundary, core_radius::Float64, loop::CurrentLoop{T},
        r::SVector{3}
    ) where {T}
    FT = float(promote_type(T, eltype(r)))
    if core_radius <= 0.0
        return SVector{3, FT}(getB_loop(r[1], r[2], r[3], loop))
    end
    (; radius, current, center, normal) = loop
    r_rel = SVector{3, FT}(r[1] - center[1], r[2] - center[2], r[3] - center[3])
    z_local = dot(r_rel, normal)
    rho_vec = r_rel - z_local * normal
    rho_sq = dot(rho_vec, rho_vec)
    rc2 = FT(core_radius * core_radius)
    if rho_sq < (FT(1.0e-8) * radius)^2
        d2 = radius^2 + z_local^2 + rc2
        b_z = FT(μ₀) * current * radius^2 / (2 * d2 * sqrt(d2))
        b_perp_fac = (3 * z_local / (2 * d2)) * b_z
        return b_perp_fac * rho_vec + b_z * normal
    end
    rho = sqrt(rho_sq)
    denom_sq = (radius + rho)^2 + z_local^2 + rc2
    denom_diff_sq = (radius - rho)^2 + z_local^2 + rc2
    k_sq = 4 * radius * rho / denom_sq
    K_val, E_val = _ellipke(k_sq)
    factor = FT(μ₀) * current / (2 * FT(π) * sqrt(denom_sq))
    b_z = factor * (K_val + (radius^2 - rho_sq - z_local^2 - rc2) / denom_diff_sq * E_val)
    b_rho = factor * (z_local / rho) *
        (-K_val + (radius^2 + rho_sq + z_local^2 + rc2) / denom_diff_sq * E_val)
    rho_hat = rho_vec / rho
    return b_rho * rho_hat + b_z * normal
end

function _solve_bs(
        ::OpenBoundary, ::Float64, mesh::SurfaceCurrentMesh{T},
        r_in::SVector{3}
    ) where {T}
    FT = float(promote_type(T, eltype(r_in)))
    r = SVector{3, FT}(r_in[1], r_in[2], r_in[3])
    B = zero(SVector{3, FT})
    @inbounds for i in eachindex(mesh.centers)
        R = r - mesh.centers[i]
        R_sq = dot(R, R)
        R_sq < FT(1.0e-40) && continue
        R_norm = sqrt(R_sq)
        K_vec = mesh.K[i][1] * mesh.tangents1[i] + mesh.K[i][2] * mesh.tangents2[i]
        B += (FT(μ0_4π) * mesh.areas[i] / (R_sq * R_norm)) * cross(K_vec, R)
    end
    return B
end

@inline function _sum_bs_tuple(
        core_radius::Float64, sources::Tuple{S}, r::SVector{3}
    ) where {S}
    return _solve_bs(OpenBoundary(), core_radius, sources[1], r)
end
@inline function _sum_bs_tuple(core_radius::Float64, sources::Tuple, r::SVector{3})
    return _solve_bs(OpenBoundary(), core_radius, sources[1], r) +
        _sum_bs_tuple(core_radius, Base.tail(sources), r)
end

function _solve_bs(
        ::OpenBoundary, core_radius::Float64, ms::MultiSource{<:Tuple}, r::SVector{3}
    )
    return _sum_bs_tuple(core_radius, ms.sources, r)
end

function _solve_bs(
        ::OpenBoundary, core_radius::Float64,
        sources::AbstractVector{<:AbstractCurrentSource}, r::SVector{3, T}
    ) where {T}
    out_T = float(promote_type(T, _source_eltype(sources)))
    B = zero(SVector{3, out_T})
    for s in sources
        B += _solve_bs(OpenBoundary(), core_radius, s, r)
    end
    return B
end

function _solve_bs(
        bc::OpenBoundary, core_radius::Float64,
        ms::MultiSource{<:AbstractVector}, r::SVector{3}
    )
    return _solve_bs(bc, core_radius, ms.sources, r)
end

# ---------------------------------------------------------------------------
# Boundary conditions (method of images, prescribed background, periodic sum)
# ---------------------------------------------------------------------------

function _solve_bs(
        bc::ConductingWall, core_radius::Float64, source, r::SVector{3}
    )
    B_direct = _solve_bs(OpenBoundary(), core_radius, source, r)
    img = image_source(source, bc)
    B_image = _solve_bs(OpenBoundary(), core_radius, img, r)
    return B_direct + B_image
end

function _solve_bs(
        bc::PrescribedBoundary, core_radius::Float64, source, r::SVector{3}
    )
    B_source = _solve_bs(OpenBoundary(), core_radius, source, r)
    B_bg = bc.field_func(r)
    return B_source + SVector{3}(B_bg)
end

function _solve_bs(
        bc::PeriodicBoundary, core_radius::Float64, source, r::SVector{3}
    )
    Lx, Ly, Lz = bc.period
    nx, ny, nz = bc.n_images
    rx = Lx > 0.0 ? (-nx:nx) : (0:0)
    ry = Ly > 0.0 ? (-ny:ny) : (0:0)
    rz = Lz > 0.0 ? (-nz:nz) : (0:0)
    ET = typeof(r[1])
    B = _solve_bs(OpenBoundary(), core_radius, source, r) -
        _solve_bs(OpenBoundary(), core_radius, source, r)
    for kz in rz, ky in ry, kx in rx
        shift = SVector{3, ET}(kx * Lx, ky * Ly, kz * Lz)
        B += _solve_bs(OpenBoundary(), core_radius, source, r - shift)
    end
    return B
end

# ---------------------------------------------------------------------------
# Public solve dispatch
# ---------------------------------------------------------------------------

"""
    solve(solver::BiotSavart, source, r)

Compute the magnetic field `B` [T] at 3D point `r` (or batch of points) generated by
`source` (e.g. a `Wire`, `CurrentLoop`, `SurfaceCurrentMesh`, `MultiSource`, or
`AbstractVector{<:AbstractCurrentSource}`).
"""
function solve(solver::BiotSavart, source::_BiotSavartSource, r::SVector{3})
    return _solve_bs(solver.bc, solver.core_radius, source, r)
end

function solve(solver::BiotSavart, source::_BiotSavartSource, r)
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(_source_eltype(source), r)
    @inbounds r_sv = SVector{3, FT}(r[1], r[2], r[3])
    return _solve_bs(solver.bc, solver.core_radius, source, r_sv)
end

function solve(
        solver::BiotSavart, source::_BiotSavartSource,
        positions::AbstractVector{<:AbstractVector}
    )
    isempty(positions) && return SVector{3, _source_eltype(source)}[]
    first_val = solve(solver, source, positions[1])
    B_out = Vector{typeof(first_val)}(undef, length(positions))
    B_out[1] = first_val
    Threads.@threads for i in 2:length(positions)
        B_out[i] = solve(solver, source, positions[i])
    end
    return B_out
end

# ---------------------------------------------------------------------------
# ConductingSphere — numerical BEM
# ---------------------------------------------------------------------------

"""
    sphere_mesh(bc::ConductingSphere) -> (centers, normals, tangents1, tangents2, areas)

Build a latitude-longitude patch mesh for the conducting sphere.
"""
function sphere_mesh(bc::ConductingSphere{T}) where {T}
    nθ, nφ = bc.n_theta, bc.n_phi
    N = nθ * nφ
    centers = Vector{SVector{3, T}}(undef, N)
    normals = Vector{SVector{3, T}}(undef, N)
    tangents1 = Vector{SVector{3, T}}(undef, N)
    tangents2 = Vector{SVector{3, T}}(undef, N)
    areas = Vector{T}(undef, N)

    dθ = π / nθ
    dφ = 2π / nφ
    R = bc.radius

    idx = 1
    for iθ in 1:nθ
        θ = (iθ - 0.5) * dθ
        sinθ, cosθ = sincos(θ)
        for iφ in 1:nφ
            φ = (iφ - 0.5) * dφ
            sinφ, cosφ = sincos(φ)

            n̂ = SVector{3, T}(sinθ * cosφ, sinθ * sinφ, cosθ)
            centers[idx] = bc.center + R * n̂
            normals[idx] = n̂
            tangents1[idx] = SVector{3, T}(cosθ * cosφ, cosθ * sinφ, -sinθ) # θ̂
            tangents2[idx] = SVector{3, T}(-sinφ, cosφ, zero(T))              # φ̂
            areas[idx] = R^2 * sinθ * dθ * dφ
            idx += 1
        end
    end
    return centers, normals, tangents1, tangents2, areas
end

"""
    compute_surface_current(bc::ConductingSphere, source) -> SurfaceCurrentMesh

Solve the BEM linear system for the induced surface current **K** on the
perfectly-conducting sphere defined by `bc` that cancels the normal component
of the incident field from `source` (`Wire`, `CurrentLoop`, `MultiSource`, or
`AbstractVector{<:AbstractCurrentSource}`).
"""
function compute_surface_current(
        bc::ConductingSphere, source::_BiotSavartSource
    )
    centers, normals, tangents1, tangents2, areas = sphere_mesh(bc)
    N = length(centers)
    T = promote_type(eltype(bc.center), _source_eltype(source))

    bare = BiotSavart()

    b = Vector{T}(undef, N)
    for i in 1:N
        B_inc = solve(bare, source, centers[i])
        b[i] = -dot(normals[i], B_inc)
    end

    M = zeros(T, N, 2N)
    for j in 1:N
        Aj = T(areas[j])
        for i in 1:N
            i == j && continue
            R = centers[i] - centers[j]
            R_norm = norm(R)
            R_norm < 1.0e-20 && continue
            pre = T(μ0_4π) * Aj / R_norm^3
            for (k, t̂) in enumerate((tangents1[j], tangents2[j]))
                dB_n = dot(normals[i], pre * cross(t̂, R))
                M[i, 2j - 2 + k] = dB_n
            end
        end
    end

    x = M \ b

    K = [SVector{2, T}(x[2j - 1], x[2j]) for j in 1:N]
    return SurfaceCurrentMesh{T}(
        centers, normals, tangents1, tangents2,
        Vector{T}(areas), K
    )
end
