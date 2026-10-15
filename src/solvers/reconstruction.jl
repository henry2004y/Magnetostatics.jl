"""
    AxisymmetricFieldReconstruction <: AbstractMagneticField

Continuous axisymmetric magnetic field reconstructed from sparse measurements
using a compactly supported Wendland C2 radial-basis expansion of the poloidal
flux `ψ(R, Z)` by default, optionally augmented with global vacuum
(homogeneous Grad–Shafranov) polynomial solutions, virtual external filaments,
and softened toroidal current elements that carry net plasma current.
The poloidal field is divergence-free by construction. The toroidal field is
fit independently when toroidal measurements are supplied.
"""
struct AxisymmetricFieldReconstruction{T} <: AbstractMagneticField
    centers::Matrix{T}
    psi_coefficients::Vector{T}
    toroidal_coefficients::Vector{T}
    length_scale::T
    domain::NTuple{4, T}
    kernel::Symbol
    center_bins::Dict{Tuple{Int, Int}, Vector{Int}}
    vacuum_coefficients::Vector{T}
    vacuum_degree::Int
    vacuum_scale::T
    filament_positions::Matrix{T}
    filament_currents::Vector{T}
    element_positions::Matrix{T}
    element_currents::Vector{T}
    element_width::T
end

function AxisymmetricFieldReconstruction(
        centers::Matrix{T}, psi_coefficients::Vector{T}, toroidal_coefficients::Vector{T},
        length_scale::T, domain::NTuple{4, T}, kernel::Symbol,
        center_bins::Dict{Tuple{Int, Int}, Vector{Int}}
    ) where {T}
    return AxisymmetricFieldReconstruction(
        centers, psi_coefficients, toroidal_coefficients, length_scale, domain,
        kernel, center_bins, T[], -1, one(T), zeros(T, 0, 2), T[],
        zeros(T, 0, 2), T[], one(T)
    )
end

_n_vacuum_terms(degree::Int) = degree < 0 ? 0 : max(1, degree)

# The basis is traversed term by term: `f(k, ψ, ∂ψ/∂R, ∂ψ/∂Z, Δ*ψ)` is called for
# every column `k` that may be nonzero at `(R, Z)`, so that fitting and evaluation
# share one non-allocating implementation.

# Values and first derivatives of the constant and the degree `2:degree`
# polynomial solutions of `Δ*ψ = ∂²ψ/∂R² - (1/R)∂ψ/∂R + ∂²ψ/∂Z² = 0`, evaluated in
# coordinates normalized by `scale`. Each degree `n` contributes
# `Σ_j a_j R^(n-j) Z^j` with `a_(j+2) = -a_j m (m - 2) / ((j + 1)(j + 2))`,
# `m = n - j`, so even degrees are up-down symmetric and odd degrees antisymmetric.
function _foreach_vacuum_term(f::F, R, Z, offset, degree::Int, scale) where {F}
    degree < 0 && return nothing
    T = promote_type(typeof(R), typeof(Z), typeof(scale))
    x, z = R / scale, Z / scale
    f(offset + 1, one(T), zero(T), zero(T), zero(T))
    for n in 2:degree
        value, dR, dZ = zero(T), zero(T), zero(T)
        j = isodd(n) ? 1 : 0
        a = one(T)
        while n - j >= 2
            m = n - j
            value += a * x^m * z^j
            dR += a * m * x^(m - 1) * z^j / scale
            j > 0 && (dZ += a * j * x^m * z^(j - 1) / scale)
            a = -a * m * (m - 2) / ((j + 1) * (j + 2))
            j += 2
        end
        f(offset + n, value, dR, dZ, zero(T))
    end
    return nothing
end

function _vacuum_terms(R, Z, degree::Int, scale)
    T = promote_type(typeof(R), typeof(Z), typeof(scale))
    n_terms = _n_vacuum_terms(degree)
    values, dR, dZ = zeros(T, n_terms), zeros(T, n_terms), zeros(T, n_terms)
    _foreach_vacuum_term(R, Z, 0, degree, scale) do k, v, d_r, d_z, _
        values[k], dR[k], dZ[k] = v, d_r, d_z
    end
    return values, dR, dZ
end

# Flux `ψ`, `B_ρ`, `B_z` and `∂B_ρ/∂z` of a coaxial loop of radius `a` carrying
# `1/μ₀` amperes, at cylindrical radius `ρ` and height `z` above the loop plane.
function _unit_loop_terms(a, ρ, z)
    α2 = (a - ρ)^2 + z^2
    β2 = (a + ρ)^2 + z^2
    β = sqrt(β2)
    m = 4 * a * ρ / β2
    K, E = ellipk(m), ellipe(m)
    s2 = a^2 + ρ^2 + z^2
    psi = sqrt(a * ρ / m) * ((1 - m / 2) * K - E) / π
    B_z = (K + (a^2 - ρ^2 - z^2) / α2 * E) / (2π * β)
    P = z / (2π * ρ * β)
    G = s2 / α2
    B_ρ = P * (-K + G * E)
    dK = (E - (1 - m) * K) / (2m * (1 - m))
    dE = (E - K) / (2m)
    dm = -8 * a * ρ * z / β2^2
    dP = (1 / β - z^2 / (β * β2)) / (2π * ρ)
    dG = 2z / α2 - 2z * s2 / α2^2
    dB_ρ = dP * (-K + G * E) + P * (-dK * dm + dG * E + G * dE * dm)
    return psi, B_ρ, B_z, dB_ρ
end

# Poloidal flux and its derivatives for coaxial filaments carrying `1/μ₀` amperes,
# so that the fitted coefficients are `μ₀ I` and comparable in scale to the flux.
# Filaments lie outside the domain, so `Δ*ψ = 0`.
function _foreach_filament_term(f::F, R, Z, offset, filaments) where {F}
    T = promote_type(typeof(R), typeof(Z), eltype(filaments))
    for k in axes(filaments, 1)
        psi, B_ρ, B_z, _ = _unit_loop_terms(filaments[k, 1], R, Z - filaments[k, 2])
        f(offset + k, psi, R * B_z, -R * B_ρ, zero(T))
    end
    return nothing
end

# Flux, its first derivatives and `Δ*ψ` for softened toroidal current elements
# carrying `1/μ₀` amperes. Each element is the loop flux evaluated at the shifted
# height `Zₖ + √((Z - Zₖ)² + ε²)`, which is smooth everywhere, carries the full
# current and reduces to the loop field far from the element. Near the element,
# `J_ϕ ≈ (I/π) ε² / (d² + ε²)²`.
function _foreach_element_term(f::F, R, Z, offset, elements, width) where {F}
    ε2 = width^2
    for k in axes(elements, 1)
        ΔZ = Z - elements[k, 2]
        u = sqrt(ΔZ^2 + ε2)
        psi, B_ρ, B_z, dB_ρ = _unit_loop_terms(elements[k, 1], R, u)
        f(
            offset + k, psi, R * B_z, -R * B_ρ * ΔZ / u,
            ε2 * R / u^2 * (dB_ρ - B_ρ / u)
        )
    end
    return nothing
end

function _element_terms(R, Z, elements, width)
    T = promote_type(typeof(R), typeof(Z), eltype(elements), typeof(width))
    n = size(elements, 1)
    values, dR, dZ, delta = zeros(T, n), zeros(T, n), zeros(T, n), zeros(T, n)
    _foreach_element_term(R, Z, 0, elements, width) do k, v, d_r, d_z, gs
        values[k], dR[k], dZ[k], delta[k] = v, d_r, d_z, gs
    end
    return values, dR, dZ, delta
end

function _build_center_bins(centers, domain, length_scale)
    Rmin, _, Zmin, _ = domain
    bins = Dict{Tuple{Int, Int}, Vector{Int}}()
    for i in axes(centers, 1)
        cell = (
            floor(Int, (centers[i, 1] - Rmin) / length_scale),
            floor(Int, (centers[i, 2] - Zmin) / length_scale),
        )
        push!(get!(bins, cell, Int[]), i)
    end
    return bins
end

@inline function _rbf_kernel_terms(dr, dz, R, length_scale, inv_scale2, kernel::Symbol)
    T = promote_type(typeof(dr), typeof(dz), typeof(R), typeof(inv_scale2))
    distance = hypot(dr, dz)
    if kernel === :gaussian
        value = exp(-0.5 * distance^2 * inv_scale2)
        d_r = -dr * inv_scale2 * value
        d_z = -dz * inv_scale2 * value
        dRR = (dr^2 * inv_scale2^2 - inv_scale2) * value
        dZZ = (dz^2 * inv_scale2^2 - inv_scale2) * value
        return true, T(value), T(d_r), T(d_z), T(dRR + dZZ - d_r / R)
    elseif kernel === :wendland_c2
        q = distance / length_scale
        q >= 1 && return false, zero(T), zero(T), zero(T), zero(T)
        one_minus_q = 1 - q
        value = one_minus_q^4 * (4q + 1)
        phi_prime_over_q = -20 * one_minus_q^3
        d_r = phi_prime_over_q * dr * inv_scale2
        d_z = phi_prime_over_q * dz * inv_scale2
        radial_hessian = -20 * one_minus_q^2 * (1 - 4q) * inv_scale2
        tangential_hessian = phi_prime_over_q * inv_scale2
        delta = radial_hessian + tangential_hessian - d_r / R
        return true, T(value), T(d_r), T(d_z), T(delta)
    end
    throw(ArgumentError("kernel must be :wendland_c2 or :gaussian"))
end

@inline function _rbf_term(
        f::F, i, R, Z, centers, length_scale, inv_scale2, kernel
    ) where {F}
    supported, value, d_r, d_z, delta = _rbf_kernel_terms(
        R - centers[i, 1], Z - centers[i, 2], R, length_scale, inv_scale2, kernel
    )
    supported && f(i, value, d_r, d_z, delta)
    return nothing
end

function _foreach_rbf_term(
        f::F, R, Z, centers, length_scale, kernel,
        domain = nothing, center_bins = nothing
    ) where {F}
    inv_scale2 = inv(length_scale^2)
    if kernel === :wendland_c2 && center_bins !== nothing
        Rmin, _, Zmin, _ = domain
        cell_r = floor(Int, (R - Rmin) / length_scale)
        cell_z = floor(Int, (Z - Zmin) / length_scale)
        for dr in -1:1, dz in -1:1
            bin = get(center_bins, (cell_r + dr, cell_z + dz), nothing)
            bin === nothing && continue
            for i in bin
                _rbf_term(f, i, R, Z, centers, length_scale, inv_scale2, kernel)
            end
        end
    else
        for i in axes(centers, 1)
            _rbf_term(f, i, R, Z, centers, length_scale, inv_scale2, kernel)
        end
    end
    return nothing
end

function _rbf_terms(
        R, Z, centers, length_scale, kernel, domain = nothing, center_bins = nothing
    )
    buffer = _TermBuffer{
        promote_type(typeof(R), typeof(Z), eltype(centers), typeof(length_scale)),
    }()
    _foreach_rbf_term(buffer, R, Z, centers, length_scale, kernel, domain, center_bins)
    n = size(centers, 1)
    sparse_terms(values) = sparsevec(buffer.indices, values, n)
    return sparse_terms(buffer.values), sparse_terms(buffer.dR),
        sparse_terms(buffer.dZ), sparse_terms(buffer.delta)
end

# `basis` is either an `AxisymmetricFieldReconstruction` or a named tuple with the
# same basis property names. Columns are ordered as RBF centers, vacuum terms,
# filaments and current elements. With `flux = false` only `Δ*ψ` is required, so
# the current-free vacuum and filament terms are skipped.
function _foreach_basis_term(f::F, basis, R, Z; flux::Bool = true) where {F}
    centers = basis.centers
    _foreach_rbf_term(
        f, R, Z, centers, basis.length_scale, basis.kernel, basis.domain,
        basis.center_bins
    )
    offset = size(centers, 1)
    if flux
        _foreach_vacuum_term(f, R, Z, offset, basis.vacuum_degree, basis.vacuum_scale)
    end
    offset += _n_vacuum_terms(basis.vacuum_degree)
    flux && _foreach_filament_term(f, R, Z, offset, basis.filament_positions)
    offset += size(basis.filament_positions, 1)
    _foreach_element_term(
        f, R, Z, offset, basis.element_positions, basis.element_width
    )
    return nothing
end

struct _TermBuffer{T}
    indices::Vector{Int}
    values::Vector{T}
    dR::Vector{T}
    dZ::Vector{T}
    delta::Vector{T}
end

_TermBuffer{T}() where {T} = _TermBuffer{T}(Int[], T[], T[], T[], T[])

function (buffer::_TermBuffer)(k, value, dR, dZ, delta)
    push!(buffer.indices, k)
    push!(buffer.values, value)
    push!(buffer.dR, dR)
    push!(buffer.dZ, dZ)
    push!(buffer.delta, delta)
    return nothing
end

function _collect_basis_terms!(buffer::_TermBuffer, basis, R, Z; flux::Bool = true)
    for v in (buffer.indices, buffer.values, buffer.dR, buffer.dZ, buffer.delta)
        empty!(v)
    end
    _foreach_basis_term(buffer, basis, R, Z; flux)
    return buffer
end

@inline function _psi_coefficient(field::AxisymmetricFieldReconstruction, k)
    n = length(field.psi_coefficients)
    k <= n && return field.psi_coefficients[k]
    k -= n
    n = length(field.vacuum_coefficients)
    k <= n && return field.vacuum_coefficients[k]
    k -= n
    n = length(field.filament_currents)
    k <= n && return μ₀ * field.filament_currents[k]
    return μ₀ * field.element_currents[k - n]
end

# Sums of `(ψ, ∂ψ/∂R, ∂ψ/∂Z, Δ*ψ)` with the poloidal coefficients and of
# `(F, ∂F/∂R, ∂F/∂Z)` with the toroidal RBF coefficients, where `B_ϕ = F`.
function _basis_sums(field::AxisymmetricFieldReconstruction, R, Z; flux::Bool = true)
    T = promote_type(eltype(field.psi_coefficients), typeof(R), typeof(Z))
    toroidal = field.toroidal_coefficients
    ncenters = length(toroidal)
    sums = Ref(zero(SVector{7, T}))
    _foreach_basis_term(field, R, Z; flux) do k, value, dR, dZ, delta
        c = _psi_coefficient(field, k)
        t = k <= ncenters ? toroidal[k] : zero(eltype(toroidal))
        sums[] += SVector(c * value, c * dR, c * dZ, c * delta, t * value, t * dR, t * dZ)
    end
    return sums[]
end

# Weighted least-squares rows `A[i, J] = V / σᵢ`, `b[i] = targetᵢ / σᵢ` in
# coordinate format.
struct _FitRows{T}
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{T}
    b::Vector{T}
end

_FitRows{T}() where {T} = _FitRows{T}(Int[], Int[], T[], T[])

function _push_row!(rows::_FitRows, indices, terms, factor, target, ncolumns)
    row = length(rows.b) + 1
    for (k, term) in zip(indices, terms)
        k <= ncolumns || continue
        value = factor * term
        iszero(value) && continue
        push!(rows.I, row)
        push!(rows.J, k)
        push!(rows.V, value)
    end
    push!(rows.b, target)
    return nothing
end

function _append_weighted_row!(
        rows::_FitRows, indices, terms, factor, target, uncertainty,
        ncolumns = typemax(Int)
    )
    isfinite(target) || return nothing
    isfinite(uncertainty) && uncertainty > 0 ||
        throw(ArgumentError("uncertainties must be finite and positive"))
    weight = inv(uncertainty)
    _push_row!(rows, indices, terms, factor * weight, target * weight, ncolumns)
    return nothing
end

function _append_fit_row!(rows, targets, uncertainties, row, target, uncertainty)
    isfinite(target) || return nothing
    isfinite(uncertainty) && uncertainty > 0 ||
        throw(ArgumentError("uncertainties must be finite and positive"))
    push!(rows, sparsevec(row))
    push!(targets, target)
    push!(uncertainties, uncertainty)
    return nothing
end

# Sparse QR is slower and uses more memory than dense QR once the design matrix is
# mostly filled, e.g. for the global Gaussian kernel or many current elements.
function _qr_design_matrix(A::SparseMatrixCSC, algorithm)
    algorithm isa QRFactorization && 3 * nnz(A) > length(A) || return A
    return Matrix(A)
end

function _fit_rbf(
        rows::_FitRows{T}, regularization, ncoefficients, exact::_FitRows{T}, algorithm
    ) where {T}
    nrows = length(rows.b)
    regularization_weight = sqrt(regularization)
    I = append!(copy(rows.I), (nrows + 1):(nrows + ncoefficients))
    J = append!(copy(rows.J), 1:ncoefficients)
    V = append!(copy(rows.V), Iterators.repeated(regularization_weight, ncoefficients))
    A = sparse(I, J, V, nrows + ncoefficients, ncoefficients)
    A = _qr_design_matrix(A, algorithm)
    b = vcat(rows.b, zeros(T, ncoefficients))

    isempty(exact.b) &&
        return LinearSolve.solve(LinearProblem(A, b), algorithm).u

    C = zeros(T, length(exact.b), ncoefficients)
    for (i, j, v) in zip(exact.I, exact.J, exact.V)
        C[i, j] = v
    end
    d = exact.b
    particular = pinv(C) * d
    residual = norm(C * particular - d, Inf)
    tolerance = 100 * eps(T) * max(norm(d, Inf), one(T))
    residual <= tolerance ||
        throw(ArgumentError("exact boundary flux constraints are inconsistent"))

    nullspace_basis = nullspace(C)
    size(nullspace_basis, 2) == 0 && return particular
    reduced = A * nullspace_basis
    correction = LinearSolve.solve(
        LinearProblem(reduced, b - A * particular), algorithm
    ).u
    return particular + nullspace_basis * correction
end

function _validate_rz_positions(positions, domain, name)
    size(positions, 2) == 2 ||
        throw(ArgumentError("$name must have columns (R, Z)"))
    all(isfinite, positions) || throw(ArgumentError("$name must be finite"))
    Rmin, Rmax, Zmin, Zmax = domain
    all(Rmin .<= positions[:, 1] .<= Rmax) ||
        throw(ArgumentError("$name radii must lie inside domain"))
    all(Zmin .<= positions[:, 2] .<= Zmax) ||
        throw(ArgumentError("$name heights must lie inside domain"))
    return nothing
end

function _median_spacing(points, fallback)
    n = size(points, 1)
    n < 2 && return fallback
    nearest = map(1:n) do i
        minimum(
            hypot(points[i, 1] - points[j, 1], points[i, 2] - points[j, 2])
                for j in 1:n if j != i
        )
    end
    spacing = sort(nearest)[cld(n, 2)]
    return spacing > 0 ? spacing : fallback
end

function _default_rbf_centers(domain, nsamples, ::Type{T}) where {T}
    Rmin, Rmax, Zmin, Zmax = domain
    nside = max(2, ceil(Int, sqrt(nsamples)))
    Rs = range(Rmin, Rmax; length = nside)
    Zs = range(Zmin, Zmax; length = nside)
    return reduce(vcat, [reshape(T[R, Z], 1, 2) for R in Rs for Z in Zs])
end

"""
    reconstruct_axisymmetric_field(positions, B_observed; domain, length_scale,
                                   regularization, uncertainties, algorithm, ...)

Fit a continuous axisymmetric field to sparse magnetic-field samples. `positions`
is an `n × 2` matrix of `(R, Z)` coordinates and `B_observed` is an `n × 3`
matrix of cylindrical `(B_R, B_ϕ, B_Z)` components. Non-finite field components
are treated as missing measurements. The default `:wendland_c2` kernel has
support radius `length_scale` and produces a sparse design matrix; `:gaussian`
selects a global Gaussian basis. The default iterative algorithm is LSMR and can
be customized with `algorithm`.

`vacuum_degree = N` augments the local basis with the constant and the
degree `2:N` polynomial solutions of `Δ*ψ = 0` (current-free fields such as a
uniform `B_Z`, quadrupole, and hexapole components, normalized by `Rmax`).
These global terms capture the smooth field from distant coils that a compact
basis alone represents poorly; the radial basis then models local currents.

`filament_positions` is an `m × 2` matrix of `(R, Z)` locations, outside the
domain, of virtual coaxial current filaments whose currents are fit
(method of fundamental solutions). Each filament field is current-free inside the
domain, so they efficiently represent fields from nearby external coils; place
them between the domain and the real sources, e.g. on a curve surrounding the
domain. Fitted currents are stored in `filament_currents` [A].

`current_element_positions` is an `m × 2` matrix of `(R, Z)` locations of
softened toroidal current elements whose currents are fit and stored in
`element_currents` [A]. Each element is a coaxial loop with a smooth core of
width `current_element_width` (default: median nearest-neighbor spacing of the
elements), with current density `J_ϕ ≈ (I/π) ε² / (d² + ε²)²` and the vacuum loop
field far away. Unlike radial basis functions, elements carry net toroidal current,
so `sum(field.element_currents)` estimates the plasma current. Place them on a
grid covering the region where plasma current may flow, with spacing comparable
to the smallest current structure.

`boundary_positions` and `boundary_flux` impose exact poloidal-flux constraints.
Supply `boundary_uncertainties` to treat them as weighted observations instead.
`current_positions` and `current_density` add toroidal-current samples
through the Grad–Shafranov equation, weighted by `current_uncertainties` in the
same units as `current_density` (A/m²). Setting `current_density` to zero at
collocation points imposes a vacuum prior; a loose uncertainty (comparable to the
expected plasma current density) instead acts as a minimum-current regularization.
Compactly supported `ψ` bumps carry no net toroidal current, so for fields with
a distributed plasma current add current elements (or use the `:gaussian`
kernel, whose wide tails can place return currents outside the domain). The fitted
current density is available through `current_density(field, pos)`. The domain is
`(Rmin, Rmax, Zmin, Zmax)`;
evaluation outside it throws `DomainError`. Singular sources such as coil
filaments cannot be represented and should lie outside the domain. The default
current uncertainty is `1` A/m². Exact boundary constraints currently
use a dense null-space reduction; pass boundary uncertainties for the scalable
iterative path.
"""
function reconstruct_axisymmetric_field(
        positions::AbstractMatrix{<:Real}, B_observed::AbstractMatrix{<:Real};
        domain, length_scale, regularization = 1.0e-8,
        uncertainties = ones(size(B_observed)),
        algorithm = KrylovJL_LSMR(),
        kernel = :wendland_c2,
        basis_centers = nothing,
        boundary_positions = nothing, boundary_flux = nothing,
        boundary_uncertainties = nothing,
        current_positions = nothing, current_density = nothing,
        current_uncertainties = nothing,
        vacuum_degree::Union{Nothing, Integer} = nothing,
        filament_positions = nothing,
        current_element_positions = nothing,
        current_element_width = nothing
    )
    size(positions, 2) == 2 ||
        throw(ArgumentError("positions must have columns (R, Z)"))
    size(B_observed) == (size(positions, 1), 3) ||
        throw(ArgumentError("B_observed must have columns (B_R, B_ϕ, B_Z)"))
    size(uncertainties) == size(B_observed) ||
        throw(ArgumentError("uncertainties must have the same size as B_observed"))
    length(domain) == 4 ||
        throw(ArgumentError("domain must be (Rmin, Rmax, Zmin, Zmax)"))

    T = promote_type(Float64, eltype(positions), eltype(B_observed), typeof(length_scale))
    domain_t = T.(domain)
    Rmin, Rmax, Zmin, Zmax = domain_t
    all(isfinite, domain_t) && Rmin > 0 && Rmax > Rmin && Zmax > Zmin ||
        throw(ArgumentError("domain must have Rmin > 0, Rmax > Rmin, Zmax > Zmin"))
    isfinite(length_scale) && length_scale > 0 ||
        throw(ArgumentError("length_scale must be finite and positive"))
    isfinite(regularization) && regularization > 0 ||
        throw(ArgumentError("regularization must be finite and positive"))
    kernel in (:wendland_c2, :gaussian) ||
        throw(ArgumentError("kernel must be :wendland_c2 or :gaussian"))

    sample_positions = Matrix{T}(positions)
    observed = Matrix{T}(B_observed)
    uncertainty_values = Matrix{T}(uncertainties)
    domain_typed = (Rmin, Rmax, Zmin, Zmax)
    _validate_rz_positions(sample_positions, domain_typed, "positions")

    centers = basis_centers === nothing ?
        _default_rbf_centers(domain_typed, size(sample_positions, 1), T) :
        Matrix{T}(basis_centers)
    _validate_rz_positions(centers, domain_typed, "basis_centers")
    size(centers, 1) > 0 || throw(ArgumentError("at least one basis center is required"))
    center_bins = _build_center_bins(centers, domain_typed, T(length_scale))
    vacuum_degree_int = vacuum_degree === nothing ? -1 : Int(vacuum_degree)
    vacuum_degree_int >= 0 || vacuum_degree === nothing ||
        throw(ArgumentError("vacuum_degree must be nonnegative"))
    vacuum_scale = Rmax
    filaments = filament_positions === nothing ? zeros(T, 0, 2) :
        Matrix{T}(filament_positions)
    size(filaments, 2) == 2 ||
        throw(ArgumentError("filament_positions must have columns (R, Z)"))
    all(isfinite, filaments) && all(>(0), view(filaments, :, 1)) ||
        throw(ArgumentError("filament_positions must be finite with R > 0"))
    any(
        Rmin <= filaments[k, 1] <= Rmax && Zmin <= filaments[k, 2] <= Zmax
            for k in axes(filaments, 1)
    ) && throw(ArgumentError("filament_positions must lie outside the domain"))
    elements = current_element_positions === nothing ? zeros(T, 0, 2) :
        Matrix{T}(current_element_positions)
    size(elements, 2) == 2 ||
        throw(ArgumentError("current_element_positions must have columns (R, Z)"))
    all(isfinite, elements) && all(>(0), view(elements, :, 1)) ||
        throw(ArgumentError("current_element_positions must be finite with R > 0"))
    width = current_element_width === nothing ?
        _median_spacing(elements, T(length_scale)) : current_element_width
    isfinite(width) && width > 0 ||
        throw(ArgumentError("current_element_width must be finite and positive"))
    element_width = T(width)
    ncenters = size(centers, 1)
    n_vacuum = _n_vacuum_terms(vacuum_degree_int)
    n_filaments = size(filaments, 1)
    basis = (;
        centers, length_scale = T(length_scale), kernel, domain = domain_typed,
        center_bins, vacuum_degree = vacuum_degree_int, vacuum_scale,
        filament_positions = filaments, element_positions = elements, element_width,
    )
    terms = _TermBuffer{T}()
    psi_rows = _FitRows{T}()
    exact_boundary_rows = _FitRows{T}()
    toroidal_rows = _FitRows{T}()

    for i in axes(sample_positions, 1)
        R, Z = sample_positions[i, 1], sample_positions[i, 2]
        _collect_basis_terms!(terms, basis, R, Z)
        for component in 1:3
            target = observed[i, component]
            isfinite(target) || continue
            sigma = uncertainty_values[i, component]
            if component == 1
                _append_weighted_row!(
                    psi_rows, terms.indices, terms.dZ, -inv(R), target, sigma
                )
            elseif component == 2
                _append_weighted_row!(
                    toroidal_rows, terms.indices, terms.values, one(T), target, sigma,
                    ncenters
                )
            else
                _append_weighted_row!(
                    psi_rows, terms.indices, terms.dR, inv(R), target, sigma
                )
            end
        end
    end

    if boundary_positions !== nothing || boundary_flux !== nothing
        boundary_positions !== nothing && boundary_flux !== nothing ||
            throw(
            ArgumentError(
                "boundary_positions and boundary_flux must be given together"
            )
        )
        size(boundary_positions, 1) == length(boundary_flux) ||
            throw(ArgumentError("boundary_positions and boundary_flux lengths differ"))
        _validate_rz_positions(boundary_positions, domain_typed, "boundary_positions")
        sigma = boundary_uncertainties === nothing ? nothing : T.(boundary_uncertainties)
        sigma === nothing || length(sigma) == length(boundary_flux) ||
            throw(ArgumentError("boundary_uncertainties length differs"))
        for i in eachindex(boundary_flux)
            _collect_basis_terms!(
                terms, basis, T(boundary_positions[i, 1]), T(boundary_positions[i, 2])
            )
            if sigma === nothing
                isfinite(boundary_flux[i]) ||
                    throw(ArgumentError("boundary_flux must be finite"))
                _push_row!(
                    exact_boundary_rows, terms.indices, terms.values, one(T),
                    T(boundary_flux[i]), typemax(Int)
                )
            else
                _append_weighted_row!(
                    psi_rows, terms.indices, terms.values, one(T),
                    T(boundary_flux[i]), sigma[i]
                )
            end
        end
    end

    if current_positions !== nothing || current_density !== nothing
        current_positions !== nothing && current_density !== nothing ||
            throw(
            ArgumentError(
                "current_positions and current_density must be given together"
            )
        )
        size(current_positions, 1) == length(current_density) ||
            throw(ArgumentError("current_positions and current_density lengths differ"))
        _validate_rz_positions(current_positions, domain_typed, "current_positions")
        sigma = current_uncertainties === nothing ?
            ones(T, length(current_density)) : T.(current_uncertainties)
        length(sigma) == length(current_density) ||
            throw(ArgumentError("current_uncertainties length differs"))
        for i in eachindex(current_density)
            R, Z = T(current_positions[i, 1]), T(current_positions[i, 2])
            _collect_basis_terms!(terms, basis, R, Z; flux = false)
            target = -T(μ₀) * R * T(current_density[i])
            _append_weighted_row!(
                psi_rows, terms.indices, terms.delta, one(T), target,
                T(μ₀) * R * sigma[i]
            )
        end
    end

    isempty(psi_rows.b) && isempty(exact_boundary_rows.b) &&
        throw(ArgumentError("poloidal or boundary constraints are required"))
    all_psi_coefficients = _fit_rbf(
        psi_rows, T(regularization),
        ncenters + n_vacuum + n_filaments + size(elements, 1),
        exact_boundary_rows, algorithm
    )
    psi_coefficients = all_psi_coefficients[1:ncenters]
    vacuum_coefficients = all_psi_coefficients[(ncenters + 1):(ncenters + n_vacuum)]
    n_external = ncenters + n_vacuum + n_filaments
    filament_currents =
        all_psi_coefficients[(ncenters + n_vacuum + 1):n_external] ./ T(μ₀)
    element_currents = all_psi_coefficients[(n_external + 1):end] ./ T(μ₀)
    toroidal_coefficients = isempty(toroidal_rows.b) ? zeros(T, ncenters) :
        _fit_rbf(toroidal_rows, T(regularization), ncenters, _FitRows{T}(), algorithm)

    return AxisymmetricFieldReconstruction(
        centers, psi_coefficients, toroidal_coefficients,
        T(length_scale), domain_typed, kernel, center_bins,
        vacuum_coefficients, vacuum_degree_int, vacuum_scale,
        filaments, filament_currents, elements, element_currents, element_width
    )
end

function poloidal_flux(field::AxisymmetricFieldReconstruction, pos)
    R = length(pos) == 2 ? pos[1] : hypot(pos[1], pos[2])
    Z = length(pos) == 2 ? pos[2] : pos[3]
    Rmin, Rmax, Zmin, Zmax = field.domain
    Rmin <= R <= Rmax && Zmin <= Z <= Zmax ||
        throw(DomainError(pos, "position is outside the reconstruction domain"))
    return _basis_sums(field, R, Z)[1]
end

function _cylindrical_field(field::AxisymmetricFieldReconstruction, R, Z)
    sums = _basis_sums(field, R, Z)
    return -sums[3] / R, sums[5], sums[2] / R
end

function (field::AxisymmetricFieldReconstruction)(pos)
    x, y, z = pos[1], pos[2], pos[3]
    R = hypot(x, y)
    Rmin, Rmax, Zmin, Zmax = field.domain
    Rmin <= R <= Rmax && Zmin <= z <= Zmax ||
        throw(DomainError(pos, "position is outside the reconstruction domain"))

    B_R, B_phi, B_Z = _cylindrical_field(field, R, z)
    cos_phi, sin_phi = x / R, y / R
    return SVector(
        B_R * cos_phi - B_phi * sin_phi,
        B_R * sin_phi + B_phi * cos_phi,
        B_Z
    )
end

function current_density(field::AxisymmetricFieldReconstruction, pos)
    x, y, z = pos[1], pos[2], pos[3]
    R = hypot(x, y)
    Rmin, Rmax, Zmin, Zmax = field.domain
    Rmin <= R <= Rmax && Zmin <= z <= Zmax ||
        throw(DomainError(pos, "position is outside the reconstruction domain"))

    sums = _basis_sums(field, R, z; flux = false)
    B_phi = sums[5]
    J_R = -sums[7] / μ₀
    J_phi = -sums[4] / (μ₀ * R)
    J_Z = (sums[6] + B_phi / R) / μ₀
    cos_phi, sin_phi = x / R, y / R
    return SVector(
        J_R * cos_phi - J_phi * sin_phi,
        J_R * sin_phi + J_phi * cos_phi,
        J_Z
    )
end
