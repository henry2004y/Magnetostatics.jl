"""
    AxisymmetricFieldReconstruction <: AbstractMagneticField

Continuous axisymmetric magnetic field reconstructed from sparse measurements
using a compactly supported Wendland C2 radial-basis expansion of the poloidal
flux `ψ(R, Z)` by default, optionally augmented with global vacuum
(homogeneous Grad–Shafranov) polynomial solutions.
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
end

function AxisymmetricFieldReconstruction(
        centers::Matrix{T}, psi_coefficients::Vector{T}, toroidal_coefficients::Vector{T},
        length_scale::T, domain::NTuple{4, T}, kernel::Symbol,
        center_bins::Dict{Tuple{Int, Int}, Vector{Int}}
    ) where {T}
    return AxisymmetricFieldReconstruction(
        centers, psi_coefficients, toroidal_coefficients, length_scale, domain,
        kernel, center_bins, T[], -1, one(T), zeros(T, 0, 2), T[]
    )
end

_n_vacuum_terms(degree::Int) = degree < 0 ? 0 : max(1, degree)

# Values and first derivatives of the constant and the degree `2:degree`
# polynomial solutions of `Δ*ψ = ∂²ψ/∂R² - (1/R)∂ψ/∂R + ∂²ψ/∂Z² = 0`, evaluated in
# coordinates normalized by `scale`. Each degree `n` contributes
# `Σ_j a_j R^(n-j) Z^j` with `a_(j+2) = -a_j m (m - 2) / ((j + 1)(j + 2))`,
# `m = n - j`, so even degrees are up-down symmetric and odd degrees antisymmetric.
function _vacuum_terms(R, Z, degree::Int, scale)
    T = promote_type(typeof(R), typeof(Z), typeof(scale))
    n_terms = _n_vacuum_terms(degree)
    values = zeros(T, n_terms)
    dR = zeros(T, n_terms)
    dZ = zeros(T, n_terms)
    n_terms == 0 && return values, dR, dZ
    x, z = R / scale, Z / scale
    values[1] = one(T)
    for n in 2:degree
        k = n
        j = isodd(n) ? 1 : 0
        a = one(T)
        while n - j >= 2
            m = n - j
            values[k] += a * x^m * z^j
            dR[k] += a * m * x^(m - 1) * z^j / scale
            j > 0 && (dZ[k] += a * j * x^m * z^(j - 1) / scale)
            a = -a * m * (m - 2) / ((j + 1) * (j + 2))
            j += 2
        end
    end
    return values, dR, dZ
end

# Poloidal flux and its derivatives for coaxial filaments carrying `1/μ₀` amperes,
# so that the fitted coefficients are `μ₀ I` and comparable in scale to the flux.
function _filament_terms(R, Z, filaments)
    T = promote_type(typeof(R), typeof(Z), eltype(filaments))
    n = size(filaments, 1)
    values, dR, dZ = zeros(T, n), zeros(T, n), zeros(T, n)
    pos = SVector(R, zero(R), Z)
    for k in 1:n
        loop = CurrentLoop(
            filaments[k, 1], inv(T(μ₀)), SVector(zero(T), zero(T), T(filaments[k, 2])),
            SVector(zero(T), zero(T), one(T))
        )
        B = getB_loop(pos, loop)
        values[k] = R * solve(VectorPotential(), loop, pos)[2]
        dR[k] = R * B[3]
        dZ[k] = -R * B[1]
    end
    return values, dR, dZ
end

function _basis_terms(
        R, Z, centers, length_scale, kernel, domain, center_bins,
        vacuum_degree, vacuum_scale, filaments
    )
    values, dR, dZ, delta = _rbf_terms(
        R, Z, centers, length_scale, kernel, domain, center_bins
    )
    vacuum_degree < 0 && size(filaments, 1) == 0 && return values, dR, dZ, delta
    v_values, v_dR, v_dZ = _vacuum_terms(R, Z, vacuum_degree, vacuum_scale)
    f_values, f_dR, f_dZ = _filament_terms(R, Z, filaments)
    n_global = length(v_values) + length(f_values)
    return vcat(values, sparse(vcat(v_values, f_values))),
        vcat(dR, sparse(vcat(v_dR, f_dR))), vcat(dZ, sparse(vcat(v_dZ, f_dZ))),
        vcat(delta, spzeros(eltype(delta), n_global))
end

function _basis_terms(field::AxisymmetricFieldReconstruction, R, Z)
    return _basis_terms(
        R, Z, field.centers, field.length_scale, field.kernel, field.domain,
        field.center_bins, field.vacuum_degree, field.vacuum_scale,
        field.filament_positions
    )
end

function _psi_coefficient_vector(field::AxisymmetricFieldReconstruction)
    return vcat(
        field.psi_coefficients, field.vacuum_coefficients,
        μ₀ .* field.filament_currents
    )
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

function _center_candidates(R, Z, centers, length_scale, domain, center_bins)
    center_bins === nothing && return axes(centers, 1)
    Rmin, _, Zmin, _ = domain
    cell_r = floor(Int, (R - Rmin) / length_scale)
    cell_z = floor(Int, (Z - Zmin) / length_scale)
    candidates = Int[]
    for dr in -1:1, dz in -1:1
        append!(candidates, get(center_bins, (cell_r + dr, cell_z + dz), Int[]))
    end
    return candidates
end

function _rbf_terms(
        R, Z, centers, length_scale, kernel, domain = nothing, center_bins = nothing
    )
    n = size(centers, 1)
    T = promote_type(typeof(R), typeof(Z), eltype(centers), typeof(length_scale))
    value_indices, value_values = Int[], T[]
    dR_indices, dR_values = Int[], T[]
    dZ_indices, dZ_values = Int[], T[]
    delta_indices, delta_values = Int[], T[]
    inv_scale2 = inv(length_scale^2)

    candidate_centers = if kernel === :wendland_c2
        _center_candidates(R, Z, centers, length_scale, domain, center_bins)
    else
        axes(centers, 1)
    end

    @inbounds for i in candidate_centers
        dr = R - centers[i, 1]
        dz = Z - centers[i, 2]
        distance = hypot(dr, dz)
        if kernel === :gaussian
            value = exp(-0.5 * distance^2 * inv_scale2)
            d_r = -dr * inv_scale2 * value
            d_z = -dz * inv_scale2 * value
            dRR = (dr^2 * inv_scale2^2 - inv_scale2) * value
            dZZ = (dz^2 * inv_scale2^2 - inv_scale2) * value
            delta = dRR + dZZ - d_r / R
        elseif kernel === :wendland_c2
            q = distance / length_scale
            if q >= 1
                continue
            else
                one_minus_q = 1 - q
                value = one_minus_q^4 * (4q + 1)
                phi_prime_over_q = -20 * one_minus_q^3
                d_r = phi_prime_over_q * dr * inv_scale2
                d_z = phi_prime_over_q * dz * inv_scale2
                radial_hessian = -20 * one_minus_q^2 * (1 - 4q) * inv_scale2
                tangential_hessian = phi_prime_over_q * inv_scale2
                delta = radial_hessian + tangential_hessian - d_r / R
            end
        else
            throw(ArgumentError("kernel must be :wendland_c2 or :gaussian"))
        end
        if !iszero(value)
            push!(value_indices, i)
            push!(value_values, value)
        end
        if !iszero(d_r)
            push!(dR_indices, i)
            push!(dR_values, d_r)
        end
        if !iszero(d_z)
            push!(dZ_indices, i)
            push!(dZ_values, d_z)
        end
        if !iszero(delta)
            push!(delta_indices, i)
            push!(delta_values, delta)
        end
    end

    return sparsevec(value_indices, value_values, n),
        sparsevec(dR_indices, dR_values, n), sparsevec(dZ_indices, dZ_values, n),
        sparsevec(delta_indices, delta_values, n)
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

function _fit_rbf(
        rows, targets, uncertainties, regularization, ncenters,
        exact_rows, exact_targets, algorithm
    )
    nrows = length(rows)
    T = promote_type(eltype(targets), eltype(uncertainties))
    row_indices = Int[]
    column_indices = Int[]
    matrix_values = T[]
    b = zeros(T, nrows + ncenters)

    for i in 1:nrows
        inv_uncertainty = inv(uncertainties[i])
        indices, values = findnz(rows[i])
        for j in eachindex(indices)
            push!(row_indices, i)
            push!(column_indices, indices[j])
            push!(matrix_values, values[j] * inv_uncertainty)
        end
        b[i] = targets[i] * inv_uncertainty
    end
    for i in 1:ncenters
        push!(row_indices, nrows + i)
        push!(column_indices, i)
        push!(matrix_values, sqrt(regularization))
    end
    A = sparse(row_indices, column_indices, matrix_values, nrows + ncenters, ncenters)

    isempty(exact_rows) &&
        return LinearSolve.solve(LinearProblem(A, b), algorithm).u

    C = zeros(T, length(exact_rows), ncenters)
    for i in eachindex(exact_rows)
        C[i, :] .= exact_rows[i]
    end
    d = exact_targets
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

`boundary_positions` and `boundary_flux` impose exact poloidal-flux constraints.
Supply `boundary_uncertainties` to treat them as weighted observations instead.
`current_positions` and `current_density` add toroidal-current samples
through the Grad–Shafranov equation, weighted by `current_uncertainties` in the
same units as `current_density` (A/m²). Setting `current_density` to zero at
collocation points imposes a vacuum prior; a loose uncertainty (comparable to the
expected plasma current density) instead acts as a minimum-current regularization.
Compactly supported `ψ` bumps carry no net toroidal current, so for fields with
a distributed plasma current prefer the `:gaussian` kernel. The fitted current
density is available through `current_density(field, pos)`. The domain is
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
        filament_positions = nothing
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
    ncenters = size(centers, 1)
    n_vacuum = _n_vacuum_terms(vacuum_degree_int)
    basis_terms(R, Z) = _basis_terms(
        R, Z, centers, T(length_scale), kernel, domain_typed, center_bins,
        vacuum_degree_int, vacuum_scale, filaments
    )
    psi_rows = SparseVector{T, Int}[]
    psi_targets = T[]
    psi_uncertainties = T[]
    exact_boundary_rows = SparseVector{T, Int}[]
    exact_boundary_targets = T[]
    toroidal_rows = SparseVector{T, Int}[]
    toroidal_targets = T[]
    toroidal_uncertainty_values = T[]

    for i in axes(sample_positions, 1)
        R, Z = sample_positions[i, 1], sample_positions[i, 2]
        values, dR, dZ, _ = basis_terms(R, Z)
        for component in 1:3
            target = observed[i, component]
            isfinite(target) || continue
            sigma = uncertainty_values[i, component]
            if component == 1
                row = -dZ / R
                _append_fit_row!(
                    psi_rows, psi_targets, psi_uncertainties,
                    row, target, sigma
                )
            elseif component == 2
                _append_fit_row!(
                    toroidal_rows, toroidal_targets,
                    toroidal_uncertainty_values, values[1:ncenters], target, sigma
                )
            else
                row = dR / R
                _append_fit_row!(
                    psi_rows, psi_targets, psi_uncertainties,
                    row, target, sigma
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
            values = first(basis_terms(boundary_positions[i, 1], boundary_positions[i, 2]))
            if sigma === nothing
                isfinite(boundary_flux[i]) ||
                    throw(ArgumentError("boundary_flux must be finite"))
                push!(exact_boundary_rows, values)
                push!(exact_boundary_targets, T(boundary_flux[i]))
            else
                _append_fit_row!(
                    psi_rows, psi_targets, psi_uncertainties,
                    values, T(boundary_flux[i]), sigma[i]
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
            R, Z = current_positions[i, 1], current_positions[i, 2]
            delta_star = last(basis_terms(R, Z))
            target = -T(μ₀) * R * T(current_density[i])
            _append_fit_row!(
                psi_rows, psi_targets, psi_uncertainties,
                delta_star, target, T(μ₀) * R * sigma[i]
            )
        end
    end

    isempty(psi_rows) && isempty(exact_boundary_rows) &&
        throw(ArgumentError("poloidal or boundary constraints are required"))
    all_psi_coefficients = _fit_rbf(
        psi_rows, psi_targets, psi_uncertainties, T(regularization),
        ncenters + n_vacuum + size(filaments, 1),
        exact_boundary_rows, exact_boundary_targets, algorithm
    )
    psi_coefficients = all_psi_coefficients[1:ncenters]
    vacuum_coefficients = all_psi_coefficients[(ncenters + 1):(ncenters + n_vacuum)]
    filament_currents = all_psi_coefficients[(ncenters + n_vacuum + 1):end] ./ T(μ₀)
    toroidal_coefficients = isempty(toroidal_rows) ? zeros(T, ncenters) :
        _fit_rbf(
            toroidal_rows, toroidal_targets, toroidal_uncertainty_values,
            T(regularization), ncenters, SparseVector{T, Int}[], T[], algorithm
        )

    return AxisymmetricFieldReconstruction(
        centers, psi_coefficients, toroidal_coefficients,
        T(length_scale), domain_typed, kernel, center_bins,
        vacuum_coefficients, vacuum_degree_int, vacuum_scale,
        filaments, filament_currents
    )
end

function poloidal_flux(field::AxisymmetricFieldReconstruction, pos)
    R = length(pos) == 2 ? pos[1] : hypot(pos[1], pos[2])
    Z = length(pos) == 2 ? pos[2] : pos[3]
    Rmin, Rmax, Zmin, Zmax = field.domain
    Rmin <= R <= Rmax && Zmin <= Z <= Zmax ||
        throw(DomainError(pos, "position is outside the reconstruction domain"))
    values = first(_basis_terms(field, R, Z))
    return dot(values, _psi_coefficient_vector(field))
end

function _cylindrical_field(field::AxisymmetricFieldReconstruction, R, Z)
    values, dR, dZ, _ = _basis_terms(field, R, Z)
    coefficients = _psi_coefficient_vector(field)
    B_R = -dot(dZ, coefficients) / R
    B_phi = dot(values[1:length(field.toroidal_coefficients)], field.toroidal_coefficients)
    B_Z = dot(dR, coefficients) / R
    return B_R, B_phi, B_Z
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

    values, dR, dZ, gs = _basis_terms(field, R, z)
    n = length(field.toroidal_coefficients)
    c = field.toroidal_coefficients
    B_phi = dot(values[1:n], c)
    J_R = -dot(dZ[1:n], c) / μ₀
    J_phi = -dot(gs, _psi_coefficient_vector(field)) / (μ₀ * R)
    J_Z = (dot(dR[1:n], c) + B_phi / R) / μ₀
    cos_phi, sin_phi = x / R, y / R
    return SVector(
        J_R * cos_phi - J_phi * sin_phi,
        J_R * sin_phi + J_phi * cos_phi,
        J_Z
    )
end
