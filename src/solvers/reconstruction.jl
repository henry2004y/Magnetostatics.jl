"""
    AxisymmetricFieldReconstruction <: AbstractMagneticField

Continuous axisymmetric magnetic field reconstructed from sparse measurements
using a Gaussian radial-basis expansion of the poloidal flux `ψ(R, Z)`.
The poloidal field is divergence-free by construction. The toroidal field is
fit independently when toroidal measurements are supplied.
"""
struct AxisymmetricFieldReconstruction{T} <: AbstractMagneticField
    centers::Matrix{T}
    psi_coefficients::Vector{T}
    toroidal_coefficients::Vector{T}
    length_scale::T
    domain::NTuple{4, T}
end

function _rbf_terms(R, Z, centers, length_scale)
    n = size(centers, 1)
    T = promote_type(typeof(R), typeof(Z), eltype(centers), typeof(length_scale))
    values = Vector{T}(undef, n)
    dR = similar(values)
    dZ = similar(values)
    delta_star = similar(values)
    inv_scale2 = inv(length_scale^2)

    @inbounds for i in 1:n
        dr = R - centers[i, 1]
        dz = Z - centers[i, 2]
        value = exp(-0.5 * (dr^2 + dz^2) * inv_scale2)
        values[i] = value
        dR[i] = -dr * inv_scale2 * value
        dZ[i] = -dz * inv_scale2 * value
        dRR = (dr^2 * inv_scale2^2 - inv_scale2) * value
        dZZ = (dz^2 * inv_scale2^2 - inv_scale2) * value
        delta_star[i] = dRR + dZZ - dR[i] / R
    end

    return values, dR, dZ, delta_star
end

function _append_fit_row!(rows, targets, uncertainties, row, target, uncertainty)
    isfinite(target) || return nothing
    isfinite(uncertainty) && uncertainty > 0 ||
        throw(ArgumentError("uncertainties must be finite and positive"))
    push!(rows, row)
    push!(targets, target)
    push!(uncertainties, uncertainty)
    return nothing
end

function _fit_rbf(
        rows, targets, uncertainties, regularization, ncenters,
        exact_rows, exact_targets
    )
    nrows = length(rows)
    T = promote_type(eltype(targets), eltype(uncertainties))
    A = zeros(T, nrows + ncenters, ncenters)
    b = zeros(T, nrows + ncenters)

    for i in 1:nrows
        inv_uncertainty = inv(uncertainties[i])
        A[i, :] .= rows[i] .* inv_uncertainty
        b[i] = targets[i] * inv_uncertainty
    end
    for i in 1:ncenters
        A[nrows + i, i] = sqrt(regularization)
    end

    isempty(exact_rows) && return A \ b

    C = reduce(vcat, [reshape(row, 1, ncenters) for row in exact_rows])
    d = exact_targets
    particular = pinv(C) * d
    residual = norm(C * particular - d, Inf)
    tolerance = 100 * eps(T) * max(norm(d, Inf), one(T))
    residual <= tolerance ||
        throw(ArgumentError("exact boundary flux constraints are inconsistent"))

    nullspace_basis = nullspace(C)
    size(nullspace_basis, 2) == 0 && return particular
    reduced = A * nullspace_basis
    correction = reduced \ (b - A * particular)
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
                                   regularization, uncertainties, ...)

Fit a continuous axisymmetric field to sparse magnetic-field samples. `positions`
is an `n × 2` matrix of `(R, Z)` coordinates and `B_observed` is an `n × 3`
matrix of cylindrical `(B_R, B_ϕ, B_Z)` components. Non-finite field components
are treated as missing measurements.

`boundary_positions` and `boundary_flux` impose exact poloidal-flux constraints.
Supply `boundary_uncertainties` to treat them as weighted observations instead.
`current_positions` and `current_density` add weighted toroidal-current samples
through the Grad–Shafranov equation. The domain is `(Rmin, Rmax, Zmin, Zmax)`;
evaluation outside it throws `DomainError`.
"""
function reconstruct_axisymmetric_field(
        positions::AbstractMatrix{<:Real}, B_observed::AbstractMatrix{<:Real};
        domain, length_scale, regularization = 1.0e-8,
        uncertainties = ones(size(B_observed)),
        basis_centers = nothing,
        boundary_positions = nothing, boundary_flux = nothing,
        boundary_uncertainties = nothing,
        current_positions = nothing, current_density = nothing,
        current_uncertainties = nothing
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
    psi_rows = Vector{Vector{T}}()
    psi_targets = T[]
    psi_uncertainties = T[]
    exact_boundary_rows = Vector{Vector{T}}()
    exact_boundary_targets = T[]
    toroidal_rows = Vector{Vector{T}}()
    toroidal_targets = T[]
    toroidal_uncertainty_values = T[]

    for i in axes(sample_positions, 1)
        R, Z = sample_positions[i, 1], sample_positions[i, 2]
        values, dR, dZ, _ = _rbf_terms(R, Z, centers, T(length_scale))
        for component in 1:3
            target = observed[i, component]
            isfinite(target) || continue
            sigma = uncertainty_values[i, component]
            if component == 1
                row = -dZ / R
                _append_fit_row!(psi_rows, psi_targets, psi_uncertainties,
                    row, target, sigma)
            elseif component == 2
                _append_fit_row!(toroidal_rows, toroidal_targets,
                    toroidal_uncertainty_values, values, target, sigma)
            else
                row = dR / R
                _append_fit_row!(psi_rows, psi_targets, psi_uncertainties,
                    row, target, sigma)
            end
        end
    end

    if boundary_positions !== nothing || boundary_flux !== nothing
        boundary_positions !== nothing && boundary_flux !== nothing ||
            throw(ArgumentError(
                "boundary_positions and boundary_flux must be given together"
            ))
        size(boundary_positions, 1) == length(boundary_flux) ||
            throw(ArgumentError("boundary_positions and boundary_flux lengths differ"))
        _validate_rz_positions(boundary_positions, domain_typed, "boundary_positions")
        sigma = boundary_uncertainties === nothing ? nothing : T.(boundary_uncertainties)
        sigma === nothing || length(sigma) == length(boundary_flux) ||
            throw(ArgumentError("boundary_uncertainties length differs"))
        for i in eachindex(boundary_flux)
            values = first(_rbf_terms(
                boundary_positions[i, 1], boundary_positions[i, 2],
                centers, T(length_scale)
            ))
            if sigma === nothing
                isfinite(boundary_flux[i]) ||
                    throw(ArgumentError("boundary_flux must be finite"))
                push!(exact_boundary_rows, values)
                push!(exact_boundary_targets, T(boundary_flux[i]))
            else
                _append_fit_row!(psi_rows, psi_targets, psi_uncertainties,
                    values, T(boundary_flux[i]), sigma[i])
            end
        end
    end

    if current_positions !== nothing || current_density !== nothing
        current_positions !== nothing && current_density !== nothing ||
            throw(ArgumentError(
                "current_positions and current_density must be given together"
            ))
        size(current_positions, 1) == length(current_density) ||
            throw(ArgumentError("current_positions and current_density lengths differ"))
        _validate_rz_positions(current_positions, domain_typed, "current_positions")
        sigma = current_uncertainties === nothing ?
            ones(T, length(current_density)) : T.(current_uncertainties)
        length(sigma) == length(current_density) ||
            throw(ArgumentError("current_uncertainties length differs"))
        for i in eachindex(current_density)
            R, Z = current_positions[i, 1], current_positions[i, 2]
            delta_star = last(_rbf_terms(R, Z, centers, T(length_scale)))
            target = -T(μ₀) * R * T(current_density[i])
            _append_fit_row!(psi_rows, psi_targets, psi_uncertainties,
                delta_star, target, sigma[i])
        end
    end

    isempty(psi_rows) && isempty(exact_boundary_rows) &&
        throw(ArgumentError("poloidal or boundary constraints are required"))
    psi_coefficients = _fit_rbf(
        psi_rows, psi_targets, psi_uncertainties, T(regularization), size(centers, 1),
        exact_boundary_rows, exact_boundary_targets
    )
    toroidal_coefficients = isempty(toroidal_rows) ? zeros(T, size(centers, 1)) :
        _fit_rbf(toroidal_rows, toroidal_targets, toroidal_uncertainty_values,
            T(regularization), size(centers, 1), Vector{Vector{T}}(), T[])

    return AxisymmetricFieldReconstruction(
        centers, psi_coefficients, toroidal_coefficients,
        T(length_scale), domain_typed
    )
end

function poloidal_flux(field::AxisymmetricFieldReconstruction, pos)
    R = length(pos) == 2 ? pos[1] : hypot(pos[1], pos[2])
    Z = length(pos) == 2 ? pos[2] : pos[3]
    Rmin, Rmax, Zmin, Zmax = field.domain
    Rmin <= R <= Rmax && Zmin <= Z <= Zmax ||
        throw(DomainError(pos, "position is outside the reconstruction domain"))
    values = first(_rbf_terms(R, Z, field.centers, field.length_scale))
    return dot(values, field.psi_coefficients)
end

function _cylindrical_field(field::AxisymmetricFieldReconstruction, R, Z)
    values, dR, dZ, _ = _rbf_terms(R, Z, field.centers, field.length_scale)
    B_R = -dot(dZ, field.psi_coefficients) / R
    B_phi = dot(values, field.toroidal_coefficients)
    B_Z = dot(dR, field.psi_coefficients) / R
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