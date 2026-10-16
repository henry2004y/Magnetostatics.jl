"""
    GradShafranovReconstruction <: AbstractMagneticField

Translationally invariant magnetohydrostatic field reconstructed from a spacecraft
crossing. In its local frame, `B = (A_y, -A_x, B_z(A))` and
`laplacian(A) = -μ₀ * dP_t/dA`. This is planar spacecraft GSR, not the
cylindrical flux interpolation of [`AxisymmetricFieldReconstruction`](@ref).

Evaluate at Cartesian positions in metres; fields are returned in tesla in the
input measurement frame. The finite domain applies to local x and y only.
The local z coordinate is unrestricted. `diagnostics` describes track residuals,
source-function extrapolation and the sampled GS residual, not uncertainty.
The interpolated flux is C1; its analytic current is piecewise continuous in y.
"""
struct GradShafranovReconstruction{T, D} <: AbstractMagneticField
    frame::SMatrix{3, 3, T, 9}
    origin::SVector{3, T}
    x_halfwidth::T
    y::Vector{T}
    flux_coefficients::Matrix{T}
    normal_coefficients::Matrix{T}
    flux_offset::T
    flux_scale::T
    transverse_coefficients::Vector{T}
    axial_coefficients::Vector{T}
    diagnostics::D
end

function _chebyshev_sum(coefficients, u)
    value = coefficients[1]
    derivative, second = zero(value), zero(value)
    length(coefficients) == 1 && return value, derivative, second
    previous, term = one(u), u
    dprevious, dterm = zero(u), one(u)
    ddprevious, ddterm = zero(u), zero(u)
    value += coefficients[2] * term
    derivative += coefficients[2]
    for k in 3:length(coefficients)
        previous, term = term, 2u * term - previous
        ddprevious, ddterm = ddterm, 4dterm + 2u * ddterm - ddprevious
        dprevious, dterm = dterm, 2previous + 2u * dterm - dprevious
        value += coefficients[k] * term
        derivative += coefficients[k] * dterm
        second += coefficients[k] * ddterm
    end
    return value, derivative, second
end

function _chebyshev_matrices(coordinates, degree)
    T = eltype(coordinates)
    values = Matrix{T}(undef, length(coordinates), degree + 1)
    second = similar(values)
    coefficient = zeros(T, degree + 1)
    for k in eachindex(coefficient)
        coefficient[k] = one(T)
        for (i, u) in enumerate(coordinates)
            values[i, k], _, second[i, k] = _chebyshev_sum(coefficient, u)
        end
        coefficient[k] = zero(T)
    end
    return values, second
end

function _gsr_terms(field::GradShafranovReconstruction, pos)
    length(pos) == 3 && all(isfinite, pos) ||
        throw(ArgumentError("position must contain three finite coordinates"))
    local_position = field.frame' * (SVector{3}(pos) - field.origin)
    x, y = local_position[1], local_position[2]
    tolerance = 16eps(eltype(field.frame)) *
        max(field.x_halfwidth, maximum(abs, field.y))
    abs(x) <= field.x_halfwidth + tolerance &&
        first(field.y) - tolerance <= y <= last(field.y) + tolerance ||
        throw(DomainError(pos, "position is outside the GSR reconstruction strip"))
    x = clamp(x, -field.x_halfwidth, field.x_halfwidth)
    y = clamp(y, first(field.y), last(field.y))
    j = min(searchsortedlast(field.y, y), length(field.y) - 1)
    h = field.y[j + 1] - field.y[j]
    t = (y - field.y[j]) / h
    h0, h1 = 2t^3 - 3t^2 + 1, -2t^3 + 3t^2
    g0, g1 = t^3 - 2t^2 + t, t^3 - t^2
    dh0, dh1 = 6t^2 - 6t, -6t^2 + 6t
    dg0, dg1 = 3t^2 - 4t + 1, 3t^2 - 2t
    ddh0, ddh1 = 12t - 6, -12t + 6
    ddg0, ddg1 = 6t - 4, 6t - 2
    u = x / field.x_halfwidth
    a0 = _chebyshev_sum(view(field.flux_coefficients, :, j), u)
    a1 = _chebyshev_sum(view(field.flux_coefficients, :, j + 1), u)
    q0 = _chebyshev_sum(view(field.normal_coefficients, :, j), u)
    q1 = _chebyshev_sum(view(field.normal_coefficients, :, j + 1), u)
    interpolate(k) = h0 * a0[k] + h1 * a1[k] + h * (g0 * q0[k] + g1 * q1[k])
    A = interpolate(1)
    Ax = interpolate(2) / field.x_halfwidth
    Axx = interpolate(3) / field.x_halfwidth^2
    Ay = (dh0 * a0[1] + dh1 * a1[1]) / h + dg0 * q0[1] + dg1 * q1[1]
    Ayy = (ddh0 * a0[1] + ddh1 * a1[1]) / h^2 +
        (ddg0 * q0[1] + ddg1 * q1[1]) / h
    return A, Ax, Ay, Axx, Ayy
end

function (field::GradShafranovReconstruction)(pos)
    A, Ax, Ay, _, _ = _gsr_terms(field, pos)
    Bz, _, _ = _chebyshev_sum(
        field.axial_coefficients, (A - field.flux_offset) / field.flux_scale
    )
    return field.frame * SVector(Ay, -Ax, Bz)
end

function current_density(field::GradShafranovReconstruction, pos)
    A, Ax, Ay, Axx, Ayy = _gsr_terms(field, pos)
    _, derivative, _ = _chebyshev_sum(
        field.axial_coefficients, (A - field.flux_offset) / field.flux_scale
    )
    dBz = derivative / field.flux_scale
    return field.frame * SVector(dBz * Ay, -dBz * Ax, -Axx - Ayy) / μ₀
end

"""
    flux_function(field::GradShafranovReconstruction, pos)

Return the planar flux function `A` [T m], with gauge `A = 0` at the first
track sample before spatial smoothing. Field lines projected into the
reconstruction plane follow contours of `A`.
"""
flux_function(field::GradShafranovReconstruction, pos) = first(_gsr_terms(field, pos))

"""
    thermal_pressure(field::GradShafranovReconstruction, pos)

Return the fitted isotropic thermal pressure `P_t(A) - B_z(A)^2/(2μ₀)` [Pa].
Extrapolated pressure functions are not guaranteed to remain physical.
"""
function thermal_pressure(field::GradShafranovReconstruction, pos)
    A = flux_function(field, pos)
    u = (A - field.flux_offset) / field.flux_scale
    Pt = first(_chebyshev_sum(field.transverse_coefficients, u))
    Bz = first(_chebyshev_sum(field.axial_coefficients, u))
    return Pt - Bz^2 / (2μ₀)
end

function _gsr_source(c, values, second, projector, pressure, offset, scale)
    flux = values * c
    source = map(flux) do A
        _, derivative, _ = _chebyshev_sum(pressure, (A - offset) / scale)
        μ₀ * derivative / scale
    end
    result = -projector * (second * c + source)
    all(isfinite, result) ||
        error("GSR continuation diverged; reduce strip width or spatial degree")
    return result
end

function _march_gsr!(
        flux, normal, ys, track_index, max_step, values, second, projector,
        pressure, offset, scale
    )
    for direction in (-1, 1)
        last_index = direction == -1 ? 1 : length(ys)
        for j in (track_index + direction):direction:last_index
            previous = j - direction
            c, q = copy(flux[:, previous]), copy(normal[:, previous])
            steps = max(1, ceil(Int, abs(ys[j] - ys[previous]) / max_step))
            h = (ys[j] - ys[previous]) / steps
            for _ in 1:steps
                acceleration = _gsr_source(
                    c, values, second, projector, pressure, offset, scale
                )
                midpoint_c, midpoint_q = c + (h / 2) * q, q + (h / 2) * acceleration
                midpoint_acceleration = _gsr_source(
                    midpoint_c, values, second, projector, pressure, offset, scale
                )
                c += h * midpoint_q
                q += h * midpoint_acceleration
                all(isfinite, c) && all(isfinite, q) ||
                    error("GSR continuation diverged; reduce strip width or spatial degree")
            end
            flux[:, j], normal[:, j] = c, q
        end
    end
    return nothing
end

"""
    reconstruct_grad_shafranov(times, B_observed, pressure;
        crossing_velocity, invariant_axis, y, ...)

Reconstruct a planar magnetohydrostatic strip from a single straight crossing.
`times` [s] must increase strictly; `B_observed` is an `n × 3` Cartesian matrix
[T]; `pressure` contains nonnegative isotropic thermal pressures [Pa].
Missing/non-finite observations are not supported.

Supply the invariant axis and constant **spacecraft velocity relative to the
structure**, `crossing_velocity = v_spacecraft - V_HT` [m/s], in the measurement
frame. Local x follows its projection perpendicular to the invariant axis;
local y completes a right-handed frame. `origin` (default: zero) is the
Cartesian position at the midpoint of the first and last sample times.
The component of crossing velocity along the invariant axis has no effect.
This function does not estimate the HT frame, axis, or equilibrium validity.

`y` is an increasing vector of strip coordinates [m], including zero and at
least one point on each side. Trapezoidal integration gives `A(x,0) = -∫B_y dx`.
Chebyshev least squares fits `P_t(A) = pressure + B_z²/(2μ₀)` and `B_z(A)`,
with `pressure_degree = 3` and `axial_degree = 3` by default.
`spatial_degree = 8` truncates the along-track Chebyshev expansion and projects
the nonlinear GS source onto the same modes at every step. This is the spatial
regularization, not a physical uniqueness condition.

Second-order midpoint continuation marches in both y directions, with internal
steps no larger than `max_step` (default: half the minimum track spacing).
Cubic Hermite interpolation of the flux makes field evaluation divergence-free
between steps; the GS equation is satisfied only approximately. Currents use
the analytic curl of this interpolated field, not the fitted source alone.

Elliptic Cauchy continuation is ill-posed. Keep a narrow strip, test smaller y
steps and spatial degrees, and validate off-track when possible. A strip wider
than `0.2 * track_length` on either side generates a warning (a heuristic, not
a stability bound). `diagnostics` reports field and pressure-function residuals,
sampled GS residual [T/m], and the continued flux range. Polynomial source
extrapolation beyond 5% of the measured flux span generates a warning; it is
not an exponential-tail model or a reliable extrapolation prescription.
Multi-valued pressure/axial-field data produce fit residuals, not an axis search.
"""
function reconstruct_grad_shafranov(
        times::AbstractVector{<:Real}, B_observed::AbstractMatrix{<:Real},
        pressure::AbstractVector{<:Real};
        crossing_velocity, invariant_axis, y, origin = SVector(0.0, 0.0, 0.0),
        spatial_degree::Integer = 8, pressure_degree::Integer = 3,
        axial_degree::Integer = 3, max_step = nothing
    )
    n = length(times)
    n >= 5 && all(isfinite, times) && all(>(0), diff(times)) ||
        throw(ArgumentError("times must contain at least five finite increasing samples"))
    size(B_observed) == (n, 3) && all(isfinite, B_observed) ||
        throw(ArgumentError("B_observed must be a finite n × 3 matrix"))
    length(pressure) == n && all(isfinite, pressure) && all(>=(0), pressure) ||
        throw(ArgumentError("pressure must contain n finite nonnegative values"))
    for (name, vector) in (
            ("crossing_velocity", crossing_velocity), ("invariant_axis", invariant_axis),
            ("origin", origin),
        )
        length(vector) == 3 && all(isfinite, vector) ||
            throw(ArgumentError("$name must contain three finite coordinates"))
    end
    1 <= spatial_degree <= min(20, n - 1) ||
        throw(ArgumentError("spatial_degree must be between 1 and min(20, n - 1)"))
    1 <= pressure_degree <= min(20, n - 1) ||
        throw(ArgumentError("pressure_degree must be between 1 and min(20, n - 1)"))
    0 <= axial_degree <= min(20, n - 1) ||
        throw(ArgumentError("axial_degree must be between 0 and min(20, n - 1)"))
    y isa AbstractVector{<:Real} && length(y) >= 3 &&
        all(isfinite, y) && all(>(0), diff(y)) &&
        first(y) < 0 < last(y) && 0 in y ||
        throw(ArgumentError("y must increase, include zero, and extend to both sides"))
    T = promote_type(Float64, eltype(times), eltype(B_observed), eltype(pressure))
    axis = SVector{3, T}(invariant_axis)
    isfinite(norm(axis)) && norm(axis) > 0 ||
        throw(ArgumentError("invariant_axis must have a finite positive norm"))
    ez = normalize(axis)
    velocity = SVector{3, T}(crossing_velocity)
    transverse_velocity = velocity - dot(velocity, ez) * ez
    speed = norm(transverse_velocity)
    isfinite(speed) && isfinite(norm(velocity)) && speed > eps(T) * norm(velocity) ||
        throw(ArgumentError("crossing_velocity must have a nonzero in-plane component"))
    ex = transverse_velocity / speed
    frame = hcat(ex, cross(ez, ex), ez)
    elapsed = T.(times) .- first(times)
    x = speed .* (elapsed .- last(elapsed) / 2)
    halfwidth = (last(x) - first(x)) / 2
    isfinite(halfwidth) && halfwidth > 0 ||
        throw(ArgumentError("crossing must have a finite positive spatial extent"))
    local_B = Matrix{T}(B_observed) * frame
    A = zeros(T, n)
    for i in 2:n
        A[i] = A[i - 1] - (x[i] - x[i - 1]) * (local_B[i, 2] + local_B[i - 1, 2]) / 2
    end
    Amin, Amax = extrema(A)
    offset, scale = (Amin + Amax) / 2, (Amax - Amin) / 2
    scale > 0 && isfinite(scale) ||
        throw(ArgumentError("crossing must sample a nonconstant finite flux function"))
    source_coordinates = (A .- offset) ./ scale
    Pt = T.(pressure) + local_B[:, 3] .^ 2 / (2T(μ₀))
    all(isfinite, Pt) || throw(ArgumentError("transverse pressure must be finite"))
    P, _ = _chebyshev_matrices(source_coordinates, pressure_degree)
    Z, _ = _chebyshev_matrices(source_coordinates, axial_degree)
    rank(P) == size(P, 2) && rank(Z) == size(Z, 2) ||
        throw(ArgumentError("insufficient distinct flux samples for source-function fits"))
    pressure_coefficients, axial_coefficients = P \ Pt, Z \ local_B[:, 3]
    values, second = _chebyshev_matrices(x / halfwidth, spatial_degree)
    rank(values) == size(values, 2) ||
        throw(ArgumentError("track sampling cannot resolve the requested spatial degree"))
    second /= halfwidth^2
    projector = pinv(values)
    ys = T.(collect(y))
    step = max_step === nothing ? minimum(diff(x)) / 2 : T(max_step)
    isfinite(step) && step > 0 ||
        throw(ArgumentError("max_step must be finite and positive"))
    track_index = findfirst(iszero, ys)::Int
    flux = zeros(T, spatial_degree + 1, length(ys))
    normal = similar(flux)
    flux[:, track_index] = projector * A
    normal[:, track_index] = projector * local_B[:, 1]
    _march_gsr!(
        flux, normal, ys, track_index, step, values, second, projector,
        pressure_coefficients, offset, scale
    )
    arguments = (
        frame, SVector{3, T}(origin), halfwidth, ys, flux, normal, offset, scale,
        pressure_coefficients, axial_coefficients,
    )
    provisional = GradShafranovReconstruction(arguments..., (;))
    predicted = reduce(
        hcat, [
            frame' * provisional(provisional.origin + frame * SVector(xi, zero(T), zero(T)))
                for xi in x
        ]
    )'
    residuals = T[]
    continued = T[]
    pressure_samples = T[]
    y_samples = sort(vcat(ys, (ys[1:(end - 1)] + ys[2:end]) / 2))
    for xi in x, yi in y_samples
        point = provisional.origin + frame * SVector(xi, yi, zero(T))
        a, _, _, axx, ayy = _gsr_terms(provisional, point)
        push!(continued, a)
        _, derivative, _ = _chebyshev_sum(pressure_coefficients, (a - offset) / scale)
        push!(residuals, axx + ayy + T(μ₀) * derivative / scale)
        push!(pressure_samples, thermal_pressure(provisional, point))
    end
    all(isfinite, residuals) && all(isfinite, pressure_samples) ||
        error("GSR interpolation became non-finite; reduce strip width or spatial degree")
    continued_range = extrema(continued)
    diagnostics = (
        track_field_rms = norm(predicted - local_B) / sqrt(3n),
        transverse_pressure_rms = norm(P * pressure_coefficients - Pt) / sqrt(n),
        axial_field_rms = norm(Z * axial_coefficients - local_B[:, 3]) / sqrt(n),
        gs_residual_rms = norm(residuals) / sqrt(length(residuals)),
        measured_flux_range = (Amin, Amax), continued_flux_range = continued_range,
        extrapolated_fraction = count(a -> a < Amin || a > Amax, continued) /
            length(continued),
        minimum_pressure = minimum(pressure_samples),
        spatial_degree, pressure_degree, axial_degree,
    )
    if maximum(abs, ys) > T(0.4) * halfwidth
        @warn "GSR strip exceeds 0.2 times the track length; continuation may be unstable"
    end
    if continued_range[1] < Amin - T(0.1) * scale ||
            continued_range[2] > Amax + T(0.1) * scale
        message = "GSR source functions extrapolate beyond 5% of the measured flux span"
        @warn message measured = (Amin, Amax) continued = continued_range
    end
    if diagnostics.minimum_pressure < 0
        message = "GSR fitted thermal pressure becomes negative"
        @warn message minimum_pressure = diagnostics.minimum_pressure
    end
    return GradShafranovReconstruction(arguments..., diagnostics)
end
