module ReconstructionUncertaintyDemo

using Magnetostatics
using StaticArrays
using LinearAlgebra
using LinearSolve: QRFactorization
using Random

const box = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)

point_matrix(points) = [p[k] for p in vec(collect(points)), k in 1:3]
samples(field, positions) = [field(p)[k] for p in eachrow(positions), k in 1:3]

function null_field(; kernel = :gaussian, radius = 0.5, scale = 1.5, amplitude = 1.0e-3)
    centers = point_matrix(Iterators.product(ntuple(_ -> (-radius, radius), 3)...))
    coefficients = zeros(3 + 3size(centers, 1))
    for (j, p) in enumerate(eachrow(centers))
        sx, sy, sz = sign.(p)
        coefficients[3j + 2] = 2amplitude * sx * sz
        coefficients[3j + 3] = amplitude * sx * sy
    end
    bins = Dict{NTuple{3, Int}, Vector{Int}}()
    for (j, p) in enumerate(eachrow(centers))
        cell = ntuple(k -> floor(Int, (p[k] - box[2k - 1]) / scale), 3)
        push!(get!(bins, cell, Int[]), j)
    end
    return FieldReconstruction3D(
        centers, coefficients, scale, box, kernel, bins, NullField(), (;)
    )
end

function boundary_samples()
    positions, normals = SVector{3, Float64}[], SVector{3, Float64}[]
    for axis in 1:3, side in (-1.0, 1.0), a in (-0.35, 0.35), b in (-0.35, 0.35)
        others = filter(!=(axis), 1:3)
        p = zeros(3)
        p[axis] = side
        p[others] .= (a, b)
        push!(positions, SVector{3}(p))
        push!(normals, SVector{3}(ntuple(k -> k == axis ? side : 0.0, 3)))
    end
    return point_matrix(positions), point_matrix(normals)
end

function problem()
    truth = null_field()
    positions = point_matrix(
        Iterators.product((-0.8, 0.1, 0.8), (-0.7, 0.4), (-0.6, 0.2, 0.6))
    )
    heldout = point_matrix(
        Iterators.product((-0.43, 0.13, 0.47), (-0.31, 0.17), (-0.23, 0.29))
    )
    surface, normals = boundary_samples()
    boundary = [dot(n, truth(p)) for (p, n) in zip(eachrow(surface), eachrow(normals))]
    return (; truth, positions, heldout, surface, normals, boundary)
end

function fit(data, observed, boundary; sigma, constrained)
    constraints = constrained ? (
            boundary_positions = data.surface, boundary_normals = data.normals,
            boundary_values = boundary,
            boundary_uncertainties = fill(sigma / 2, length(boundary)),
        ) : (;)
    return reconstruct_field_3d(
        data.positions, observed; domain = box, basis_centers = data.truth.centers,
        length_scale = data.truth.length_scale, kernel = :gaussian,
        uncertainties = fill(sigma, size(observed)), regularization = 1.0e-4 / sigma^2,
        regularization_mode = :coefficients, abstol = 1.0e-13, reltol = 1.0e-13,
        algorithm = QRFactorization(),
        constraints...
    )
end

function jacobian(field, point; h = 1.0e-5)
    return hcat(
        (
            let step = SVector{3}(ntuple(j -> j == k ? h : 0.0, 3))
                    (field(point + step) - field(point - step)) / (2h)
            end for k in 1:3
        )...
    )
end

function locate_null(field; start = SVector(0.0, 0.0, 0.0), h = 1.0e-5)
    point = start
    for _ in 1:20
        residual = norm(field(point))
        derivative = jacobian(field, point; h)
        singular_values = svdvals(derivative)
        if minimum(singular_values) <= 1.0e-8 * maximum(singular_values)
            return (; point, residual, status = :degenerate, positive = 0)
        end
        if residual < 1.0e-11
            positive = count(v -> real(v) > 0, eigvals(Matrix(derivative)))
            return (; point, residual, status = :converged, positive)
        end
        candidate = point - derivative \ field(point)
        if any(abs.(candidate) .> 0.4)
            return (; point, residual, status = :outside_search, positive = 0)
        end
        point = candidate
    end
    return (; point, residual = norm(field(point)), status = :maxiters, positive = 0)
end

function prediction_summary(fields, positions)
    values = cat((samples(field, positions) for field in fields)...; dims = 3)
    average = sum(values; dims = 3)[:, :, 1] / length(fields)
    variance = sum(
        (values[:, :, j] - average) .^ 2 for j in axes(values, 3)
    ) / length(fields)
    return (; average, variance, rms_spread = sqrt(sum(variance) / length(variance)))
end

function prediction_covariance(data, point; sigma, constrained)
    ncoefficients = length(data.truth.coefficients)
    basis = map(1:ncoefficients) do j
        coefficients = zeros(ncoefficients)
        coefficients[j] = 1
        t = data.truth
        FieldReconstruction3D(
            t.centers, coefficients, t.length_scale, t.domain, t.kernel,
            t.center_bins, NullField(), (;)
        )
    end
    operator = hcat((vec(samples(field, data.positions)) for field in basis)...) / sigma
    if constrained
        boundary_operator = [
            dot(n, field(p)) for (p, n) in zip(
                    eachrow(data.surface), eachrow(data.normals)
                ), field in basis
        ] / (sigma / 2)
        operator = vcat(operator, boundary_operator)
    end
    n = size(operator, 1)
    ridge = (1.0e-2 / sigma) * Matrix{Float64}(I, ncoefficients, ncoefficients)
    response = qr(vcat(operator, ridge)) \ vcat(
        Matrix{Float64}(I, n, n), zeros(ncoefficients, n)
    )
    evaluation = hcat((field(point) for field in basis)...)
    sensitivity = evaluation * response
    return sensitivity * sensitivity'
end

function run_ensemble(; pairs = 64, sigma = 1.0e-5, seed = 20261016)
    pairs >= 2 || throw(ArgumentError("at least two antithetic pairs are required"))
    sigma > 0 && isfinite(sigma) ||
        throw(ArgumentError("sigma must be positive and finite"))
    data = problem()
    rng = Xoshiro(seed)
    observed = samples(data.truth, data.positions) +
        sigma * randn(rng, size(data.positions))
    boundary = data.boundary + (sigma / 2) * randn(rng, length(data.boundary))
    perturbations = [
        (randn(rng, size(observed)), randn(rng, length(boundary))) for _ in 1:pairs
    ]
    results = map((false, true)) do constrained
        nominal = fit(data, observed, boundary; sigma, constrained)
        fields = [
            fit(
                    data, observed + side * sigma * noise,
                    boundary + side * (sigma / 2) * boundary_noise; sigma, constrained
                )
                for (noise, boundary_noise) in perturbations for side in (-1, 1)
        ]
        summary = prediction_summary(fields, data.heldout)
        nulls = [locate_null(field) for field in fields]
        return (; nominal, fields, summary, nulls)
    end
    return (; data, observed, boundary, sigma, free = results[1], bounded = results[2])
end

struct PerturbedField{F, G} <: AbstractMagneticField
    baseline::F
    perturbation::G
end

(field::PerturbedField)(p) = field.baseline(p) + field.perturbation(p)

function ambiguity_case()
    data = problem()
    hidden = null_field(;
        kernel = :wendland_c2, radius = 0.04, scale = 0.25, amplitude = -1.0e-4
    )
    alternative = PerturbedField(data.truth, hidden)
    positive_root = locate_null(alternative; start = SVector(0.2, 0.0, 0.0))
    negative_root = locate_null(alternative; start = SVector(-0.2, 0.0, 0.0))
    return (; data, hidden, alternative, positive_root, negative_root)
end

end
