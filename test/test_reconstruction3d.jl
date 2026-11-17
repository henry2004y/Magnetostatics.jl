module test_reconstruction3d

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra
using LinearSolve

const domain = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)

function point_matrix(points)
    return reduce(vcat, [reshape(collect(p), 1, 3) for p in points])
end

function samples(field, positions)
    return point_matrix([field(SVector{3}(p)) for p in eachrow(positions)])
end

function divergence(field, p; h = 1.0e-5)
    return sum(1:3) do k
        step = SVector{3}(ntuple(i -> i == k ? h : 0.0, 3))
        (field(p + step)[k] - field(p - step)[k]) / (2h)
    end
end

function manufactured(; kernel = :gaussian, scale = 2.0)
    centers = [-0.4 0.1 0.2; 0.3 -0.2 -0.3]
    coefficients = [0.01, 0.02, -0.01, 0.02, -0.03, 0.01, -0.01, 0.01, 0.04]
    return FieldReconstruction3D(
        centers, coefficients, scale, domain, kernel,
        Magnetostatics._bins_3d(centers, domain, scale), NullField(), (;)
    )
end

function cube_surface(coordinates)
    points, normals = SVector{3, Float64}[], SVector{3, Float64}[]
    for axis in 1:3, sign in (-1.0, 1.0), a in coordinates, b in coordinates
        other_axes = filter(!=(axis), 1:3)
        point = zeros(3)
        point[axis] = sign
        point[other_axes] .= (a, b)
        push!(points, SVector{3}(point))
        push!(normals, SVector{3}(ntuple(k -> k == axis ? sign : 0.0, 3)))
    end
    return point_matrix(points), point_matrix(normals)
end

struct InvalidReference <: AbstractMagneticField end
(::InvalidReference)(pos) = SVector(NaN, 0.0, 0.0)

@testset "Uniform field and reference" begin
    positions = [0.1 0.2 0.3; -0.5 0.3 0.2; 0.7 -0.3 -0.2; 0.2 -0.7 0.8]
    true_field = UniformField(SVector(0.01, -0.02, 0.03))
    centers = zeros(1, 3)
    reconstructed = reconstruct_field_3d(
        positions, samples(true_field, positions);
        domain, length_scale = 2.0, basis_centers = centers,
        regularization = 1.0e-14
    )
    @test reconstructed(SVector(-0.8, 0.6, -0.4)) ≈ true_field(zeros(3)) atol = 1.0e-10
    @test norm(current_density(reconstructed, zeros(3))) < 1.0e-6
    @test reconstructed(0.1, 0.2, 0.3) ≈ true_field(zeros(3))
    output = zeros(3)
    @test reconstructed(output, zeros(3)) ≈ true_field(zeros(3))
    @test reconstructed.diagnostics.weighted_rms < 1.0e-10
    reference_fit = reconstruct_field_3d(
        positions, samples(true_field, positions);
        domain, length_scale = 2.0, basis_centers = centers, reference_field = true_field
    )
    @test norm(reference_fit.coefficients) < 1.0e-12
    @test reference_fit(zeros(3)) == true_field(zeros(3))
    @test current_density(reference_fit, zeros(3)) == zeros(3)
    @test_throws DomainError reconstructed([1.1, 0.0, 0.0])
    @test_throws DomainError current_density(reconstructed, [0.0, NaN, 0.0])
    @test_throws ArgumentError reconstructed([0.0, 0.0])
end

@testset "Manufactured curl field" begin
    centers = [-0.4 0.1 0.2; 0.3 -0.2 -0.3]
    scale = 2.0
    coefficients = [0.01, 0.02, -0.01, 0.02, -0.03, 0.01, -0.01, 0.01, 0.04]
    for kernel in (:gaussian, :wendland_c2)
        truth = FieldReconstruction3D(
            centers, coefficients, scale, domain, kernel,
            Magnetostatics._bins_3d(centers, domain, scale), NullField(), (;)
        )
        positions = point_matrix(
            Iterators.product((-0.7, 0.0, 0.7), (-0.6, 0.2, 0.6), (-0.5, 0.1, 0.7))
        )
        reconstructed = reconstruct_field_3d(
            positions, samples(truth, positions);
            domain, length_scale = scale, basis_centers = centers, kernel,
            regularization = 1.0e-14, abstol = 1.0e-13, reltol = 1.0e-13
        )
        @test reconstructed.diagnostics.unsupported_grid_points == 0
        for p in (SVector(0.13, -0.17, 0.29), SVector(-0.81, 0.41, 0.08))
            @test reconstructed(p) ≈ truth(p) atol = 1.0e-9
            @test abs(divergence(reconstructed, p)) < 1.0e-8
            @test current_density(reconstructed, p) ≈
                current_density(truth, p) rtol = 1.0e-7
            h = 1.0e-5
            gradient = hcat(
                [
                    begin
                        step = SVector{3}(ntuple(i -> i == k ? h : 0.0, 3))
                        (reconstructed(p + step) - reconstructed(p - step)) / (2h)
                    end for k in 1:3
                ]...
            )
            curl = SVector(
                gradient[3, 2] - gradient[2, 3],
                gradient[1, 3] - gradient[3, 1], gradient[2, 1] - gradient[1, 2]
            )
            @test current_density(reconstructed, p) ≈ curl / Magnetostatics.μ₀ rtol = 1.0e-7
        end
    end
end

@testset "Normal boundary observations" begin
    truth = manufactured()
    positions = point_matrix(Iterators.product((-0.7, 0.2), (-0.6, 0.4), (-0.5, 0.6)))
    surface, normals = cube_surface((-0.6, 0.0, 0.6))
    boundary = [dot(n, truth(p)) for (p, n) in zip(eachrow(surface), eachrow(normals))]
    for exact in (true, false)
        fitted = reconstruct_field_3d(
            positions, samples(truth, positions);
            domain, length_scale = truth.length_scale, basis_centers = truth.centers,
            kernel = :gaussian, regularization = 1.0e-14,
            boundary_positions = surface, boundary_normals = normals,
            boundary_values = boundary,
            boundary_uncertainties = exact ? nothing : fill(1.0e-4, length(boundary))
        )
        @test fitted.diagnostics.exact_boundary_max_error < 1.0e-12
        @test fitted.diagnostics.exact_constraint_count == (exact ? length(boundary) : 0)
        if !exact
            @test fitted.diagnostics.boundary_rms < 1.0e-8
        end
        heldout, heldout_normals = cube_surface((-0.43, 0.27))
        @test maximum(
            abs(dot(n, fitted(p) - truth(p)))
                for (p, n) in zip(eachrow(heldout), eachrow(heldout_normals))
        ) < 1.0e-8
    end
    conflicting = repeat(surface[1:1, :], 2, 1)
    @test_throws ArgumentError reconstruct_field_3d(
        positions, samples(truth, positions);
        domain, length_scale = 2.0, basis_centers = truth.centers,
        boundary_positions = conflicting, boundary_normals = repeat(normals[1:1, :], 2, 1),
        boundary_values = [0.0, 1.0]
    )
    boundary_only = reconstruct_field_3d(
        zeros(0, 3), zeros(0, 3);
        domain, length_scale = 2.0, basis_centers = truth.centers, kernel = :gaussian,
        boundary_positions = surface, boundary_normals = normals,
        boundary_values = boundary,
        regularization = 1.0e-14
    )
    @test boundary_only.diagnostics.weighted_rms === nothing
    @test boundary_only.diagnostics.observation_count == 0
    @test boundary_only(SVector(0.21, -0.18, 0.37)) ≈
        truth(SVector(0.21, -0.18, 0.37)) atol = 1.0e-8
    reference = UniformField(SVector(0.02, -0.01, 0.03))
    uniform_boundary = [
        dot(n, reference(p)) for (p, n) in zip(eachrow(surface), eachrow(normals))
    ]
    reference_fit = reconstruct_field_3d(
        zeros(0, 3), zeros(0, 3);
        domain, length_scale = 2.0, basis_centers = truth.centers,
        reference_field = reference, boundary_positions = surface,
        boundary_normals = normals,
        boundary_values = uniform_boundary
    )
    @test norm(reference_fit.coefficients) < 1.0e-10
end

@testset "Directional, noisy and missing probes" begin
    truth = manufactured()
    positions = point_matrix(
        Iterators.product((-0.7, 0.0, 0.7), (-0.6, 0.2, 0.6), (-0.5, 0.1, 0.7))
    )
    values = samples(truth, positions)
    directions = point_matrix(
        [
            normalize(SVector(cos(i), sin(i), cos(2i))) for i in axes(positions, 1)
        ]
    )
    directional_values = [
        dot(n, b) for (n, b) in zip(eachrow(directions), eachrow(values))
    ]
    directional_fit = reconstruct_field_3d(
        zeros(0, 3), zeros(0, 3);
        domain, length_scale = truth.length_scale, basis_centers = truth.centers,
        kernel = :gaussian, directional_positions = positions,
        directional_directions = directions, directional_values,
        regularization = 1.0e-14, abstol = 1.0e-13, reltol = 1.0e-13
    )
    heldout = SVector(0.13, -0.17, 0.29)
    @test directional_fit(heldout) ≈ truth(heldout) atol = 1.0e-8
    @test directional_fit.diagnostics.directional_rms < 1.0e-9
    @test directional_fit.diagnostics.probe_rms === nothing
    missing = copy(values)
    missing[1:3:end, 1] .= NaN
    missing[2:3:end, 2] .= NaN
    sigma = fill(1.0e-4, size(values))
    sigma[isnan.(missing)] .= NaN
    missing_fit = reconstruct_field_3d(
        positions, missing;
        domain, length_scale = truth.length_scale, basis_centers = truth.centers,
        kernel = :gaussian, uncertainties = sigma, regularization = 1.0e-14
    )
    @test missing_fit(heldout) ≈ truth(heldout) atol = 1.0e-8
    @test missing_fit.diagnostics.observation_count == count(isfinite, missing)
    noisy = values + 1.0e-5 * reshape(sin.(1:length(values)), size(values))
    noisy[1, 1] += 0.02
    weighted_sigma = copy(sigma)
    weighted_sigma .= 1.0e-4
    weighted_sigma[1, 1] = 0.02
    fits = map((weighted_sigma, fill(1.0e-4, size(values)))) do uncertainty
        reconstruct_field_3d(
            positions, noisy;
            domain, length_scale = truth.length_scale, basis_centers = truth.centers,
            kernel = :gaussian, uncertainties = uncertainty, regularization = 1.0e-10
        )
    end
    @test norm(fits[1](heldout) - truth(heldout)) < norm(fits[2](heldout) - truth(heldout))
    @test norm(fits[1](heldout) - truth(heldout)) < 1.0e-4
end

@testset "Current observations and units" begin
    truth = manufactured()
    positions = point_matrix(Iterators.product((-0.7, 0.2), (-0.6, 0.4), (-0.5, 0.6)))
    currents = point_matrix([current_density(truth, p) for p in eachrow(positions)])
    noisy = samples(truth, positions)
    noisy[1, 1] += 0.01
    fits = map((1.0e-3, 1.0e9)) do uncertainty
        reconstruct_field_3d(
            positions, noisy;
            domain, length_scale = truth.length_scale, basis_centers = truth.centers,
            kernel = :gaussian, current_positions = positions, current_density = currents,
            current_uncertainties = fill(uncertainty, size(currents)),
            regularization = 1.0e-14
        )
    end
    heldout = SVector(0.13, -0.17, 0.29)
    @test norm(current_density(fits[1], heldout) - current_density(truth, heldout)) <
        norm(current_density(fits[2], heldout) - current_density(truth, heldout))
    @test fits[1].diagnostics.current_rms < 1.0e-2
    partial_currents = copy(currents)
    partial_currents[1, 2] = NaN
    reference_fit = reconstruct_field_3d(
        positions, samples(truth, positions);
        domain, length_scale = 2.0, basis_centers = truth.centers,
        reference_field = truth, kernel = :gaussian, current_positions = positions,
        current_density = partial_currents, current_uncertainties = ones(size(currents))
    )
    @test norm(reference_fit.coefficients) < 1.0e-10
    @test reference_fit.diagnostics.current_rms < 1.0e-10
    @test isapprox(
        current_density(reference_fit, heldout),
        current_density(truth, heldout);
        rtol = 1.0e-12, atol = 1.0e-10
    )
    current_only = reconstruct_field_3d(
        zeros(0, 3), zeros(0, 3);
        domain, length_scale = truth.length_scale, basis_centers = truth.centers,
        kernel = :gaussian, current_positions = positions, current_density = currents,
        current_uncertainties = fill(1.0e-3, size(currents)), regularization = 1.0e-14
    )
    @test current_density(current_only, heldout) ≈
        current_density(truth, heldout) rtol = 1.0e-6
end

@testset "Physical regularization and support diagnostics" begin
    truth = manufactured()
    energies = map((5, 9, 17)) do n
        L, gaps, count = Magnetostatics._regularizer_3d(
            truth.centers, truth.length_scale, truth.kernel, truth.domain,
            truth.center_bins, (n, n, n), :field, Float64
        )
        @test gaps == 0
        @test count == n^3
        points, weights = Magnetostatics._quadrature_3d(domain, (n, n, n), Float64)
        @test sum(weights) ≈ 1
        norm(L * truth.coefficients)^2
    end
    @test abs(energies[3] - energies[2]) < abs(energies[2] - energies[1])
    @test energies[2] ≈ energies[3] rtol = 0.03
    positions = [0.0 0.0 0.0; 0.8 0.8 0.8]
    values = repeat([0.01 -0.02 0.03], 2, 1)
    fit = @test_logs (:warn, r"RBF support has gaps") reconstruct_field_3d(
        positions, values;
        domain, length_scale = 0.2, basis_centers = zeros(1, 3),
        regularization_grid = (3, 3, 3)
    )
    @test fit.diagnostics.unsupported_observations == 1
    @test fit.diagnostics.unsupported_grid_points == 26
    @test fit.diagnostics.sampled_support_fraction ≈ 1 / 27
    @test fit(SVector(0.8, 0.8, 0.8)) ≈ SVector{3}(values[2, :]) atol = 1.0e-8
    covered = reconstruct_field_3d(
        positions, values; domain, length_scale = 1.5, regularization_grid = (5, 5, 5)
    )
    @test covered.diagnostics.sampled_support_fraction == 1
    @test covered.diagnostics.unsupported_observations == 0
    @test size(covered.centers, 1) == 4^3
    for p in (
            SVector(0.19, -0.23, 0.31), SVector(0.99, 0.99, 0.99),
            SVector(0.2 - 1.0e-4, 0.0, 0.0), SVector(0.2 + 1.0e-4, 0.0, 0.0),
        )
        @test abs(divergence(covered, p)) < 1.0e-7
        @test abs(divergence(fit, p; h = 1.0e-6)) < 1.0e-7
    end
    bins = Magnetostatics._bins_3d(zeros(1, 3), domain, 0.2)
    center_field = FieldReconstruction3D(
        zeros(1, 3), [0.0, 0.0, 0.0, 0.01, -0.02, 0.03],
        0.2, domain, :wendland_c2, bins, NullField(), (;)
    )
    @test all(isfinite, current_density(center_field, zeros(3)))
    @test abs(divergence(center_field, zeros(3); h = 1.0e-6)) < 1.0e-7
    @test abs(divergence(center_field, SVector(0.2, 0.0, 0.0); h = 1.0e-6)) < 1.0e-7
end

@testset "Independent analytical field and sensitivity" begin
    box = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
    loop = CurrentLoop(
        1.0, 1.0e4, SVector(0.2, -0.1, 2.0), normalize(SVector(0.3, 0.2, 1.0))
    )
    truth = CurrentLoopAnalytic(loop)
    positions = point_matrix(
        Iterators.product((-0.45, -0.15, 0.15, 0.45), (-0.4, 0.0, 0.4), (-0.4, 0.0, 0.4))
    )
    heldout = point_matrix(
        Iterators.product((-0.37, 0.11, 0.39), (-0.31, 0.17), (-0.23, 0.29))
    )
    reference_values = samples(truth, heldout)
    errors = Float64[]
    for (n, scale, penalty) in ((3, 1.6, 1.0e-10), (4, 1.2, 1.0e-10), (4, 1.6, 1.0e-8))
        centers = point_matrix(
            Iterators.product(ntuple(_ -> range(-0.5, 0.5; length = n), 3)...)
        )
        fitted = reconstruct_field_3d(
            positions, samples(truth, positions);
            domain = box, length_scale = scale, basis_centers = centers,
            kernel = :gaussian, regularization = penalty,
            uncertainties = fill(1.0e-5, size(positions)),
            abstol = 1.0e-12, reltol = 1.0e-12
        )
        relative_error = norm(samples(fitted, heldout) - reference_values) /
            norm(reference_values)
        push!(errors, relative_error)
        @test relative_error < 0.03
        @test maximum(abs(divergence(fitted, SVector{3}(p))) for p in eachrow(heldout)) <
            1.0e-8
    end
    @test errors[2] < errors[1]
    dipole_box = (-0.5, 0.5, -0.5, 0.5, 3.5, 4.5)
    dipole = Dipole(SVector(3.0e3, -2.0e3, 1.0e4))
    dipole_positions = copy(positions)
    dipole_positions[:, 3] .+= 4.0
    dipole_heldout = copy(heldout)
    dipole_heldout[:, 3] .+= 4.0
    dipole_centers = point_matrix(
        Iterators.product(
            range(-0.5, 0.5; length = 4), range(-0.5, 0.5; length = 4),
            range(3.5, 4.5; length = 4)
        )
    )
    dipole_fit = reconstruct_field_3d(
        dipole_positions, samples(dipole, dipole_positions);
        domain = dipole_box, length_scale = 1.2, basis_centers = dipole_centers,
        kernel = :gaussian, uncertainties = fill(1.0e-6, size(dipole_positions)),
        regularization = 1.0e-10, abstol = 1.0e-12, reltol = 1.0e-12
    )
    dipole_values = samples(dipole, dipole_heldout)
    @test norm(samples(dipole_fit, dipole_heldout) - dipole_values) /
        norm(dipole_values) < 0.03
end

@testset "Input validation" begin
    positions = [0.0 0.0 0.0]
    values = [0.01 -0.02 0.03]
    defaults = (; domain, length_scale = 2.0, basis_centers = zeros(1, 3))
    fit(; kwargs...) = reconstruct_field_3d(
        positions, values; merge(defaults, (; kwargs...))...
    )
    for kwargs in (
            (; domain = (0.0, 1.0)), (; domain = (0.0, 0.0, -1.0, 1.0, -1.0, 1.0)),
            (; length_scale = 0.0), (; length_scale = Inf), (; regularization = -1.0),
            (; abstol = NaN), (; reltol = 0.0), (; maxiters = 0), (; kernel = :invalid),
            (; regularization_mode = :invalid), (; regularization_grid = (1, 3, 3)),
            (; basis_centers = zeros(0, 3)), (; basis_centers = zeros(2, 2)),
            (; basis_centers = [2.0 0.0 0.0]), (; uncertainties = zeros(1, 3)),
            (; uncertainties = ones(2, 3)), (; reference_field = InvalidReference()),
            (; boundary_positions = positions), (; boundary_uncertainties = [1.0]),
            (; directional_uncertainties = [1.0]), (; current_uncertainties = ones(1, 3)),
            (; current_positions = positions),
            (; current_positions = positions, current_density = zeros(1, 2)),
            (;
                current_positions = positions, current_density = zeros(1, 3),
                current_uncertainties = zeros(1, 3),
            ),
            (;
                directional_positions = positions, directional_directions = [0.0 0.0 0.0],
                directional_values = [0.0],
            ),
            (;
                directional_positions = positions, directional_directions = [1.0 0.0 0.0],
                directional_values = [Inf],
            ),
        )
        @test_throws ArgumentError fit(; kwargs...)
    end
    @test_throws ArgumentError reconstruct_field_3d(zeros(1, 2), values; defaults...)
    @test_throws ArgumentError reconstruct_field_3d([NaN 0.0 0.0], values; defaults...)
    @test_throws ArgumentError reconstruct_field_3d(positions, fill(Inf, 1, 3); defaults...)
    @test_throws ArgumentError reconstruct_field_3d(positions, fill(NaN, 1, 3); defaults...)
    for kwargs in (
            (; boundary_normals = [2.0 0.0 0.0], boundary_values = [0.0]),
            (; boundary_normals = [1.0 0.0 0.0], boundary_values = [NaN]),
            (; boundary_normals = [1.0 0.0 0.0], boundary_values = [0.0, 0.0]),
            (;
                boundary_normals = [1.0 0.0 0.0], boundary_values = [0.0],
                boundary_uncertainties = [-1.0],
            ),
            (; boundary_normals = [NaN 0.0 0.0], boundary_values = [0.0]),
        )
        @test_throws ArgumentError fit(; boundary_positions = positions, kwargs...)
    end
end

@testset "Zero-normal wall, clustering and solver failures" begin
    centers = zeros(1, 3)
    scale = 0.6
    truth = FieldReconstruction3D(
        centers, [0.0, 0.0, 0.0, 0.01, -0.02, 0.03], scale, domain, :wendland_c2,
        Magnetostatics._bins_3d(centers, domain, scale), NullField(), (;)
    )
    positions = point_matrix(
        Iterators.product((-0.25, 0.1, 0.25), (-0.2, 0.15), (-0.1, 0.2))
    )
    surface, normals = cube_surface((-0.6, 0.0, 0.6))
    fit = @test_logs (:warn, r"RBF support has gaps") reconstruct_field_3d(
        positions, samples(truth, positions);
        domain, length_scale = scale, basis_centers = centers,
        boundary_positions = surface, boundary_normals = normals,
        boundary_values = zeros(size(surface, 1)), regularization = 1.0e-14
    )
    @test fit(SVector(0.1, 0.2, -0.1)) ≈ truth(SVector(0.1, 0.2, -0.1)) atol = 1.0e-8
    heldout, heldout_normals = cube_surface((-0.43, 0.27))
    @test maximum(
        abs(dot(n, fit(p))) for (p, n) in zip(eachrow(heldout), eachrow(heldout_normals))
    ) < 1.0e-12
    gaussian_truth = manufactured()
    cluster = [0.1 0.1 0.1; 0.11 0.1 0.1; 0.1 0.11 0.1]
    clustered = reconstruct_field_3d(
        cluster, samples(gaussian_truth, cluster);
        domain, length_scale = 2.0, basis_centers = gaussian_truth.centers,
        kernel = :gaussian, regularization = 1.0e-5
    )
    @test clustered.diagnostics.observation_count == 9
    @test all(isfinite, clustered(SVector(-0.8, 0.7, -0.6)))
    @test clustered.diagnostics.weighted_rms < 1.0e-4
    @test_throws ErrorException reconstruct_field_3d(
        positions, samples(gaussian_truth, positions);
        domain, length_scale = 2.0, basis_centers = gaussian_truth.centers,
        kernel = :gaussian, maxiters = 1
    )
    ridge_fit = reconstruct_field_3d(
        positions, samples(gaussian_truth, positions);
        domain, length_scale = 2.0, basis_centers = gaussian_truth.centers,
        kernel = :gaussian, regularization_mode = :coefficients,
        regularization = 1.0e-14
    )
    @test ridge_fit(SVector(0.13, -0.17, 0.29)) ≈
        gaussian_truth(SVector(0.13, -0.17, 0.29)) atol = 1.0e-7

    let
        eval_B(f, r) = (f(r); @allocated f(r))
        eval_J(f, r) = (current_density(f, r); @allocated current_density(f, r))
        rq = SVector(0.1, 0.1, 0.1)
        @test eval_B(fit, rq) == 0
        @test eval_J(fit, rq) == 0
        @test eval_B(ridge_fit, rq) == 0
        @test eval_J(ridge_fit, rq) == 0
    end
end

@testset "Curl-free harmonics, virtual dipoles, sources and covariance" begin
    let box = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
            eval_B(f, r) = (f(r); @allocated f(r)),
            eval_J(f, r) = (current_density(f, r); @allocated current_density(f, r))
        probes = point_matrix(
            Iterators.product(
                range(-0.4, 0.4; length = 4), range(-0.4, 0.4; length = 4),
                range(-0.4, 0.4; length = 4)
            )
        )
        held = point_matrix(
            Iterators.product((-0.31, 0.23), (-0.27, 0.19), (-0.29, 0.21))
        )

        # Pure curl-free solid harmonics on an off-domain dipole
        dip_ext = Dipole(SVector(1.5e4, -0.8e4, 2.2e4), SVector(0.15, -0.2, 2.1))
        cf_fit = reconstruct_field_3d(
            probes, samples(dip_ext, probes);
            domain = box, model = :curl_free, vacuum_degree = 6,
            regularization = 1.0e-12
        )
        @test cf_fit.model === :curl_free
        @test isempty(cf_fit.centers)
        @test cf_fit.vacuum_degree == 6
        @test length(cf_fit.vacuum_coefficients) == 60
        err_cf = norm(samples(cf_fit, held) - samples(dip_ext, held)) /
            norm(samples(dip_ext, held))
        @test err_cf < 1.0e-3
        pt = SVector(0.13, -0.17, 0.21)
        @test abs(divergence(cf_fit, pt; h = 1.0e-6)) / norm(cf_fit(pt)) < 1.0e-9
        @test current_density(cf_fit, pt) == SVector(0.0, 0.0, 0.0)
        @test eval_B(cf_fit, pt) == 0
        @test eval_J(cf_fit, pt) == 0
        @test_throws ArgumentError reconstruct_field_3d(
            probes, samples(dip_ext, probes);
            domain = box, model = :curl_free,
            current_positions = probes[1:2, :], current_density = zeros(2, 3)
        )

        # Virtual 3D external dipoles recover dipole moments accurately
        d_pos = [0.2 -0.3 1.9; -1.8 0.4 -0.6]
        d1 = Dipole(SVector(5.0e3, -2.0e3, 8.0e3), SVector(0.2, -0.3, 1.9))
        d2 = Dipole(SVector(-3.0e3, 4.0e3, 1.0e3), SVector(-1.8, 0.4, -0.6))
        truth_dips = d1 + d2
        dip_fit = reconstruct_field_3d(
            probes, samples(truth_dips, probes);
            domain = box, model = :curl_free, dipole_positions = d_pos,
            regularization = 1.0e-14, abstol = 1.0e-13, reltol = 1.0e-13
        )
        @test dip_fit.dipole_moments[1] ≈ d1.M rtol = 1.0e-7
        @test dip_fit.dipole_moments[2] ≈ d2.M rtol = 1.0e-7
        @test dip_fit(pt) ≈ truth_dips(pt) rtol = 1.0e-8
        @test eval_B(dip_fit, pt) == 0
        @test_throws ArgumentError reconstruct_field_3d(
            probes, samples(truth_dips, probes);
            domain = box, model = :curl_free, dipole_positions = [0.1 0.1 0.1]
        )

        # Parametric source_fields templates fit exact weights & current density
        s1 = HarrisSheet(1.0, 0.35)
        s2 = CurrentLoopAnalytic(
            CurrentLoop(0.8, 1.0, SVector(0.0, 0.0, 1.4), SVector(0.0, 0.0, 1.0))
        )
        truth_src = 0.75 * s1 + 1200.0 * s2
        src_fit = reconstruct_field_3d(
            probes, samples(truth_src, probes);
            domain = box, model = :curl_free, source_fields = (s1, s2),
            regularization = 1.0e-14, abstol = 1.0e-13, reltol = 1.0e-13
        )
        @test src_fit.source_weights[1] ≈ 0.75 rtol = 1.0e-8
        @test src_fit.source_weights[2] ≈ 1200.0 rtol = 1.0e-8
        @test current_density(src_fit, pt) ≈ current_density(truth_src, pt) rtol = 1.0e-8
        @test eval_B(src_fit, pt) == 0
        @test eval_J(src_fit, pt) == 0

        # Covariance propagation vs dense QR formula (unconstrained and constrained)
        @test_throws ArgumentError field_covariance(cf_fit, pt)
        sigma = 2.0e-5
        surf, nrms = cube_surface((-0.35, 0.35))
        surf_pts = 0.5 .* surf
        b_vals = [dot(n, dip_ext(p)) for (p, n) in zip(eachrow(surf_pts), eachrow(nrms))]
        fit_cov = reconstruct_field_3d(
            probes, samples(dip_ext, probes);
            domain = box, model = :curl_free, vacuum_degree = 3,
            uncertainties = fill(sigma, size(probes)),
            regularization = 1.0e-4 / sigma^2,
            regularization_mode = :coefficients,
            boundary_positions = surf_pts, boundary_normals = nrms,
            boundary_values = b_vals,
            boundary_uncertainties = fill(sigma / 2, length(b_vals)),
            store_covariance = true
        )
        cov_pt = field_covariance(fit_cov, pt)
        @test issymmetric(cov_pt)
        @test minimum(eigvals(Symmetric(cov_pt))) > 0
        # Double noise and rescale reg -> variance scales by 4x at matched reg
        fit_cov2 = reconstruct_field_3d(
            probes, samples(dip_ext, probes);
            domain = box, model = :curl_free, vacuum_degree = 3,
            uncertainties = fill(2 * sigma, size(probes)),
            regularization = 1.0e-4 / (2 * sigma)^2,
            regularization_mode = :coefficients,
            boundary_positions = surf_pts, boundary_normals = nrms,
            boundary_values = b_vals,
            boundary_uncertainties = fill(sigma, length(b_vals)),
            store_covariance = true
        )
        @test field_covariance(fit_cov2, pt) ≈ 4 * cov_pt rtol = 1.0e-10
        @test_throws DomainError field_covariance(fit_cov, SVector(1.0, 0.0, 0.0))

        # Exact boundary constraint covariance check (zero normal variance along normal)
        fit_exact_cov = reconstruct_field_3d(
            probes, samples(dip_ext, probes);
            domain = box, model = :curl_free, vacuum_degree = 3,
            uncertainties = fill(sigma, size(probes)),
            boundary_positions = surf_pts[1:1, :],
            boundary_normals = nrms[1:1, :],
            boundary_values = b_vals[1:1],
            regularization = 1.0e-8,
            regularization_mode = :coefficients,
            store_covariance = true
        )
        n1 = SVector{3}(nrms[1, :])
        p1 = SVector{3}(surf_pts[1, :])
        cov_bdy = field_covariance(fit_exact_cov, p1)
        @test dot(n1, cov_bdy * n1) < 1.0e-20

        # Directional and boundary constraints combined with virtual dipoles & sources
        truth_ext = d1 + 500.0 * s2
        dir_vecs = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        dir_pts = [0.2 -0.1 0.05; -0.15 0.25 -0.1; 0.05 0.1 -0.25]
        dir_vals = [
            dot(SVector{3}(d), truth_ext(SVector{3}(p)))
                for (p, d) in zip(eachrow(dir_pts), eachrow(dir_vecs))
        ]
        bdy_vals_ext = [
            dot(SVector{3}(n), truth_ext(SVector{3}(p)))
                for (p, n) in zip(eachrow(surf_pts), eachrow(nrms))
        ]
        fit_ext_dir = reconstruct_field_3d(
            probes, samples(truth_ext, probes);
            domain = box, model = :curl_free, vacuum_degree = 3,
            dipole_positions = d_pos[1:1, :], source_fields = (s2,),
            directional_positions = dir_pts, directional_directions = dir_vecs,
            directional_values = dir_vals, directional_uncertainties = fill(1.0e-3, 3),
            boundary_positions = surf_pts, boundary_normals = nrms,
            boundary_values = bdy_vals_ext,
            boundary_uncertainties = fill(1.0e-3, length(bdy_vals_ext)),
            regularization = 1.0e-14, algorithm = QRFactorization()
        )
        @test fit_ext_dir(pt) ≈ truth_ext(pt) rtol = 1.0e-4

        fit_def_vac = reconstruct_field_3d(
            probes, samples(dip_ext, probes); domain = box, model = :curl_free
        )
        @test fit_def_vac.vacuum_degree == 3
        @test_throws ArgumentError reconstruct_field_3d(
            probes, samples(dip_ext, probes); domain = box
        )
        @test_throws ArgumentError reconstruct_field_3d(
            probes, samples(dip_ext, probes);
            domain = box, model = :curl_free, basis_centers = probes[1:2, :]
        )
        @test_throws ArgumentError reconstruct_field_3d(
            probes, samples(dip_ext, probes); domain = box, model = :invalid
        )
    end
end

end # module test_reconstruction3d
