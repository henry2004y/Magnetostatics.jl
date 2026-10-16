module test_reconstruction_cases

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra
using SparseArrays

point_matrix(points) = [p[k] for p in vec(collect(points)), k in 1:3]
samples(field, positions) = point_matrix(field(SVector{3}(p)) for p in eachrow(positions))

function lattice(box, coordinates)
    return point_matrix(Iterators.product(
        ntuple(k -> [
            box[2k - 1] + t * (box[2k] - box[2k - 1]) for t in coordinates
        ], 3)...
    ))
end

function relative_error(field, truth, positions)
    values = samples(truth, positions)
    return norm(samples(field, positions) - values) / norm(values)
end

field_checksum(field, point) = sum(field(point))
current_checksum(field, point) = sum(current_density(field, point))

function evaluation_allocations(field, point)
    field_checksum(field, point)
    current_checksum(field, point)
    return (
        @allocated(field_checksum(field, point)),
        @allocated(current_checksum(field, point)),
    )
end

@testset "Biot-Savart sources and physical boundaries" begin
    box = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
    wall_box = (-0.5, 0.5, -0.5, 0.5, 0.0, 0.6)
    sphere_box = (-0.25, 0.25, -0.25, 0.25, 0.9, 1.4)
    loop = discretize_loop(CurrentLoop(
        1.0, 1.0e4, SVector(0.2, -0.1, 2.0), normalize(SVector(0.3, 0.2, 1.0))
    ), 150)
    wire = Wire([SVector(0.0, -10.0, 1.5), SVector(0.0, 10.0, 1.5)], 1.0e4)
    background = UniformField(SVector(0.0, 0.0, 0.01))
    image = Dipole(SVector(0.0, 0.0, -2π * 0.01 / Magnetostatics.μ₀))
    sphere_truth(p) = background(p) + image(p)
    wall_points = point_matrix(
        SVector(x, y, 0.0) for x in range(-0.4, 0.4; length = 5),
            y in range(-0.4, 0.4; length = 5)
    )
    sphere_points = point_matrix(
        SVector(x, y, sqrt(1 - x^2 - y^2))
            for x in range(-0.2, 0.2; length = 5), y in range(-0.2, 0.2; length = 5)
    )
    for (name, truth, domain, boundary) in (
            (
                "Open loop", p -> solve(BiotSavart(), loop, p), box, (;),
            ),
            (
                "Conducting wall", p -> solve(BiotSavart(ConductingWall(3, 0.0)), wire, p),
                wall_box, (;
                    boundary_positions = wall_points,
                    boundary_normals = repeat([0.0 0.0 1.0], size(wall_points, 1)),
                    boundary_values = zeros(size(wall_points, 1)),
                ),
            ),
            (
                "Conducting sphere exterior", sphere_truth, sphere_box, (;
                    boundary_positions = sphere_points, boundary_normals = sphere_points,
                    boundary_values = zeros(size(sphere_points, 1)),
                    boundary_uncertainties = fill(1.0e-6, size(sphere_points, 1)),
                    reference_field = background,
                ),
            ),
        )
        @testset "$name" begin
            positions = lattice(domain, (0.05, 0.35, 0.65, 0.95))
            heldout = lattice(domain, (0.17, 0.53, 0.81))
            if name == "Conducting sphere exterior"
                positions = positions[[norm(p) > 1.02 for p in eachrow(positions)], :]
                heldout = heldout[[norm(p) > 1.02 for p in eachrow(heldout)], :]
            end
            observed = samples(truth, positions)
            fitted = reconstruct_field_3d(
                positions, observed;
                domain, length_scale = name == "Open loop" ? 1.2 : 0.8,
                kernel = :gaussian,
                basis_centers = lattice(domain, range(0, 1; length = 4)),
                uncertainties = fill(1.0e-5, size(observed)),
                regularization = 1.0e-8, regularization_grid = (5, 5, 5),
                abstol = 1.0e-13, reltol = 1.0e-13, boundary...
            )
            error = relative_error(fitted, truth, heldout)
            @test error < (name == "Conducting sphere exterior" ? 0.05 : 0.03)
            @test fitted.diagnostics.unsupported_observations == 0
            @test fitted.diagnostics.sampled_support_fraction == 1
            if name == "Conducting wall"
                @test fitted.diagnostics.exact_boundary_max_error < 1.0e-12
                @test maximum(
                    abs(fitted(SVector(x, y, 0.0))[3])
                        for x in (-0.27, 0.21), y in (-0.19, 0.31)
                ) < 1.0e-5
            elseif name == "Conducting sphere exterior"
                @test maximum(
                    abs(dot(p, fitted(p))) for p in (
                        SVector(x, y, sqrt(1 - x^2 - y^2))
                            for x in (-0.13, 0.17), y in (-0.11, 0.15)
                    )
                ) < 1.0e-4
                @test all(abs(dot(p, truth(p))) < 1.0e-16 for p in eachrow(sphere_points))
            end
        end
    end
end

@testset "Distributed volume current" begin
    domain = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
    truth = ZPinch(1.0e4, 2.0)
    positions = lattice(domain, (0.05, 0.35, 0.65, 0.95))
    heldout = lattice(domain, (0.17, 0.53, 0.81))
    observed = samples(truth, positions)
    current_positions = lattice(domain, (0.1, 0.5, 0.9))
    currents = point_matrix(current_density(truth, p) for p in eachrow(current_positions))
    fitted = reconstruct_field_3d(
        positions, observed;
        domain, length_scale = 1.2, kernel = :gaussian,
        basis_centers = lattice(domain, range(0, 1; length = 4)),
        uncertainties = fill(1.0e-5, size(observed)),
        current_positions, current_density = currents,
        current_uncertainties = fill(1.0, size(currents)),
        regularization = 1.0e-8, regularization_grid = (5, 5, 5),
        abstol = 1.0e-12, reltol = 1.0e-12
    )
    @test relative_error(fitted, truth, heldout) < 0.01
    expected = point_matrix(current_density(truth, p) for p in eachrow(heldout))
    predicted = point_matrix(current_density(fitted, p) for p in eachrow(heldout))
    @test norm(predicted - expected) / norm(expected) < 0.02
end

@testset "Direct evaluation agrees with sparse operators" begin
    domain = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
    centers = lattice(domain, range(0, 1; length = 4))
    coefficients = sin.(1:(3 + 3size(centers, 1)))
    for kernel in (:gaussian, :wendland_c2)
        field = FieldReconstruction3D(
            centers, coefficients, 0.7, domain, kernel,
            Magnetostatics._bins_3d(centers, domain, 0.7),
            UniformField(SVector(0.01, -0.02, 0.03)), (;)
        )
        for p in (
                SVector(0.17, -0.23, 0.31), SVector{3}(centers[1, :]),
                SVector(1.0, 1.0, 1.0), SVector(-0.3, -1.0, -1.0),
            )
            B, curl, _ = Magnetostatics._operators_3d(field, p)
            @test field(p) ≈ field.reference_field(p) +
                SVector{3}(map(row -> dot(row, coefficients), B))
            @test current_density(field, p) ≈
                SVector{3}(map(row -> dot(row, coefficients) / Magnetostatics.μ₀, curl))
        end
        rows = first(Magnetostatics._operators_3d(field, SVector(0.1, 0.2, 0.3)))
        @test Matrix(Magnetostatics._matrix_3d(rows, length(coefficients), Float64)) ==
            reduce(hcat, rows)'
        @test all(<(512), evaluation_allocations(field, SVector(0.1, 0.2, 0.3)))
    end
    empty_matrix = Magnetostatics._matrix_3d(SparseVector{Float64, Int}[], 4, Float64)
    @test size(empty_matrix) == (0, 4)
end

end # module test_reconstruction_cases
