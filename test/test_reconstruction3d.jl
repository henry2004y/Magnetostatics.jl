module test_reconstruction3d

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

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
    scale = 1.5
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
        for p in (SVector(0.13, -0.17, 0.29), SVector(-0.81, 0.41, 0.08))
            @test reconstructed(p) ≈ truth(p) atol = 1.0e-9
            @test abs(divergence(reconstructed, p)) < 1.0e-8
            @test current_density(reconstructed, p) ≈ current_density(truth, p) rtol = 1e-7
            h = 1.0e-5
            gradient = hcat([
                begin
                    step = SVector{3}(ntuple(i -> i == k ? h : 0.0, 3))
                    (reconstructed(p + step) - reconstructed(p - step)) / (2h)
                end for k in 1:3
            ]...)
            curl = SVector(
                gradient[3, 2] - gradient[2, 3],
                gradient[1, 3] - gradient[3, 1], gradient[2, 1] - gradient[1, 2]
            )
            @test current_density(reconstructed, p) ≈ curl / Magnetostatics.μ₀ rtol = 1e-7
        end
    end
end

end # module test_reconstruction3d
