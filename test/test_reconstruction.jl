module test_reconstruction

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

@testset "Axisymmetric field reconstruction" begin
    domain = (0.4, 1.3, -0.6, 0.6)
    centers = reduce(vcat, [reshape([R, Z], 1, 2) for R in (0.5, 0.85, 1.2)
                            for Z in (-0.45, 0.0, 0.45)])
    length_scale = 0.3
    true_field = AxisymmetricFieldReconstruction(
        centers,
        1.0e-3 .* [1.0, -0.4, 0.8, 0.5, -1.2, 0.7, -0.6, 0.3, 1.1],
        0.02 .* [0.2, -0.3, 0.5, 0.7, 0.1, -0.4, 0.6, -0.2, 0.3],
        length_scale, domain
    )
    positions = [(R, Z) for R in range(0.48, 1.22; length = 5)
                          for Z in range(-0.5, 0.5; length = 4)]
    position_matrix = reduce(vcat, [reshape(collect(pos), 1, 2) for pos in positions])
    field_matrix = reduce(vcat, [reshape(collect(true_field(SVector(R, 0.0, Z))), 1, 3)
                                 for (R, Z) in positions])
    boundary_positions = reduce(vcat, [reshape([R, Z], 1, 2)
        for R in (domain[1], domain[2]), Z in (-0.5, 0.0, 0.5)])
    boundary_values = [
        poloidal_flux(true_field, SVector(R, 0.0, Z))
        for (R, Z) in eachrow(boundary_positions)
    ]
    current_positions = [0.6 -0.3; 0.9 0.2; 1.1 -0.1]
    current_values = [begin
        delta_star = last(Magnetostatics._rbf_terms(R, Z, centers, length_scale))
        -dot(delta_star, true_field.psi_coefficients) / (Magnetostatics.μ₀ * R)
    end for (R, Z) in eachrow(current_positions)]

    reconstructed = reconstruct_axisymmetric_field(
        position_matrix, field_matrix;
        domain, length_scale, regularization = 1.0e-16,
        basis_centers = centers,
        boundary_positions, boundary_flux = boundary_values,
        current_positions, current_density = current_values,
        current_uncertainties = fill(1.0e3, length(current_values))
    )

    @test all(
        isapprox(reconstructed(SVector(R, 0.0, Z)), true_field(SVector(R, 0.0, Z));
            rtol = 1.0e-7, atol = 1.0e-9)
        for (R, Z) in positions
    )

    noisy_fields = copy(field_matrix)
    noisy_fields[1, 1] += 0.02
    weighted_uncertainties = fill(1.0e-4, size(field_matrix))
    weighted_uncertainties[1, 1] = 0.02
    weighted_fit = reconstruct_axisymmetric_field(
        position_matrix, noisy_fields;
        domain, length_scale, regularization = 1.0e-16,
        basis_centers = centers, uncertainties = weighted_uncertainties
    )
    unweighted_fit = reconstruct_axisymmetric_field(
        position_matrix, noisy_fields;
        domain, length_scale, regularization = 1.0e-16,
        basis_centers = centers, uncertainties = fill(1.0e-4, size(field_matrix))
    )
    heldout = SVector(0.77, 0.0, 0.17)
    @test norm(weighted_fit(heldout) - true_field(heldout)) <
        norm(unweighted_fit(heldout) - true_field(heldout))

    R, Z, h = 0.82, 0.13, 1.0e-5
    Bp = reconstructed(SVector(R + h, 0.0, Z))
    Bm = reconstructed(SVector(R - h, 0.0, Z))
    Bz_p = reconstructed(SVector(R, 0.0, Z + h))[3]
    Bz_m = reconstructed(SVector(R, 0.0, Z - h))[3]
    div_B = ((R + h) * Bp[1] - (R - h) * Bm[1]) / (2h * R) +
        (Bz_p - Bz_m) / (2h)
    @test abs(div_B) < 1.0e-8
    @test poloidal_flux(reconstructed, SVector(R, 0.0, Z)) ≈
        poloidal_flux(true_field, SVector(R, 0.0, Z))
    @test_throws DomainError reconstructed(SVector(1.4, 0.0, 0.0))
    @test_throws DomainError poloidal_flux(reconstructed, SVector(1.4, 0.0, 0.0))

    boundary_only = reconstruct_axisymmetric_field(
        position_matrix, fill(NaN, size(field_matrix));
        domain, length_scale, regularization = 1.0e-16,
        basis_centers = centers, boundary_positions,
        boundary_flux = boundary_values
    )
    @test all(isapprox(
        poloidal_flux(boundary_only, SVector(boundary_positions[i, 1], 0.0,
            boundary_positions[i, 2])), boundary_values[i]; atol = 1.0e-12
    ) for i in eachindex(boundary_values))

    soft_boundary = reconstruct_axisymmetric_field(
        position_matrix, field_matrix;
        domain, length_scale, regularization = 1.0e-16,
        basis_centers = centers, boundary_positions,
        boundary_flux = boundary_values,
        boundary_uncertainties = fill(1.0e-4, length(boundary_values))
    )
    @test isfinite(poloidal_flux(soft_boundary, SVector(R, 0.0, Z)))

    conflicting_boundary = repeat(boundary_positions[1:1, :], 2, 1)
    @test_throws ArgumentError reconstruct_axisymmetric_field(
        position_matrix, fill(NaN, size(field_matrix));
        domain, length_scale, regularization = 1.0e-16,
        basis_centers = centers,
        boundary_positions = conflicting_boundary,
        boundary_flux = [0.0, 1.0]
    )
end

@testset "Sparse FLARE field samples" begin
    flare = FLARECoils(; I_pf = 40.0e3, I_ef = 10.0e3, B_guide = 0.015)
    positions = [(R, Z) for R in range(0.45, 1.25; length = 6)
                          for Z in range(-0.55, 0.55; length = 5)]
    position_matrix = reduce(vcat, [reshape(collect(pos), 1, 2) for pos in positions])
    field_matrix = reduce(vcat, [reshape(collect(flare(SVector(R, 0.0, Z))), 1, 3)
                                 for (R, Z) in positions])
    reconstructed = reconstruct_axisymmetric_field(
        position_matrix, field_matrix;
        domain = (0.4, 1.3, -0.6, 0.6),
        length_scale = 0.3,
        regularization = 1.0e-12
    )

    sensor_errors = [
        norm(reconstructed(SVector(R, 0.0, Z)) - flare(SVector(R, 0.0, Z)))
        for (R, Z) in positions
    ]
    @test sqrt(sum(abs2, sensor_errors) / length(sensor_errors)) < 0.01
    @test maximum(sensor_errors) < 0.02
end

end # module test_reconstruction