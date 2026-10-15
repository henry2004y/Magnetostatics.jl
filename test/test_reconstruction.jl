module test_reconstruction

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra
using LinearSolve: QRFactorization

@testset "Axisymmetric field reconstruction" begin
    domain = (0.4, 1.3, -0.6, 0.6)
    centers = reduce(
        vcat, [
            reshape([R, Z], 1, 2) for R in (0.5, 0.85, 1.2)
                for Z in (-0.45, 0.0, 0.45)
        ]
    )
    length_scale = 0.3
    true_field = AxisymmetricFieldReconstruction(
        centers,
        1.0e-3 .* [1.0, -0.4, 0.8, 0.5, -1.2, 0.7, -0.6, 0.3, 1.1],
        0.02 .* [0.2, -0.3, 0.5, 0.7, 0.1, -0.4, 0.6, -0.2, 0.3],
        length_scale, domain, :gaussian,
        Magnetostatics._build_center_bins(centers, domain, length_scale)
    )
    positions = [
        (R, Z) for R in range(0.48, 1.22; length = 5)
            for Z in range(-0.5, 0.5; length = 4)
    ]
    position_matrix = reduce(vcat, [reshape(collect(pos), 1, 2) for pos in positions])
    field_matrix = reduce(
        vcat, [
            reshape(collect(true_field(SVector(R, 0.0, Z))), 1, 3)
                for (R, Z) in positions
        ]
    )
    boundary_positions = reduce(
        vcat, [
            reshape([R, Z], 1, 2)
                for R in (domain[1], domain[2]), Z in (-0.5, 0.0, 0.5)
        ]
    )
    boundary_values = [
        poloidal_flux(true_field, SVector(R, 0.0, Z))
            for (R, Z) in eachrow(boundary_positions)
    ]
    current_positions = [0.6 -0.3; 0.9 0.2; 1.1 -0.1]
    current_values = [
        begin
            delta_star = last(
                Magnetostatics._rbf_terms(
                    R, Z, centers, length_scale, :gaussian
                )
            )
            -dot(delta_star, true_field.psi_coefficients) / (Magnetostatics.μ₀ * R)
        end for (R, Z) in eachrow(current_positions)
    ]

    reconstructed = reconstruct_axisymmetric_field(
        position_matrix, field_matrix;
        domain, length_scale, regularization = 1.0e-16,
        basis_centers = centers, kernel = :gaussian,
        boundary_positions, boundary_flux = boundary_values,
        current_positions, current_density = current_values,
        current_uncertainties = fill(1.0e3, length(current_values))
    )

    @test all(
        isapprox(
            reconstructed(SVector(R, 0.0, Z)), true_field(SVector(R, 0.0, Z));
            rtol = 1.0e-7, atol = 1.0e-9
        )
            for (R, Z) in positions
    )

    noisy_fields = copy(field_matrix)
    noisy_fields[1, 1] += 0.02
    weighted_uncertainties = fill(1.0e-4, size(field_matrix))
    weighted_uncertainties[1, 1] = 0.02
    weighted_fit = reconstruct_axisymmetric_field(
        position_matrix, noisy_fields;
        domain, length_scale, regularization = 1.0e-16,
        basis_centers = centers, kernel = :gaussian,
        uncertainties = weighted_uncertainties
    )
    unweighted_fit = reconstruct_axisymmetric_field(
        position_matrix, noisy_fields;
        domain, length_scale, regularization = 1.0e-16,
        basis_centers = centers, kernel = :gaussian,
        uncertainties = fill(1.0e-4, size(field_matrix))
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
        basis_centers = centers, kernel = :gaussian, boundary_positions,
        boundary_flux = boundary_values
    )
    @test all(
        isapprox(
            poloidal_flux(
                boundary_only, SVector(
                    boundary_positions[i, 1], 0.0,
                    boundary_positions[i, 2]
                )
            ), boundary_values[i]; atol = 1.0e-12
        ) for i in eachindex(boundary_values)
    )

    soft_boundary = reconstruct_axisymmetric_field(
        position_matrix, field_matrix;
        domain, length_scale, regularization = 1.0e-16,
        basis_centers = centers, kernel = :gaussian, boundary_positions,
        boundary_flux = boundary_values,
        boundary_uncertainties = fill(1.0e-4, length(boundary_values))
    )
    @test isfinite(poloidal_flux(soft_boundary, SVector(R, 0.0, Z)))

    conflicting_boundary = repeat(boundary_positions[1:1, :], 2, 1)
    @test_throws ArgumentError reconstruct_axisymmetric_field(
        position_matrix, fill(NaN, size(field_matrix));
        domain, length_scale, regularization = 1.0e-16,
        basis_centers = centers, kernel = :gaussian,
        boundary_positions = conflicting_boundary,
        boundary_flux = [0.0, 1.0]
    )
end

@testset "Sparse FLARE field samples" begin
    flare = FLARECoils(; I_pf = 40.0e3, I_ef = 10.0e3, B_guide = 0.015)
    positions = [
        (R, Z) for R in range(0.45, 1.25; length = 6)
            for Z in range(-0.55, 0.55; length = 5)
    ]
    position_matrix = reduce(vcat, [reshape(collect(pos), 1, 2) for pos in positions])
    field_matrix = reduce(
        vcat, [
            reshape(collect(flare(SVector(R, 0.0, Z))), 1, 3)
                for (R, Z) in positions
        ]
    )
    reconstructed = reconstruct_axisymmetric_field(
        position_matrix, field_matrix;
        domain = (0.4, 1.3, -0.6, 0.6),
        length_scale = 0.3,
        regularization = 1.0e-12
    )
    support_values = first(
        Magnetostatics._rbf_terms(
            0.45, -0.55, reconstructed.centers, 0.3, :wendland_c2
        )
    )
    @test count(!iszero, support_values) < length(support_values)

    sensor_errors = [
        norm(reconstructed(SVector(R, 0.0, Z)) - flare(SVector(R, 0.0, Z)))
            for (R, Z) in positions
    ]
    @test sqrt(sum(abs2, sensor_errors) / length(sensor_errors)) < 0.01
    @test maximum(sensor_errors) < 0.02

    R, Z, h = 0.82, 0.13, 1.0e-5
    Bp = reconstructed(SVector(R + h, 0.0, Z))
    Bm = reconstructed(SVector(R - h, 0.0, Z))
    Bz_p = reconstructed(SVector(R, 0.0, Z + h))[3]
    Bz_m = reconstructed(SVector(R, 0.0, Z - h))[3]
    div_B = ((R + h) * Bp[1] - (R - h) * Bm[1]) / (2h * R) +
        (Bz_p - Bz_m) / (2h)
    @test abs(div_B) < 1.0e-6
end

@testset "Vacuum basis terms" begin
    let degree = 6, scale = 1.3, R = 0.83, Z = -0.27, h = 1.0e-4
        ψ(R, Z) = first(Magnetostatics._vacuum_terms(R, Z, degree, scale))
        values, dR, dZ = Magnetostatics._vacuum_terms(R, Z, degree, scale)
        @test length(values) == degree
        @test dR ≈ (ψ(R + h, Z) - ψ(R - h, Z)) / 2h rtol = 1.0e-6
        @test dZ ≈ (ψ(R, Z + h) - ψ(R, Z - h)) / 2h rtol = 1.0e-6 atol = 1.0e-10
        delta_star = (ψ(R + h, Z) - 2ψ(R, Z) + ψ(R - h, Z)) / h^2 -
            (ψ(R + h, Z) - ψ(R - h, Z)) / (2h * R) +
            (ψ(R, Z + h) - 2ψ(R, Z) + ψ(R, Z - h)) / h^2
        @test maximum(abs, delta_star) < 1.0e-5
    end
    @test isempty(first(Magnetostatics._vacuum_terms(1.0, 0.0, -1, 1.0)))
end

function rectangle_points(domain, offset, n)
    Rmin, Rmax = domain[1] - offset, domain[2] + offset
    Zmin, Zmax = domain[3] - offset, domain[4] + offset
    perimeter = 2 * ((Rmax - Rmin) + (Zmax - Zmin))
    points = zeros(n, 2)
    for (i, s) in enumerate(range(0, perimeter; length = n + 1)[1:n])
        if s < Rmax - Rmin
            points[i, :] .= (Rmin + s, Zmin)
        elseif s < (Rmax - Rmin) + (Zmax - Zmin)
            points[i, :] .= (Rmax, Zmin + s - (Rmax - Rmin))
        elseif s < 2 * (Rmax - Rmin) + (Zmax - Zmin)
            points[i, :] .= (Rmax - (s - (Rmax - Rmin) - (Zmax - Zmin)), Zmax)
        else
            points[i, :] .= (Rmin, Zmax - (s - 2 * (Rmax - Rmin) - (Zmax - Zmin)))
        end
    end
    return points
end

@testset "Coil-free FLARE reconstruction" begin
    flare = FLARECoils(; I_pf = 40.0e3, I_ef = 10.0e3, B_guide = 0.015)
    domain = (0.45, 1.25, -0.28, 0.28)
    positions = [
        (R, Z) for R in range(0.48, 1.22; length = 10)
            for Z in range(-0.25, 0.25; length = 7)
    ]
    position_matrix = reduce(vcat, [reshape(collect(pos), 1, 2) for pos in positions])
    field_matrix = reduce(
        vcat, [
            reshape(collect(flare(SVector(R, 0.0, Z))), 1, 3)
                for (R, Z) in positions
        ]
    )
    filaments = rectangle_points(domain, 0.1, 48)
    vacuum_positions = reduce(
        vcat, [
            reshape([R, Z], 1, 2) for R in range(domain[1], domain[2]; length = 9)
                for Z in range(domain[3], domain[4]; length = 6)
        ]
    )
    common = (;
        domain, length_scale = 0.3, regularization = 1.0e3,
        uncertainties = fill(5.0e-4, size(field_matrix)),
        current_positions = vacuum_positions,
        current_density = zeros(size(vacuum_positions, 1)),
        current_uncertainties = fill(1.0e4, size(vacuum_positions, 1)),
        algorithm = QRFactorization(),
    )
    reconstructed = reconstruct_axisymmetric_field(
        position_matrix, field_matrix; common..., filament_positions = filaments
    )
    rbf_only = reconstruct_axisymmetric_field(position_matrix, field_matrix; common...)

    @test length(reconstructed.filament_currents) == size(filaments, 1)
    heldout = [
        SVector(R, 0.0, Z) for R in range(0.5, 1.2; length = 13)
            for Z in range(-0.22, 0.22; length = 9)
    ]
    poloidal_error(field) = sqrt(
        sum(p -> sum(abs2, (field(p) - flare(p))[SA[1, 3]]), heldout) /
            sum(p -> sum(abs2, flare(p)[SA[1, 3]]), heldout)
    )
    @test poloidal_error(reconstructed) < 0.02
    @test poloidal_error(reconstructed) < poloidal_error(rbf_only)

    with_vacuum = reconstruct_axisymmetric_field(
        position_matrix, field_matrix; common..., vacuum_degree = 6
    )
    @test length(with_vacuum.vacuum_coefficients) == 6
    @test poloidal_error(with_vacuum) < poloidal_error(rbf_only)

    @test_throws ArgumentError reconstruct_axisymmetric_field(
        position_matrix, field_matrix; common...,
        filament_positions = [0.8 0.0]
    )
    @test_throws ArgumentError reconstruct_axisymmetric_field(
        position_matrix, field_matrix; common..., vacuum_degree = -1
    )
end

@testset "Current density uncertainties" begin
    domain = (0.4, 1.3, -0.6, 0.6)
    positions = [0.6 0.0; 1.0 0.0; 0.8 0.3; 0.8 -0.3]
    field_matrix = zeros(size(positions, 1), 3)
    current_positions = [0.85 0.0]
    centers = reduce(
        vcat, [reshape([R, Z], 1, 2) for R in (0.6, 0.85, 1.1) for Z in (-0.3, 0.0, 0.3)]
    )
    J_target = [1.0e5]
    fit(σ) = reconstruct_axisymmetric_field(
        positions, field_matrix;
        domain, length_scale = 0.3, regularization = 1.0e-12,
        basis_centers = centers,
        uncertainties = fill(1.0e-2, size(field_matrix)),
        current_positions, current_density = J_target,
        current_uncertainties = [σ]
    )
    J_phi(field) = current_density(field, SVector(0.85, 0.0, 0.0))[2]
    @test J_phi(fit(1.0)) ≈ J_target[1] rtol = 1.0e-3
    @test abs(J_phi(fit(1.0e9))) < abs(J_phi(fit(1.0)))

    field_matrix[:, 2] .= (0.01, 0.02, -0.01, 0.005)
    field = reconstruct_axisymmetric_field(
        positions, field_matrix;
        domain, length_scale = 0.3, basis_centers = centers,
        current_positions, current_density = J_target
    )
    p = SVector(0.8 * cos(0.3), 0.8 * sin(0.3), 0.1)
    h = 1.0e-5
    unit_vectors = (SA[1.0, 0, 0], SA[0, 1.0, 0], SA[0, 0, 1.0])
    dB = [(field(p + h * e) - field(p - h * e)) / 2h for e in unit_vectors]
    curl = SVector(dB[2][3] - dB[3][2], dB[3][1] - dB[1][3], dB[1][2] - dB[2][1])
    J = current_density(field, p)
    @test abs(J[3]) > 0
    @test J ≈ curl / Magnetostatics.μ₀ rtol = 1.0e-5
    @test_throws DomainError current_density(field, SVector(2.0, 0.0, 0.0))
end

struct GaussianCoreSheet
    loops::Vector{CurrentLoop{Float64}}
    σ::Float64
end

function GaussianCoreSheet(; Rc, a, δ, Ip, h = 0.01, σ = 0.01)
    nodes = [
        (R + 0.00317, Z + 0.00213) for R in (Rc - 3a):h:(Rc + 3a) for Z in (-4δ):h:(4δ)
    ]
    w = [exp(-((R - Rc) / a)^2) * sech(Z / δ)^2 for (R, Z) in nodes]
    loops = [
        CurrentLoop(R, Ip * wk / sum(w), SVector(0.0, 0.0, Z), SVector(0.0, 0.0, 1.0))
            for ((R, Z), wk) in zip(nodes, w)
    ]
    return GaussianCoreSheet(loops, σ)
end

function (s::GaussianCoreSheet)(p)
    R = hypot(p[1], p[2])
    B = zero(SVector{3, Float64})
    for loop in s.loops
        B += getB_loop(p, loop)
        dR, dZ = R - loop.radius, p[3] - loop.center[3]
        d2 = dR^2 + dZ^2
        c = -Magnetostatics.μ₀ * loop.current * exp(-d2 / (2s.σ^2)) / (2π * d2)
        B += SVector(c * dZ * p[1] / R, c * dZ * p[2] / R, -c * dR)
    end
    return B
end

function toroidal_current(s::GaussianCoreSheet, R, Z)
    return sum(s.loops) do loop
        d2 = (R - loop.radius)^2 + (Z - loop.center[3])^2
        loop.current * exp(-d2 / (2s.σ^2))
    end / (2π * s.σ^2)
end

@testset "FLARE plasma current sheet reconstruction" begin
    flare = FLARECoils(; I_pf = 45.0e3, I_ef = 12.0e3, B_guide = 0.0)
    sheet = GaussianCoreSheet(; Rc = 0.95, a = 0.15, δ = 0.05, Ip = 15.0e3)
    reference(p) = flare(p) + sheet(p)

    let R = 0.9, Z = 0.02, h = 1.0e-5
        p(R, Z) = SVector(R, 0.0, Z)
        dBR_dZ = (reference(p(R, Z + h))[1] - reference(p(R, Z - h))[1]) / 2h
        dBZ_dR = (reference(p(R + h, Z))[3] - reference(p(R - h, Z))[3]) / 2h
        J_curl = (dBR_dZ - dBZ_dR) / Magnetostatics.μ₀
        @test J_curl ≈ toroidal_current(sheet, R, Z) rtol = 1.0e-4
    end

    domain = (0.45, 1.25, -0.28, 0.28)
    positions = [
        (R, Z) for R in range(0.48, 1.22; length = 12)
            for Z in range(-0.25, 0.25; length = 11)
    ]
    position_matrix = reduce(vcat, [reshape(collect(pos), 1, 2) for pos in positions])
    field_matrix = reduce(
        vcat, [
            reshape(collect(reference(SVector(R, 0.0, Z))), 1, 3)
                for (R, Z) in positions
        ]
    )
    for i in axes(field_matrix, 1)
        field_matrix[i, 1] += 2.0e-4 * sin(i)
        field_matrix[i, 3] += 2.0e-4 * sin(2i)
    end
    ls = 0.06
    grid(h) = reduce(
        vcat, [
            reshape([R, Z], 1, 2) for R in domain[1]:h:domain[2]
                for Z in domain[3]:h:domain[4]
        ]
    )
    collocation = grid(ls / 4)
    reconstructed = reconstruct_axisymmetric_field(
        position_matrix, field_matrix;
        domain, kernel = :gaussian, length_scale = ls,
        basis_centers = grid(ls / 2), regularization = 1.0e2,
        uncertainties = fill(5.0e-4, size(field_matrix)),
        filament_positions = rectangle_points(domain, 0.1, 48),
        current_positions = collocation,
        current_density = zeros(size(collocation, 1)),
        current_uncertainties = fill(1.0e5, size(collocation, 1)),
        algorithm = QRFactorization()
    )

    Rs = range(domain[1], domain[2]; length = 41)
    Zs = range(domain[3], domain[4]; length = 29)
    heldout = [SVector(R, 0.0, Z) for R in Rs for Z in Zs]
    Bp_error = sqrt(
        sum(p -> sum(abs2, (reconstructed(p) - reference(p))[SA[1, 3]]), heldout) /
            sum(p -> sum(abs2, reference(p)[SA[1, 3]]), heldout)
    )
    @test Bp_error < 0.02

    J_ref = [toroidal_current(sheet, p[1], p[3]) for p in heldout]
    J_rec = [current_density(reconstructed, p)[2] for p in heldout]
    dA = step(Rs) * step(Zs)
    @test sum(J_rec) * dA ≈ sum(J_ref) * dA rtol = 0.03
    @test norm(J_rec - J_ref) / norm(J_ref) < 0.25

    with_elements = reconstruct_axisymmetric_field(
        position_matrix, field_matrix;
        domain, length_scale = 0.12, basis_centers = grid(0.06),
        regularization = 1.0e2,
        uncertainties = fill(5.0e-4, size(field_matrix)),
        filament_positions = rectangle_points(domain, 0.1, 48),
        current_element_positions = grid(0.05),
        current_positions = collocation,
        current_density = zeros(size(collocation, 1)),
        current_uncertainties = fill(1.0e5, size(collocation, 1))
    )
    @test with_elements.kernel === :wendland_c2
    @test with_elements.element_width ≈ 0.05
    @test sqrt(
        sum(p -> sum(abs2, (with_elements(p) - reference(p))[SA[1, 3]]), heldout) /
            sum(p -> sum(abs2, reference(p)[SA[1, 3]]), heldout)
    ) < 0.02
    J_elements = [current_density(with_elements, p)[2] for p in heldout]
    @test sum(J_elements) * dA ≈ sum(J_ref) * dA rtol = 0.03
    @test sum(with_elements.element_currents) ≈ 15.0e3 rtol = 0.05
    @test norm(J_elements - J_ref) / norm(J_ref) < 0.25
end

@testset "Current element kernel" begin
    let a = 0.9, R = 1.05, z = 0.13, h = 1.0e-6
        μ₀ = Magnetostatics.μ₀
        loop = CurrentLoop(a, inv(μ₀), SVector(0.0, 0, 0), SVector(0.0, 0, 1))
        psi, B_ρ, B_z, dB_ρ = Magnetostatics._unit_loop_terms(a, R, z)
        @test psi ≈ R * solve(VectorPotential(), loop, SVector(R, 0, z))[2] rtol = 1.0e-10
        @test SVector(B_ρ, 0, B_z) ≈ getB_loop(SVector(R, 0, z), loop) rtol = 1.0e-10
        B(z) = getB_loop(SVector(R, 0, z), loop)[1]
        @test dB_ρ ≈ (B(z + h) - B(z - h)) / 2h rtol = 1.0e-6
    end

    elements = [0.9 0.02]
    width = 0.04
    terms(R, Z) = first.(Magnetostatics._element_terms(R, Z, elements, width))
    for (R, Z) in ((1.0, 0.1), (0.9, 0.02), (0.92, 0.021), (0.6, -0.2))
        psi, dR, dZ, delta = terms(R, Z)
        ψ(R, Z) = terms(R, Z)[1]
        h = 1.0e-5
        @test dR ≈ (ψ(R + h, Z) - ψ(R - h, Z)) / 2h rtol = 1.0e-6
        @test dZ ≈ (ψ(R, Z + h) - ψ(R, Z - h)) / 2h atol = 1.0e-6
        h = 1.0e-4
        Δstar = (ψ(R + h, Z) - 2psi + ψ(R - h, Z)) / h^2 - dR / R +
            (ψ(R, Z + h) - 2psi + ψ(R, Z - h)) / h^2
        @test delta ≈ Δstar rtol = 1.0e-4
    end

    Rs = range(0.3, 1.5; length = 481)
    Zs = range(-0.6, 0.6; length = 481)
    μ₀J = [-terms(R, Z)[4] / R for R in Rs, Z in Zs]
    @test sum(μ₀J) * step(Rs) * step(Zs) ≈ 1 rtol = 0.01
    @test maximum(μ₀J) ≈ 1 / (π * width^2) rtol = 0.01

    @test Magnetostatics._median_spacing([0.5 0.0; 0.6 0.0; 0.8 0.0], 1.0) ≈ 0.1
    @test Magnetostatics._median_spacing([0.5 0.0], 1.0) == 1.0
    positions = [0.6 0.0; 1.0 0.0; 0.8 0.3; 0.8 -0.3]
    kw = (; domain = (0.4, 1.3, -0.6, 0.6), length_scale = 0.3)
    @test_throws ArgumentError reconstruct_axisymmetric_field(
        positions, zeros(4, 3); kw..., current_element_positions = [0.0 0.0]
    )
    @test_throws ArgumentError reconstruct_axisymmetric_field(
        positions, zeros(4, 3); kw...,
        current_element_positions = [0.8 0.0], current_element_width = -1.0
    )
end

end # module test_reconstruction
