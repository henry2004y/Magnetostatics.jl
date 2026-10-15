module test_fieldtracer

using Test
using FieldTracer
using Magnetostatics
using StaticArrays
using LinearAlgebra

@testset "Parker spiral callback" begin
    AU = 1.496e11
    r0 = 2.5 * 6.96e8
    Omega, v_sw = 2.865e-6, 4.0e5
    spiral = ParkerSpiral(5.0e-9 * (AU / r0)^2, r0; v_sw, Omega)
    outside(p, s) = norm(p) > 5.0 * AU
    ds = 0.002 * AU
    for latitude in (-60, -30, 0, 30, 60)
        theta, phi = deg2rad(90 - latitude), 0.7
        seed = r0 * SVector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
        result = trace(
            (p, s) -> spiral(p), collect(seed);
            ds, max_length = 20.0 * AU, direction = "forward", isoutside = outside
        )
        @test result.status == :outside
        @test first(result.points) == seed
        @test length(result.points) == length(result.s) > 1
        @test all(diff(result.s) .> 0)
        @test 0 <= 5.0 * AU - norm(last(result.points)) <= ds
        exact_point(p) = let r = norm(p)
            angle = phi - Omega / v_sw * (r - r0 - r0 * log(r / r0))
            r * SVector(sin(theta) * cos(angle), sin(theta) * sin(angle), cos(theta))
        end
        error = maximum(norm(p - exact_point(p)) / AU for p in result.points)
        @test error < 1.0e-5
        refined = trace(
            (p, s) -> spiral(p), collect(seed); ds = ds / 2,
            max_length = 20.0 * AU, direction = "forward", isoutside = outside
        )
        @test refined.status == :outside
        refined_error = maximum(norm(p - exact_point(p)) / AU for p in refined.points)
        @test refined_error < error / 4
    end
end

@testset "Reconstructed field stays inside its domain" begin
    domain = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
    positions = [-0.3 -0.2 0.1; 0.1 0.3 -0.2; 0.3 -0.1 0.2; -0.1 0.2 -0.3]
    truth = UniformField(SVector(0.01, -0.02, 0.03))
    samples = [truth(p)[k] for p in eachrow(positions), k in 1:3]
    reconstructed = reconstruct_field_3d(
        positions, samples; domain, length_scale = 2.0,
        basis_centers = zeros(1, 3), regularization = 1.0e-14
    )
    isoutside(p, s) = any(p[k] < domain[2k - 1] || p[k] > domain[2k] for k in 1:3)
    evaluated = Vector{Float64}[]
    field(p, s) = begin
        push!(evaluated, copy(p))
        reconstructed(p)
    end
    seed = SVector(0.0, 0.0, 0.0)
    ds = 0.01
    result = trace(field, collect(seed); ds, max_length = 4.0, isoutside)
    for (branch, direction) in ((result.forward, 1), (result.backward, -1))
        @test branch.status == :outside
        @test first(branch.points) == seed
        @test all(direction .* branch.s .>= 0)
        @test all(!isoutside(p, s) for (p, s) in zip(branch.points, branch.s))
        @test all(
            isapprox(p, seed + s * normalize(truth(seed)); atol = 1.0e-9)
                for (p, s) in zip(branch.points, branch.s)
        )
        endpoint = last(branch.points)
        @test minimum(
            abs(endpoint[k] - domain[j]) for k in 1:3 for j in (2k - 1, 2k)
        ) <= ds
    end
    @test !isempty(evaluated)
    @test all(!isoutside(p, 0.0) for p in evaluated)
    path = vcat(reverse(result.backward.points), result.forward.points[2:end])
    @test count(==(seed), path) == 1
    @test dot(last(path) - first(path), truth(seed)) > 0
    @test trace((p, s) -> NullField()(p), collect(seed)).forward.status == :zero_field
end

@testset "FFT component arrays" begin
    n, dx = 8, 0.1
    x = LinRange(-n * dx / 2, n * dx / 2 - dx, n)
    y, z = x, x
    J = zeros(3, n, n, n)
    set_current_wire!(
        J, x, y, z, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0), 1.0, 2.5dx
    )
    B = solve(FFTSolver(), J, dx)
    Bx, By, Bz = eachslice(B; dims = 1)
    ds = dx / 20
    for seed_z in (-0.2, 0.0, 0.2), radius in (0.08, 0.16, 0.24)
        maxstep = ceil(Int, 2pi * radius / ds) + 1
        path = trace(
            Bx, By, Bz, radius, 0.0, seed_z, x, y, z;
            ds, maxstep, gridtype = "ndgrid", direction = "forward"
        )
        @test path isa NTuple{3, Vector{Float64}}
        @test all(length(coordinates) == length(first(path)) for coordinates in path)
        @test 2 < length(first(path)) <= maxstep
        @test all(all(isfinite, coordinates) for coordinates in path)
        @test all(
            all(first(grid) <= p <= last(grid) for p in coordinates)
                for (grid, coordinates) in zip((x, y, z), path)
        )
        @test [first(coordinates) for coordinates in path] ≈ [radius, 0.0, seed_z]
        @test all(isapprox(p, seed_z; atol = 1.0e-12) for p in path[3])
        @test path[2][2] > path[2][1]
    end
end

end
