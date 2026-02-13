using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

@testset "BiotSavart Solver" begin
    # Test 1: Field from a straight wire
    # B = mu0 * I / (2 * pi * r)

    # Create a long wire along z-axis
    L = 100.0
    points = [SVector(0.0, 0.0, z) for z in range(-L, L, length = 1000)]
    wire = Wire(points, 1.0)
    solver = BiotSavart()

    r = SVector(1.0, 0.0, 0.0) # 1m away
    B = solve(solver, wire, r)

    # Analytical expected value for infinite wire
    μ0 = 4π * 1.0e-7
    B_expected = μ0 * 1.0 / (2 * π * 1.0)

    # The wire is finite, so result will be slightly less
    # exact finite wire: B = mu0 I / (4 pi r) * (cos theta1 + cos theta2)
    # theta1 approx 0, theta2 approx 0 -> cos terms approx 2 => infinite wire result

    @test isapprox(B[2], B_expected, rtol = 1.0e-2)
    @test abs(B[1]) < 1.0e-10
    @test abs(B[3]) < 1.0e-10

    # Test 2: Field from a circular loop (on-axis)
    R = 1.0
    I = 1.0
    loop_wire = discretize_loop(R, 100, I) # 100 segments

    z = 0.5
    r_axis = SVector(0.0, 0.0, z)
    B_loop = solve(solver, loop_wire, r_axis)

    B_analytical = μ0 * I * R^2 / (2 * (R^2 + z^2)^(3 / 2))

    @test isapprox(B_loop[3], B_analytical, rtol = 1.0e-2)
    @test abs(B_loop[1]) < 1.0e-10
    @test abs(B_loop[2]) < 1.0e-10
end
