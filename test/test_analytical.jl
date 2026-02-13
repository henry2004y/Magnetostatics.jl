using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

@testset "Analytical Fields" begin
    @testset "Harris Sheet" begin
        B0 = 1.0
        L = 2.0
        harris = HarrisSheet(B0, L)

        @test harris(SVector(0.0, 0.0, 0.0)) == SVector(0.0, 0.0, 0.0)
        @test isapprox(harris(SVector(0.0, 0.0, 100.0))[1], B0, rtol = 1.0e-5)
        @test isapprox(harris(SVector(0.0, 0.0, -100.0))[1], -B0, rtol = 1.0e-5)
    end

    @testset "Dipole" begin
        M = SVector(0.0, 0.0, 1.0)
        dipole = Dipole(M)

        # On z-axis, B = (mu0/4pi) * 2M / z^3
        z = 2.0
        B_z = dipole(SVector(0.0, 0.0, z))
        μ0_4π = 1.0e-7
        expected = μ0_4π * 2 * 1.0 / z^3

        @test isapprox(B_z[3], expected, rtol = 1.0e-5)

        # On x-axis, B = -(mu0/4pi) * M / x^3
        x = 2.0
        B_x = dipole(SVector(x, 0.0, 0.0))
        expected_x = -μ0_4π * 1.0 / x^3

        @test isapprox(B_x[3], expected_x, rtol = 1.0e-5)
    end

    @testset "CurrentLoop" begin
        # Test on-axis field strength
        R = 1.0
        I = 1.0
        center = SVector(0.0, 0.0, 0.0)
        normal = SVector(0.0, 0.0, 1.0)
        loop = CurrentLoop(R, I, center, normal)

        # At center
        B_center = getB_loop(center, loop)
        μ0 = 4π * 1.0e-7
        expected_center = μ0 * I / (2 * R)
        @test isapprox(B_center[3], expected_center, rtol = 1.0e-5)
        @test abs(B_center[1]) < 1.0e-10
        @test abs(B_center[2]) < 1.0e-10

        # At z = R
        z = R
        r = SVector(0.0, 0.0, z)
        B_z = getB_loop(r, loop)
        expected_z = μ0 * I * R^2 / (2 * (R^2 + z^2)^(3 / 2))
        @test isapprox(B_z[3], expected_z, rtol = 1.0e-5)

        # Test off-axis (symmetry check)
        r_pos = SVector(0.5, 0.0, 0.5)
        r_neg = SVector(-0.5, 0.0, 0.5)
        B_pos = getB_loop(r_pos, loop)
        B_neg = getB_loop(r_neg, loop)

        @test isapprox(B_pos[1], -B_neg[1], rtol = 1.0e-10) # B_rho component reverses in x
        @test isapprox(B_pos[3], B_neg[3], rtol = 1.0e-10)  # B_z component symmetric
    end
end
