module test_equilibria

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra
import Magnetostatics: besselj0, besselj1


@testset "Plasma Equilibria" begin
    @testset "SolovevEquilibrium" begin
        B0 = 2.0
        R0 = 3.0
        kappa = 1.5
        q0 = 1.1

        eq = SolovevEquilibrium(B0, R0, kappa, q0)

        # On-axis: R = R0, Z = 0, phi = 0 -> (R0, 0, 0)
        r_axis = SVector(R0, 0.0, 0.0)
        @test poloidal_flux(eq, r_axis) ≈ 0.0
        B_axis = eq(r_axis)
        # B_R = 0, B_Z = 0, B_phi = B0. At x=R0, y=0: B_x = 0, B_y = B0, B_z = 0
        @test isapprox(B_axis, SVector(0.0, B0, 0.0), atol = 1.0e-10)

        # Flux should be positive away from axis
        r_off = SVector(R0 + 0.5, 0.0, 0.5)
        @test poloidal_flux(eq, r_off) > 0.0

        # In-place evaluation
        dr = zeros(3)
        eq(dr, r_off)
        @test dr ≈ eq(r_off)

        # Origin handling
        @test eq(SVector(0.0, 0.0, 0.0)) == SVector(0.0, 0.0, 0.0)
    end

    @testset "TaylorState" begin
        B0 = 1.0
        alpha = 2.0
        taylor = TaylorState(B0, alpha)

        # On-axis: B = (0, 0, B0)
        r_axis = SVector(0.0, 0.0, 1.0)
        @test isapprox(taylor(r_axis), SVector(0.0, 0.0, B0), atol = 1.0e-12)

        # Off-axis point
        r_pt = SVector(0.5, 0.0, 0.0)
        B_pt = taylor(r_pt)
        # J1(alpha * r) for By, J0(alpha * r) for Bz
        @test B_pt[1] ≈ 0.0 atol = 1.0e-12
        @test B_pt[2] ≈ B0 * besselj1(alpha * 0.5)
        @test B_pt[3] ≈ B0 * besselj0(alpha * 0.5)

        # Current density: J = alpha / mu0 * B
        J_pt = current_density(taylor, r_pt)
        @test J_pt ≈ (alpha / Magnetostatics.μ₀) * B_pt

        # Vector potential: A = 1 / alpha * B
        A_pt = vector_potential(taylor, r_pt)
        @test A_pt ≈ (1.0 / alpha) * B_pt

        # RFP reversal point: J0(x0) = 0 where x0 ≈ 2.4048255577
        r_rev = 2.4048255577 / alpha
        B_rev = taylor(SVector(r_rev, 0.0, 0.0))
        @test abs(B_rev[3]) < 1.0e-6
    end

    @testset "ScrewPinch" begin
        I_val = 50.0
        a_val = 0.5
        Bz0 = 1.5
        pinch = ScrewPinch(I_val, a_val, Bz0)

        # On axis
        @test pinch(SVector(0.0, 0.0, 0.0)) ≈ SVector(0.0, 0.0, Bz0)

        # Inside pinch
        r_in = SVector(0.2, 0.0, 0.0)
        B_in = pinch(r_in)
        expected_By = Magnetostatics.μ₀ * I_val * 0.2 / (2π * a_val^2)
        @test isapprox(B_in[2], expected_By, rtol = 1.0e-5)
        @test B_in[3] == Bz0

        # Current density
        J_in = current_density(pinch, r_in)
        @test J_in ≈ SVector(0.0, 0.0, I_val / (π * a_val^2))
        J_out = current_density(pinch, SVector(1.0, 0.0, 0.0))
        @test J_out == SVector(0.0, 0.0, 0.0)
    end
end

end # module test_equilibria
