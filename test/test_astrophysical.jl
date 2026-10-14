module test_astrophysical

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

@testset "Astrophysical Magnetic Fields" begin
    @testset "ParkerSpiral" begin
        B0 = 5.0e-9 # 5 nT at 1 AU
        r0 = 1.496e11 # 1 AU in meters
        v_sw = 4.0e5
        Omega = 2.865e-6
        spiral = ParkerSpiral(B0, r0; v_sw = v_sw, Omega = Omega)

        # At reference radius on the equator (x = r0, y = 0, z = 0)
        # r = r0 -> B_phi = 0, Br = B0
        r_ref = SVector(r0, 0.0, 0.0)
        B_ref = spiral(r_ref)
        @test isapprox(B_ref, SVector(B0, 0.0, 0.0), rtol = 1.0e-5)

        # At 2 AU on the equator (x = 2*r0, y = 0, z = 0)
        r_2au = SVector(2 * r0, 0.0, 0.0)
        B_2au = spiral(r_2au)
        # Br should decrease as (r0/r)^2 = 1/4
        @test isapprox(B_2au[1], B0 / 4, rtol = 1.0e-5)
        # Bphi is in -y direction
        @test B_2au[2] < 0.0

        # In-place call
        dr = zeros(3)
        spiral(dr, r_2au)
        @test dr ≈ B_2au
    end

    @testset "GuideFieldHarrisSheet" begin
        B0 = 20.0e-9
        L = 100.0e3
        Bg = 10.0e-9
        Bn = 2.0e-9
        sheet = GuideFieldHarrisSheet(B0, L; Bg = Bg, Bn = Bn)

        r_mid = SVector(0.0, 0.0, 0.0)
        B_mid = sheet(r_mid)
        @test isapprox(B_mid, SVector(0.0, Bg, Bn), atol = 1.0e-15)

        # Far from sheet (z >> L)
        r_far = SVector(0.0, 0.0, 10 * L)
        B_far = sheet(r_far)
        @test isapprox(B_far[1], B0, rtol = 1.0e-5)
        @test B_far[2] == Bg
        @test B_far[3] == Bn

        # Current density at center: Jy = B0 / (mu0 * L)
        J_mid = current_density(sheet, r_mid)
        @test isapprox(J_mid[2], B0 / (Magnetostatics.μ₀ * L), rtol = 1.0e-5)

        # Vector potential
        A_mid = vector_potential(sheet, r_mid)
        @test length(A_mid) == 3
    end

    @testset "GoldHoyleFluxTube" begin
        B0 = 0.1
        tau = 0.5 # 1/m
        loop = GoldHoyleFluxTube(B0, tau)

        # On-axis: B = (0, 0, B0)
        @test loop(SVector(0.0, 0.0, 0.0)) ≈ SVector(0.0, 0.0, B0)

        # Off-axis point on x-axis (r = 1.0, x = 1.0, y = 0.0)
        r_val = 1.0
        B_pt = loop(SVector(r_val, 0.0, 0.0))
        # Bz = B0 / (1 + (tau*r)^2)
        expected_Bz = B0 / (1 + (tau * r_val)^2)
        expected_By = (B0 * tau * r_val) / (1 + (tau * r_val)^2)
        @test B_pt[3] ≈ expected_Bz
        @test B_pt[2] ≈ expected_By
        @test B_pt[1] ≈ 0.0 atol = 1.0e-12

        # Verify constant pitch angle: B_theta / Bz = tau * r
        @test (B_pt[2] / B_pt[3]) ≈ tau * r_val
    end
end

end # module test_astrophysical
