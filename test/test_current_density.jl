module test_current_density
using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

@testset "Current Density" begin
    let B0 = 1.0, L = 1.0, d = 2.0
        sheet = BifurcatedHarrisSheet(B0, L, d)
        r = SA[0.0, 0.0, 0.0]
        j = current_density(sheet, r)

        # Expected value at z=0: B0 / (2 * μ₀ * L) * 2 * sech(d / L)^2 = B0 / (μ₀ * L) * sech(d)^2
        # μ₀ is approximately 1.2566e-6
        expected_jy = B0 / (Magnetostatics.μ₀ * L) * sech(d / L)^2
        @test j[2] ≈ expected_jy
        @test j[1] == 0.0
        @test j[3] == 0.0
    end

    let B0 = 1.2, L = 0.5
        sheet = HarrisSheet(B0, L)
        r = SA[0.0, 0.0, 0.5]
        j = current_density(sheet, r)
        expected_jy = B0 / (Magnetostatics.μ₀ * L) * sech(r[3] / L)^2
        @test j[2] ≈ expected_jy
    end

    let B1 = 2.0, B2 = 1.0, L = 1.5
        sheet = AsymmetricHarrisSheet(B1, B2, L)
        r = SA[0.0, 0.0, -1.0]
        j = current_density(sheet, r)
        expected_jy = (B1 - B2) / (2 * Magnetostatics.μ₀ * L) * sech(r[3] / L)^2
        @test j[2] ≈ expected_jy
    end

    let B0 = 1.0, L = 2.0
        sheet = ForceFreeHarrisSheet(B0, L)
        r = SA[0.0, 0.0, 1.0]
        j = current_density(sheet, r)
        expected_jx = B0 / (Magnetostatics.μ₀ * L) * sech(r[3] / L) * tanh(r[3] / L)
        expected_jy = B0 / (Magnetostatics.μ₀ * L) * sech(r[3] / L)^2
        @test j[1] ≈ expected_jx
        @test j[2] ≈ expected_jy
    end

    let B0 = 1.0, L = 1.0, Lx = 1.0, ε = 0.5
        fadeev = FadeevIsland(B0, L, Lx, ε)
        r = SA[0.0, 0.0, 0.0] # Center of an island (or between)
        j = current_density(fadeev, r)

        # At r=0: f = cosh(0) + ε cos(0) = 1 + ε
        # term1 = (1/L) * (1 + ε)
        # term2 = (L ε / Lx^2) * (ε + 1)
        # For L=Lx=1: term1 = term2, so j_y = B0 / (μ₀ (1+ε)^2) * (1 + ε - ε(1 + ε)/1) ? No.
        # Wait. term1 = 1 + ε. term2 = ε(ε + 1).
        # j_y = B0 / (μ₀ (1+ε)^2) * (1 + ε - ε^2 - ε) = B0 / (μ₀ (1+ε)^2) * (1 - ε^2)
        # = B0 (1-ε) (1+ε) / (μ₀ (1+ε)^2) = B0 (1-ε) / (μ₀ (1+ε))

        expected_jy = B0 * (1 - ε) / (Magnetostatics.μ₀ * (1 + ε))
        @test j[2] ≈ expected_jy
    end
end

end # module test_current_density
