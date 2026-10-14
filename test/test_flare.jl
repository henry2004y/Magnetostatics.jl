module test_flare

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

@testset "FLARE Experiment Coils" begin
    # Default FLARE configuration (vacuum reconnection field)
    flare = FLARECoils(;
        R_pf = 0.75,
        z_pf = 0.4,
        I_pf = 40.0e3,
        R_ef = 1.3,
        z_ef = 0.7,
        I_ef = 10.0e3,
        B_guide = 0.0
    )

    # On midplane (z = 0), radial field B_R must be 0 by up-down symmetry
    r_mid = SVector(0.75, 0.0, 0.0)
    B_mid = flare(r_mid)
    # Bx is along R on x-axis, so Bx should be zero
    @test abs(B_mid[1]) < 1.0e-12

    # In-place call
    dr = zeros(3)
    flare(dr, r_mid)
    @test dr ≈ B_mid

    # Guide field test: on x-axis, phi_hat = (0, 1, 0)
    Bg_val = 0.05 # 500 Gauss
    flare_gf = FLARECoils(;
        R_pf = 0.75,
        z_pf = 0.4,
        I_pf = 40.0e3,
        R_ef = 1.3,
        z_ef = 0.7,
        I_ef = 10.0e3,
        B_guide = Bg_val
    )
    B_gf = flare_gf(r_mid)
    @test isapprox(B_gf[2], Bg_val, rtol = 1.0e-5)
    @test abs(B_gf[1]) < 1.0e-12

    # Vector potential and poloidal flux
    A_mid = vector_potential(flare, r_mid)
    @test length(A_mid) == 3
    psi_mid = poloidal_flux(flare, r_mid)
    @test psi_mid isa Real

    # Plasma current layer addition
    flare_plasma = FLARECoils(;
        R_pf = 0.75,
        z_pf = 0.4,
        I_pf = 40.0e3,
        R_ef = 1.3,
        z_ef = 0.7,
        I_ef = 10.0e3,
        I_p = 5.0e3
    )
    B_plasma = flare_plasma(r_mid)
    # Plasma current alters Bz at the midplane
    @test B_plasma != B_mid
end

end # module test_flare
