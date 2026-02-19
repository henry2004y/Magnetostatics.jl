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
        expected = Magnetostatics.μ0_4π * 2 * 1.0 / z^3

        @test isapprox(B_z[3], expected, rtol = 1.0e-5)

        # On x-axis, B = -(mu0/4pi) * M / x^3
        x = 2.0
        B_x = dipole(SVector(x, 0.0, 0.0))
        expected_x = -Magnetostatics.μ0_4π * 1.0 / x^3

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
        expected_center = Magnetostatics.μ₀ * I / (2 * R)
        @test isapprox(B_center[3], expected_center, rtol = 1.0e-5)
        @test abs(B_center[1]) < 1.0e-10 && abs(B_center[2]) < 1.0e-10

        # At z = R
        z = R
        r = SVector(0.0, 0.0, z)
        B_z = getB_loop(r, loop)
        expected_z = Magnetostatics.μ₀ * I * R^2 / (2 * (R^2 + z^2)^(3 / 2))
        @test isapprox(B_z[3], expected_z, rtol = 1.0e-5)

        # Test off-axis (symmetry check)
        r_pos = SVector(0.5, 0.0, 0.5)
        r_neg = SVector(-0.5, 0.0, 0.5)
        B_pos = getB_loop(r_pos, loop)
        B_neg = getB_loop(r_neg, loop)

        @test isapprox(B_pos[1], -B_neg[1], rtol = 1.0e-10) &&
            isapprox(B_pos[3], B_neg[3], rtol = 1.0e-10)
    end

    @testset "Configurations" begin
        @testset "getB_tokamak_coil" begin
            # User provided test
            # x, y, z, a, b, ICoils, IPlasma
            val = getB_tokamak_coil(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)[1]
            @test val ≈ -1.2473997957952395e-6
        end

        @testset "getB_tokamak_profile" begin
            q_profile(nr) = nr^2 + 2 * nr + 0.5
            val = getB_tokamak_profile(1.0, 1.0, 1.0, q_profile, 2.0, 1.0, 1.0)[1]
            @test val ≈ -0.7666260799054282
        end

        @testset "getB_mirror" begin
            # x, y, z, distance, a, I1
            B = getB_mirror(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
            @test B[3] ≈ 8.99176285573213e-7
        end

        @testset "getB_bottle" begin
            # x, y, z, distance, a, b, I1, I2
            B = getB_bottle(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0)
            @test B[3] ≈ 1.5274948162911718e-6
        end

        @testset "getB_zpinch" begin
            # Check B field at distance r from wire with current I
            # B = mu0 * I / (2 * pi * r)
            x = 1.0
            y = 0.0
            z = 0.0
            I = 1.0
            a = 0.1 # radius

            # Outside wire
            B = getB_zpinch(x, y, z, I, a)
            expected_mag = Magnetostatics.μ₀ * I / (2 * π * x)
            # Direction should be tangent to circle, here B_y
            @test isapprox(B[2], expected_mag, rtol = 1.0e-5)

            # Inside wire
            x_in = 0.05
            B_in = getB_zpinch(x_in, y, z, I, a)
            # B = mu0 * I * r / (2 * pi * a^2)
            expected_in_mag = Magnetostatics.μ₀ * I * x_in / (2 * π * a^2)
            @test isapprox(B_in[2], expected_in_mag, rtol = 1.0e-5)
        end
    end
    @testset "SVector Dispatch" begin
        r = SVector(1.0, 1.0, 1.0)

        # getB_mirror
        B_mirror_s = getB_mirror(r, 1.0, 1.0, 1.0)
        B_mirror_c = getB_mirror(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        @test B_mirror_s == B_mirror_c

        # getB_bottle
        B_bottle_s = getB_bottle(r, 1.0, 1.0, 1.0, 1.0, 1.0)
        B_bottle_c = getB_bottle(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        @test B_bottle_s == B_bottle_c

        # getB_tokamak_coil
        B_tok_coil_s = getB_tokamak_coil(r, 1.0, 1.0, 1.0, 1.0)
        B_tok_coil_c = getB_tokamak_coil(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        @test B_tok_coil_s == B_tok_coil_c

        # getB_tokamak_profile
        q_profile(nr) = nr^2 + 2 * nr + 0.5
        B_tok_prof_s = getB_tokamak_profile(r, q_profile, 2.0, 1.0, 1.0)
        B_tok_prof_c = getB_tokamak_profile(1.0, 1.0, 1.0, q_profile, 2.0, 1.0, 1.0)
        @test B_tok_prof_s == B_tok_prof_c

        # getB_zpinch
        B_zpinch_s = getB_zpinch(r, 1.0, 0.1)
        B_zpinch_c = getB_zpinch(1.0, 1.0, 1.0, 1.0, 0.1)
        @test B_zpinch_s == B_zpinch_c

        # HarrisSheet and Dipole with AbstractVector
        harris = HarrisSheet(1.0, 2.0)
        r_vec = [0.0, 0.0, 0.0]
        @test harris(r_vec) == SVector(0.0, 0.0, 0.0)

        dipole = Dipole(SVector(0.0, 0.0, 1.0))
        @test dipole(r_vec) == SVector(0.0, 0.0, 0.0)

        r_vec_z = [0.0, 0.0, 2.0]
        B_dipole = dipole(r_vec_z)
        @test B_dipole isa SVector
        @test isapprox(B_dipole[3], Magnetostatics.μ0_4π * 2 / 8, rtol = 1.0e-5)
    end
end
