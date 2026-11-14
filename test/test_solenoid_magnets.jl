module test_solenoid_magnets

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

function _fd_curl(f, r; h = 1.0e-5)
    ex = SVector(h, 0.0, 0.0)
    ey = SVector(0.0, h, 0.0)
    ez = SVector(0.0, 0.0, h)
    dF_dx = (f(r + ex) - f(r - ex)) / (2h)
    dF_dy = (f(r + ey) - f(r - ey)) / (2h)
    dF_dz = (f(r + ez) - f(r - ez)) / (2h)
    return SVector(
        dF_dy[3] - dF_dz[2],
        dF_dz[1] - dF_dx[3],
        dF_dx[2] - dF_dy[1]
    )
end

function _fd_div(f, r; h = 1.0e-5)
    ex = SVector(h, 0.0, 0.0)
    ey = SVector(0.0, h, 0.0)
    ez = SVector(0.0, 0.0, h)
    return (
        f(r + ex)[1] - f(r - ex)[1] +
            f(r + ey)[2] - f(r - ey)[2] +
            f(r + ez)[3] - f(r - ez)[3]
    ) / (2h)
end

function _zero_alloc_call(field, r)
    field(r)
    return @allocated field(r)
end

@testset "Solenoid, HelmholtzCoil and AntiHelmholtzCoil" begin
    let R = 0.5, L = 1.2, I_turn = 2.0, turns = 500
        sol = Solenoid(R, L, I_turn; turns = turns)
        NI = turns * I_turn
        @test sol.current ≈ NI

        # Compare against 2000 discrete CurrentLoops
        n_loops = 2000
        dz = L / n_loops
        loops = [
            CurrentLoop(
                    R, NI / n_loops, SVector(0.0, 0.0, -L / 2 + (i - 0.5) * dz),
                    SVector(0.0, 0.0, 1.0)
                ) for i in 1:n_loops
        ]
        ms = MultiSource(loops)
        for pt in (
                SVector(0.0, 0.0, 0.25), SVector(0.2, 0.0, 0.15),
                SVector(0.8, -0.3, 0.35), SVector(0.5, 0.0, 0.0),
            )
            # On-sheet point (0.5, 0, 0) skips B_z jump, test off-sheet points against stack
            if norm(SA[pt[1], pt[2]]) != R
                @test sol(pt) ≈ solve(BiotSavart(), ms, pt) rtol = 1.0e-5
                @test vector_potential(sol, pt) ≈
                    solve(VectorPotential(), ms, pt) rtol = 1.0e-5
            end
        end

        # Verify curl(A) == B and zero current outside sheet
        pt_off = SVector(0.22, -0.17, 0.31)
        @test _fd_curl(p -> vector_potential(sol, p), pt_off) ≈ sol(pt_off) rtol = 1.0e-5
        @test current_density(sol, pt_off) == SVector(0.0, 0.0, 0.0)
        @test poloidal_flux(sol, pt_off) ≈
            hypot(pt_off[1], pt_off[2]) *
            dot(vector_potential(sol, pt_off), normalize(SA[-pt_off[2], pt_off[1], 0.0]))
        @test _zero_alloc_call(sol, pt_off) == 0

        # Long solenoid limit
        sol_long = Solenoid(0.1, 40.0, 100.0)
        B0 = Magnetostatics.μ₀ * 100.0 / 40.0
        @test sol_long(SVector(0.05, 0.0, 0.0)) ≈ SVector(0.0, 0.0, B0) rtol = 1.0e-4

        # Translated and tilted solenoid
        c = SVector(0.3, -0.4, 0.2)
        n_tilde = normalize(SA[1.0, -2.0, 3.0])
        sol_tilt = Solenoid(R, L, NI; center = c, normal = n_tilde)
        pt_t = c + SVector(0.15, 0.25, -0.2)
        @test _fd_curl(p -> vector_potential(sol_tilt, p), pt_t) ≈
            sol_tilt(pt_t) rtol = 1.0e-4

        # Thick solenoid: curl(B) == mu0 * J inside winding & outside
        sol_thick = Solenoid(0.5, 1.0, 400.0; thickness = 0.2)
        pt_in_w = SVector(0.52, 0.0, 0.1)
        pt_core = SVector(0.15, 0.1, 0.1)
        @test _fd_curl(sol_thick, pt_in_w) ≈
            Magnetostatics.μ₀ * current_density(sol_thick, pt_in_w) rtol = 1.0e-3
        @test norm(_fd_curl(sol_thick, pt_core)) < 1.0e-11
        @test _fd_curl(p -> vector_potential(sol_thick, p), pt_in_w) ≈
            sol_thick(pt_in_w) rtol = 1.0e-4
        @test _fd_curl(p -> vector_potential(sol_thick, p), pt_core) ≈
            sol_thick(pt_core) rtol = 1.0e-4
        @test _zero_alloc_call(sol_thick, pt_in_w) == 0
    end

    let hcoil = HelmholtzCoil(1.0, 10.0)
        @test hcoil isa MagneticMirror
        B0 = hcoil(SVector(0.0, 0.0, 0.0))
        B_near = hcoil(SVector(0.0, 0.0, 0.02))
        @test isapprox(B0, B_near; rtol = 1.0e-6)
        p = SVector(0.3, -0.2, 0.1)
        @test _fd_curl(q -> vector_potential(hcoil, q), p) ≈ hcoil(p) rtol = 1.0e-4
    end

    let ah = AntiHelmholtzCoil(1.0, 10.0)
        @test ah.loop1.current == -ah.loop2.current
        @test ah(SVector(0.0, 0.0, 0.0)) == SVector(0.0, 0.0, 0.0)
        h = 1.0e-4
        dBz = (ah(SVector(0.0, 0.0, h))[3] - ah(SVector(0.0, 0.0, -h))[3]) / (2h)
        dBx = (ah(SVector(h, 0.0, 0.0))[1] - ah(SVector(-h, 0.0, 0.0))[1]) / (2h)
        @test dBz ≈ -2 * dBx rtol = 1.0e-6
        p = SVector(0.3, -0.2, 0.1)
        @test _fd_curl(q -> vector_potential(ah, q), p) ≈ ah(p) rtol = 1.0e-4
    end
end

@testset "Permanent Magnets: Sphere, Cylinder, Cuboid" begin
    let M = SVector(1.0e5, -2.0e5, 3.0e5), R = 0.4, c = SVector(0.2, -0.1, 0.3)
        sph = SphereMagnet(M, R, c)
        pt_in = c + SVector(0.1, -0.05, 0.08)
        @test sph(pt_in) ≈ (2 / 3) * Magnetostatics.μ₀ * M
        @test _fd_curl(q -> vector_potential(sph, q), pt_in) ≈ sph(pt_in) rtol = 1.0e-5
        pt_out = c + SVector(0.8, 0.3, -0.2)
        dip = Dipole((4 * π / 3) * R^3 * M, c)
        @test sph(pt_out) ≈ dip(pt_out)
        @test _fd_curl(q -> vector_potential(sph, q), pt_out) ≈ sph(pt_out) rtol = 1.0e-5
        # Continuity of normal B across surface r = R
        n_hat = normalize(SA[0.3, -0.4, 0.8])
        @test dot(n_hat, sph(c + (R - 1.0e-6) * n_hat)) ≈
            dot(n_hat, sph(c + (R + 1.0e-6) * n_hat)) rtol = 1.0e-4
        @test current_density(sph, pt_in) == SVector(0.0, 0.0, 0.0)
    end

    let cyl = CylinderMagnet(
            0.3, 0.8, 5.0e5; center = SA[0.1, 0.2, -0.1], normal = SA[0.0, 1.0, 1.0]
        )
        pt = SVector(0.25, 0.3, 0.1)
        @test abs(_fd_div(cyl, pt)) / norm(cyl(pt)) < 1.0e-8
        @test norm(_fd_curl(cyl, pt)) / norm(cyl(pt)) < 1.0e-8
        @test _fd_curl(q -> vector_potential(cyl, q), pt) ≈ cyl(pt) rtol = 1.0e-4
        @test current_density(cyl, pt) == SVector(0.0, 0.0, 0.0)
    end

    let dims = SVector(0.6, 1.0, 1.4),
            M = SVector(2.0e5, -1.0e5, 4.0e5),
            c = SVector(0.1, -0.2, 0.3)
        cube = CuboidMagnet(SVector(0.8, 0.8, 0.8), M)
        @test cube(SVector(0.0, 0.0, 0.0)) ≈ (2 / 3) * Magnetostatics.μ₀ * M rtol = 1.0e-12

        # Demagnetization factor trace sum Nx + Ny + Nz = 1 for any aspect ratio
        box = CuboidMagnet(dims, M; center = c)
        B_c = box(c)
        # B/mu0 = (1 - N) M => sum_k (B_k / (mu0 M_k)) = 3 - (Nx + Ny + Nz) = 2
        trace_check = sum(i -> B_c[i] / (Magnetostatics.μ₀ * M[i]), 1:3)
        @test trace_check ≈ 2.0 rtol = 1.0e-12

        # Far-field dipole limit
        pt_far = c + SVector(20.0, -15.0, 25.0)
        dip = Dipole(prod(dims) * M, c)
        @test box(pt_far) ≈ dip(pt_far) rtol = 1.0e-3

        # Zero divergence and zero curl both inside and outside
        pt_out = c + SVector(0.7, -0.4, 0.9)
        pt_in = c + SVector(0.05, -0.1, 0.2)
        @test abs(_fd_div(box, pt_out)) / norm(box(pt_out)) < 1.0e-8
        @test norm(_fd_curl(box, pt_out)) / norm(box(pt_out)) < 1.0e-8
        @test abs(_fd_div(box, pt_in)) / norm(box(pt_in)) < 1.0e-8
        @test norm(_fd_curl(box, pt_in)) / norm(box(pt_in)) < 1.0e-8
        @test _zero_alloc_call(box, pt_out) == 0

        # Rotated cuboid matches rotate_field
        θ = 0.65
        R = @SMatrix(
            [
                cos(θ) 0.0 sin(θ);
                0.0 1.0 0.0;
                -sin(θ) 0.0 cos(θ)
            ]
        )
        box_rot = CuboidMagnet(dims, M; center = c, rotation = R)
        box_tf = rotate_field(CuboidMagnet(dims, M), R)
        pt_test = c + SVector(0.45, 0.12, -0.6)
        @test box_rot(pt_test) ≈ box_tf(pt_test - c) rtol = 1.0e-12
    end
end

end # module test_solenoid_magnets
