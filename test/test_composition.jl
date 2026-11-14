module test_composition

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

function _zero_alloc_eval(field, r)
    field(r)
    return @allocated field(r)
end

function _zero_alloc_j(field, r)
    current_density(field, r)
    return @allocated current_density(field, r)
end

@testset "Field Algebra and Transforms" begin
    let d1 = Dipole(SA[0.0, 0.0, 1.0], SA[0.5, 0.0, 0.0]),
            u1 = UniformField(SA[0.0, 0.0, 2.0e-6]),
            h1 = HarrisSheet(1.5e-5, 0.4),
            r = SVector(0.8, -0.3, 0.25)
        f = d1 + 2.0 * u1 - 0.5 * h1
        @test f isa CompositeField
        @test f(r) ≈ d1(r) + 2.0 * u1(r) - 0.5 * h1(r)
        @test current_density(f, r) ≈
            current_density(d1, r) + 2.0 * current_density(u1, r) -
            0.5 * current_density(h1, r)
        @test vector_potential(f, r) ≈
            vector_potential(d1, r) + 2.0 * vector_potential(u1, r) -
            0.5 * vector_potential(h1, r)
        @test _zero_alloc_eval(f, r) == 0
        @test _zero_alloc_j(f, r) == 0
        @test NullField() + d1 === d1
        @test d1 + NullField() === d1
        @test 3.0 * NullField() === NullField()
        @test (2.0 * (3.0 * u1))(r) ≈ 6.0 * u1(r)
        @test (u1 / 2.0)(r) ≈ 0.5 * u1(r)
    end

    let h = HarrisSheet(2.0, 0.5),
            offset = SVector(0.2, -0.4, 0.3),
            θ = 0.7,
            R = @SMatrix(
            [
                cos(θ) -sin(θ) 0.0;
                sin(θ) cos(θ) 0.0;
                0.0 0.0 1.0
            ]
        )
        r = SVector(0.6, 0.9, -0.1)
        f_tr = rotate_field(translate_field(h, offset), R)
        r_check = SVector(1.1, -0.3, 0.45)
        @test f_tr(r_check) ≈ R * h(R' * r_check - offset)
        @test current_density(f_tr, r_check) ≈
            R * current_density(h, R' * r_check - offset)
        @test _fd_curl(p -> vector_potential(f_tr, p), r_check) ≈
            f_tr(r_check) rtol = 1.0e-5
        @test _fd_curl(f_tr, r_check) ≈
            Magnetostatics.μ₀ * current_density(f_tr, r_check) rtol = 1.0e-5
    end
end

@testset "Analytical Vector Potential and Current Density Completeness" begin
    let r_gen = SVector(1.2, -0.5, 0.4),
            r_cyl = SVector(3.4, 0.6, 0.35),
            r_wire_in = SVector(0.15, -0.1, 0.2),
            r_wire_out = SVector(0.8, 0.3, -0.2),
            r_parker = SVector(2.0e11, 1.0e11, 0.8e11)
        loop = CurrentLoop(1.0, 2.5, SA[0.0, 0.0, 0.1], SA[0.1, 0.2, 1.0])
        models = [
            (Dipole(SA[1.0, -2.0, 3.0], SA[0.1, 0.2, -0.1]), r_gen),
            (AsymmetricHarrisSheet(2.0, -1.0, 0.8), r_gen),
            (ForceFreeHarrisSheet(1.5, 0.7), r_gen),
            (BifurcatedHarrisSheet(1.8, 0.6, 0.9), r_gen),
            (FadeevIsland(1.2, 0.8, 1.1, 0.35), r_gen),
            (GuideFieldHarrisSheet(1.4, 0.6; Bg = 0.5, Bn = 0.2), r_gen),
            (CurrentLoopAnalytic(loop), r_gen),
            (MagneticMirror(1.2, 0.8, 1.5), r_gen),
            (MagneticBottle(1.6, 0.7, 1.0, 2.0, 0.8), r_gen),
            (ZPinch(20.0, 0.4; center = SA[0.1, -0.1]), r_wire_in),
            (ZPinch(20.0, 0.4; center = SA[0.1, -0.1]), r_wire_out),
            (ScrewPinch(25.0, 0.4, 0.8; center = SA[0.05, -0.05]), r_wire_in),
            (ScrewPinch(25.0, 0.4, 0.8; center = SA[0.05, -0.05]), r_wire_out),
            (InfiniteWire(15.0, 0.4, SA[0.0, 0.0, 0.0], SA[0.2, 0.1, 1.0]), r_wire_in),
            (InfiniteWire(15.0, 0.4, SA[0.0, 0.0, 0.0], SA[0.2, 0.1, 1.0]), r_wire_out),
            (TaylorState(1.2, 1.7; center = SA[0.05, 0.1]), r_gen),
            (SolovevEquilibrium(2.0, 3.0, 1.5, 1.2), r_cyl),
            (GoldHoyleFluxTube(1.3, 0.9; center = SA[0.1, -0.2]), r_gen),
            (FLARECoils(; I_ef = 20.0e3, I_p = 5.0e3), r_gen),
            (ParkerSpiral(5.0e-9, 1.496e11), r_parker),
        ]
        for (field, pt) in models
            h_step = field isa ParkerSpiral ? 1.0e5 : 1.0e-5
            @test _fd_curl(p -> vector_potential(field, p), pt; h = h_step) ≈
                field(pt) rtol = 1.0e-4
            if !(field isa FadeevIsland)
                @test _fd_curl(field, pt; h = h_step) ≈
                    Magnetostatics.μ₀ * current_density(field, pt) rtol = 1.0e-3 atol =
                    1.0e-7 * norm(field(pt))
            end
        end
    end
end

@testset "GridMagneticField" begin
    let xs = range(-1.0, 1.0; length = 21),
            ys = range(-1.0, 1.0; length = 21),
            zs = range(-1.0, 1.0; length = 21)
        B_lin = zeros(3, length(xs), length(ys), length(zs))
        for k in eachindex(zs), j in eachindex(ys), i in eachindex(xs)
            x, y, z = xs[i], ys[j], zs[k]
            B_lin[1, i, j, k] = 2x - 3y + 0.5z
            B_lin[2, i, j, k] = -x + 1.5y + 4z
            B_lin[3, i, j, k] = 0.25x - y + 2z
        end
        gf_lin = GridMagneticField(B_lin, xs, ys, zs)
        pt = SVector(0.37, -0.42, 0.61)
        @test gf_lin(pt) ≈ SVector(
            2pt[1] - 3pt[2] + 0.5pt[3],
            -pt[1] + 1.5pt[2] + 4pt[3],
            0.25pt[1] - pt[2] + 2pt[3]
        ) atol = 1.0e-14
        J_exp = SVector(-5.0, 0.25, 2.0) / Magnetostatics.μ₀
        @test current_density(gf_lin, pt) ≈ J_exp rtol = 1.0e-12
        @test _zero_alloc_eval(gf_lin, pt) == 0
        @test _zero_alloc_j(gf_lin, pt) == 0

        @test_throws DomainError gf_lin(SVector(2.0, 0.0, 0.0))
        gf_zero = GridMagneticField(B_lin, xs, ys, zs; extrapolate = :zero)
        @test gf_zero(SVector(2.0, 0.0, 0.0)) == SVector(0.0, 0.0, 0.0)
        gf_clamp = GridMagneticField(B_lin, xs, ys, zs; extrapolate = :clamp)
        @test gf_clamp(SVector(5.0, 0.0, 0.0)) ≈ gf_lin(SVector(1.0, 0.0, 0.0))

        # Smooth nonlinear field comparing cubic vs linear accuracy
        f_smooth = r -> SVector(sin(r[2]), cos(r[3]), exp(0.3 * r[1]))
        B_nl = zeros(3, length(xs), length(ys), length(zs))
        for k in eachindex(zs), j in eachindex(ys), i in eachindex(xs)
            B_nl[:, i, j, k] = f_smooth(SVector(xs[i], ys[j], zs[k]))
        end
        gf_c = GridMagneticField(B_nl, xs, ys, zs; interpolation = :cubic)
        gf_l = GridMagneticField(B_nl, xs, ys, zs; interpolation = :linear)
        p_test = SVector(0.13, -0.27, 0.31)
        err_c = norm(gf_c(p_test) - f_smooth(p_test))
        err_l = norm(gf_l(p_test) - f_smooth(p_test))
        @test err_c < 0.05 * err_l
        @test _zero_alloc_eval(gf_c, p_test) == 0
    end
end

end # module test_composition
