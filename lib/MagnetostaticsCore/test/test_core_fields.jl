module test_core_fields

using LinearAlgebra: cross, dot, norm
using MagnetostaticsCore
using StaticArrays: @SMatrix, SVector
using Test

function finite_diff_curl(f, r; h::Float64 = 1.0e-6)
    ex = SVector(h, 0.0, 0.0)
    ey = SVector(0.0, h, 0.0)
    ez = SVector(0.0, 0.0, h)
    dA_dx = (f(r + ex) - f(r - ex)) / (2 * h)
    dA_dy = (f(r + ey) - f(r - ey)) / (2 * h)
    dA_dz = (f(r + ez) - f(r - ez)) / (2 * h)
    return SVector(
        dA_dy[3] - dA_dz[2],
        dA_dz[1] - dA_dx[3],
        dA_dx[2] - dA_dy[1]
    )
end

@testset "Dependencies remain lightweight" begin
    proj = joinpath(@__DIR__, "..", "Project.toml")
    txt = read(proj, String)
    @test !occursin("FFTW", txt)
    @test !occursin("LinearSolve", txt)
    @test !occursin("SparseArrays", txt)
    @test !occursin("SpecialFunctions", txt)
end

@testset "NullField, UniformField, and composition" begin
    nf = NullField()
    r = SVector(0.3, -0.4, 0.5)
    @test nf(r) == SVector(0.0, 0.0, 0.0)
    @test vector_potential(nf, r) == SVector(0.0, 0.0, 0.0)
    B0, A0 = field_and_potential(nf, r)
    @test B0 == SVector(0.0, 0.0, 0.0)
    @test A0 == SVector(0.0, 0.0, 0.0)

    uf = UniformField([1.0, -2.0, 3.0])
    @test uf(r) == SVector(1.0, -2.0, 3.0)
    @test uf(r[1], r[2], r[3]) == SVector(1.0, -2.0, 3.0)
    @test finite_diff_curl(x -> vector_potential(uf, x), r) ≈ uf(r) atol = 1.0e-9

    uf2 = UniformField([0.0, 1.0, -1.0])
    cf = uf + 2.0 * uf2 - uf2
    @test cf(r) ≈ uf(r) + uf2(r)
    @test vector_potential(cf, r) ≈ vector_potential(uf, r) + vector_potential(uf2, r)

    rot_mat = @SMatrix [
        0.0 -1.0 0.0;
        1.0 0.0 0.0;
        0.0 0.0 1.0
    ]
    tf = rotate_field(translate_field(uf, SVector(1.0, 0.0, 0.0)), rot_mat)
    @test tf(r) ≈ SVector(2.0, 1.0, 3.0)
    @test finite_diff_curl(x -> vector_potential(tf, x), r) ≈ tf(r) atol = 1.0e-9

    uf32 = UniformField(SVector(1.0f0, -2.0f0, 3.0f0))
    r32 = SVector(0.3f0, -0.4f0, 0.5f0)
    rot32 = Float32.(rot_mat)
    @test eltype(TransformedField(uf32, SVector(1.0f0, 0.0f0, 0.0f0))(r32)) === Float32
    @test eltype(rotate_field(uf32, rot32)(r32)) === Float32
end

@testset "Dipole and LineDipole2D" begin
    d = Dipole([0.0, 0.0, 1.0]; center = [0.1, -0.2, 0.0], r_min = 0.25, k = 1.0)
    r_far = SVector(1.5, 0.4, -0.3)
    B, A = field_and_potential(d, r_far)
    @test B ≈ d(r_far)
    @test A ≈ vector_potential(d, r_far)
    @test finite_diff_curl(x -> vector_potential(d, x), r_far) ≈ B atol = 1.0e-9
    @test current_density(d, r_far) == SVector(0.0, 0.0, 0.0)

    # Inside softening radius r_min = 0.25: bounded field at center
    @test norm(d(SVector(0.1, -0.2, 0.0))) < 100.0
    r_in = SVector(0.12, -0.18, 0.05)
    @test norm(d(r_in)) == norm(field_and_potential(d, r_in)[1])

    ld = LineDipole2D([0.0, 1.5], [0.2, -0.1]; r_min = 0.2, k = 1.0)
    r2d = SVector(1.2, 0.7, 0.3)
    B2, A2 = field_and_potential(ld, r2d)
    @test B2 ≈ ld(r2d)
    @test A2 ≈ vector_potential(ld, r2d)
    @test B2[3] == 0.0
    @test finite_diff_curl(x -> vector_potential(ld, x), r2d) ≈ B2 atol = 1.0e-9
    @test current_density(ld, r2d) == SVector(0.0, 0.0, 0.0)
    dr = zeros(3)
    @test ld(dr, r2d) === dr
    @test dr ≈ Vector(B2)
end

@testset "InfiniteWire, WireSegment, CurrentLoopAnalytic, ToroidalField" begin
    r = SVector(0.8, -0.5, 0.4)

    # InfiniteWire with Plummer core and k = 1.0
    w_core = InfiniteWire(
        2.5, 0.0, [0.1, 0.2, 0.0], [0.0, 0.0, 1.0];
        core_radius = 0.2, k = 1.0
    )
    Bw, Aw = field_and_potential(w_core, r)
    @test Bw ≈ w_core(r)
    @test Aw ≈ vector_potential(w_core, r)
    @test finite_diff_curl(x -> vector_potential(w_core, x), r) ≈ Bw atol = 1.0e-9
    @test current_density(w_core, SVector(0.1, 0.2, 0.0)) == SVector(0.0, 0.0, 0.0)

    # WireSegment with Plummer core
    seg = WireSegment(
        [0.0, 0.0, -1.0], [0.0, 0.0, 1.0], 3.0;
        core_radius = 0.15, k = 1.0
    )
    Bs, As = field_and_potential(seg, r)
    @test Bs ≈ seg(r)
    @test As ≈ vector_potential(seg, r)
    @test finite_diff_curl(x -> vector_potential(seg, x), r) ≈ Bs atol = 1.0e-9
    @test current_density(seg, r) == SVector(0.0, 0.0, 0.0)

    # CurrentLoopAnalytic with Plummer core and k = 1.0
    lp = CurrentLoop(1.0, 2.0, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
    la_soft = CurrentLoopAnalytic(lp; core_radius = 0.1, k = 1.0)
    Bl, Al = field_and_potential(la_soft, r)
    @test Bl ≈ la_soft(r)
    @test Al ≈ vector_potential(la_soft, r)
    @test finite_diff_curl(x -> vector_potential(la_soft, x), r) ≈ Bl atol = 1.0e-8
    @test finite_diff_curl(
        x -> vector_potential(CurrentLoopAnalytic(lp), x), r
    ) ≈ getB_loop(r, lp) atol = 1.0e-12

    # Near-axis CurrentLoopAnalytic stays finite
    @test all(isfinite, la_soft(SVector(0.0, 0.0, 0.5)))

    # ToroidalField with Plummer core
    tf = ToroidalField(1.8; center = [0.1, -0.1, 0.0], core_radius = 0.25)
    Bt, At = field_and_potential(tf, r)
    @test Bt ≈ tf(r)
    @test At ≈ vector_potential(tf, r)
    @test finite_diff_curl(x -> vector_potential(tf, x), r) ≈ Bt atol = 1.0e-9
end

@testset "SphericalHarmonicField and Workspace" begin
    # Axial dipole g[2, 1] = g_1^0 should match Dipole along z
    g = zeros(4, 4)
    h = zeros(4, 4)
    g[2, 1] = 1.0
    r_ref = 1.0
    sh = SphericalHarmonicField(g, h, r_ref)
    ws = SphericalHarmonicWorkspace(sh.nmax)

    # Schmidt g_1^0 corresponds to dipole M = [0, 0, +g_1^0 * r_ref^3 / k] with B = (3(M⋅r̂)r̂ - M)/r^3
    d_eq = Dipole([0.0, 0.0, 1.0]; k = 1.0)
    r_out = SVector(1.4, -0.7, 0.9)
    B_sh = sh(r_out)
    @test B_sh ≈ d_eq(r_out) rtol = 1.0e-12
    @test evaluate!(ws, sh, r_out) ≈ B_sh

    # Mixed degree n=1,2,3 (tilted dipole + quadrupole)
    g[2, 2] = 0.3
    h[2, 2] = -0.2
    g[3, 1] = 0.4
    g[3, 2] = -0.15
    h[3, 3] = 0.25
    sh2 = SphericalHarmonicField(g, h, r_ref; center = [0.2, -0.1, 0.05])
    B_ws = evaluate!(ws, sh2, r_out)
    @test B_ws ≈ sh2(r_out) rtol = 1.0e-13

    # Clamping inside r_ref keeps magnitude finite
    @test all(isfinite, evaluate!(ws, sh2, SVector(0.4, 0.1, 0.1)))

    # Zero-allocation check on preallocated workspace
    allocs = @allocated evaluate!(ws, sh2, r_out)
    @test allocs == 0
end

end
