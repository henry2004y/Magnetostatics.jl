module test_boundaries

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

# Shared geometry: a straight wire along the z-axis
function make_wire(; I = 1.0)
    pts = [SVector(0.0, 0.0, -2.0), SVector(0.0, 0.0, 2.0)]
    return Wire(pts, I)
end

@testset "OpenBoundary" begin
    wire = make_wire()
    solver = BiotSavart(OpenBoundary())
    r = SVector(1.0, 0.0, 0.0)

    B_open = solve(solver, wire, r)
    B_bare = solve(BiotSavart(), wire, r)   # default is also OpenBoundary
    @test B_open ≈ B_bare atol = 1.0e-14
end

@testset "ConductingWall" begin
    # Wire at x=0 parallel to z-axis; conducting wall at x=2
    # The image wire is at x=4 with reversed current.
    wire = make_wire()
    bc = ConductingWall(1, 2.0)
    solver = BiotSavart(bc)

    let r = SVector(2.0, 0.5, 0.0)
        # At the wall surface (x=2) the contribution from the real wire and
        # its image to Bx (normal to the wall) must cancel.
        B_total = solve(solver, wire, r)
        # Bx is the component normal to the wall
        @test abs(B_total[1]) < 1.0e-10
    end

    let
        img = Magnetostatics.image_source(wire, bc)
        # Image wire should be at x = 2*2 - 0 = 4 and have reversed current
        @test img.points[1][1] ≈ 4.0
        @test img.current ≈ -wire.current
    end
end

@testset "PrescribedBoundary" begin
    wire = make_wire()
    B_bg_val = SVector(0.0, 0.0, 1.0e-3)        # uniform background field
    bg_func = r -> B_bg_val
    bc = PrescribedBoundary(bg_func)
    solver = BiotSavart(bc)
    r = SVector(1.0, 0.0, 0.0)

    B_total = solve(solver, wire, r)
    B_source = solve(BiotSavart(), wire, r)

    @test B_total ≈ B_source + B_bg_val atol = 1.0e-14
end

@testset "ConductingSphere" begin
    # Wire along z-axis; conducting sphere centred at origin, radius 0.5.
    # The wire is outside the sphere (closest approach at x=0, y=0, R=0.5<1).
    wire = make_wire()
    bc = ConductingSphere(SVector(0.0, 0.0, 0.0), 0.5, 16, 32)

    mesh = compute_surface_current(bc, wire)
    @test length(mesh.centers) == 16 * 32

    # Verify BC: n̂ · B_total ≈ 0 at every patch centre
    bare = BiotSavart()
    max_residual = maximum(eachindex(mesh.centers)) do i
        r_p = mesh.centers[i]
        B_inc = solve(bare, wire, r_p)
        B_K = solve(bare, mesh, r_p)
        abs(dot(mesh.normals[i], B_inc + B_K))
    end

    # Residual should be small relative to the peak incident field magnitude
    B_peak = maximum(
        i -> norm(solve(bare, wire, mesh.centers[i])),
        eachindex(mesh.centers)
    )
    @test max_residual / B_peak < 0.05
end

@testset "MultiSource, Loop Wall Image, PeriodicBoundary and Filament Softening" begin
    let l1 = CurrentLoop(0.8, 12.0, SA[0.0, 0.0, 0.4], SA[0.0, 0.0, 1.0]),
            l2 = CurrentLoop(0.8, 12.0, SA[0.0, 0.0, -0.4], SA[0.0, 0.0, 1.0]),
            ms = MultiSource(l1, l2),
            w1 = discretize_loop(l1, 160),
            w2 = discretize_loop(l2, 160),
            r = SVector(0.25, -0.15, 0.1)
        @test solve(BiotSavart(), ms, r) ≈ getB_loop(r, l1) + getB_loop(r, l2)
        @test solve(BiotSavart(), [l1, l2], r) ≈ getB_loop(r, l1) + getB_loop(r, l2)
        @test solve(VectorPotential(), ms, r) ≈
            solve(VectorPotential(), l1, r) + solve(VectorPotential(), l2, r)

        # ConductingWall image for CurrentLoop matches discretized Wire image
        wall = ConductingWall(3, 0.0)
        solver_w = BiotSavart(wall)
        @test abs(solve(solver_w, l1, SVector(0.3, -0.2, 0.0))[3]) < 1.0e-12
        @test solve(solver_w, l1, SVector(0.3, -0.2, 0.7)) ≈
            solve(solver_w, w1, SVector(0.3, -0.2, 0.7)) rtol = 1.0e-3

        # Tilted loop perpendicular to wall
        l_tilt = CurrentLoop(0.5, 5.0, SA[0.0, 0.0, 0.8], SA[1.0, 0.2, 0.3])
        @test abs(solve(solver_w, l_tilt, SVector(0.1, 0.4, 0.0))[3]) < 1.0e-12

        # PeriodicBoundary sums explicit lattice images
        pbc = PeriodicBoundary((2.0, 2.0, 0.0); n_images = (1, 1, 0))
        B_pbc = solve(BiotSavart(pbc), ms, r)
        B_manual = sum(
            solve(BiotSavart(), ms, r - SVector(2.0 * kx, 2.0 * ky, 0.0))
                for kx in -1:1, ky in -1:1
        )
        @test B_pbc ≈ B_manual

        # Finite-core filament softening eliminates on-wire singularity
        straight = Wire([SA[0.0, 0.0, -5.0], SA[0.0, 0.0, 5.0]], 10.0)
        bs_soft = BiotSavart(; core_radius = 0.05)
        @test solve(bs_soft, straight, SA[0.0, 0.0, 0.0]) == SVector(0.0, 0.0, 0.0)
        B_core = solve(bs_soft, straight, SA[0.05, 0.0, 0.0])
        @test isapprox(
            B_core[2], Magnetostatics.μ₀ * 10.0 / (4 * π * 0.05); rtol = 1.0e-3
        )
        @test all(isfinite, solve(bs_soft, l1, SA[0.8, 0.0, 0.4]))

        # ConductingSphere BEM accepts MultiSource directly
        bc_s = ConductingSphere(SA[0.0, 0.0, 0.0], 0.2, 8, 16)
        mesh_ms = compute_surface_current(bc_s, [w1, w2])
        @test length(mesh_ms.centers) == 8 * 16
    end
end

end # module test_boundaries
