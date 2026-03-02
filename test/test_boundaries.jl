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

@testset "PeriodicBoundary" begin
    wire = make_wire()
    period = SVector(4.0, 4.0, 8.0)
    bc = PeriodicBoundary(period)
    solver = BiotSavart(bc)
    r = SVector(1.0, 0.0, 0.0)

    B_periodic = solve(solver, wire, r; nreplica = 1)

    # Manual reference: sum over the (2*1+1)^3 = 27 image sources
    images = Magnetostatics.periodic_sources(wire, bc; nreplica = 1)
    B_ref = sum(img -> solve(BiotSavart(), img, r), images)

    @test B_periodic ≈ B_ref atol = 1.0e-14
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

end # module test_boundaries
