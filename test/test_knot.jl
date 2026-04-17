module test_knot

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

@testset "Torus Knot and Trefoil Knot" begin
    # Test creation
    R, r = 1.0, 0.3
    current = 1.0
    knot = TorusKnot(R, r, 2, 3, current)
    @test knot.R == R && knot.r == r && knot.p == 2 && knot.q == 3 && knot.current == current

    # Test trefoil convenience constructor
    trefoil = TrefoilKnot(R, r, current)
    @test trefoil.p == 2 && trefoil.q == 3

    # Test discretization
    n_segments = 100
    wire = discretize_knot(knot, n_segments)
    @test length(wire.points) == n_segments + 1 && wire.current == current

    # Test closed loop
    @test wire.points[1] ≈ wire.points[end]

    # Test magnetic field calculation (basic sanity check)
    solver = BiotSavart()
    # At the center (0,0,0), for a symmetric torus knot, B should be along the normal (z-axis)
    B = solver(wire, [0.0, 0.0, 0.0])
    @test B[1] ≈ 0.0 atol = 1.0e-10
    @test B[2] ≈ 0.0 atol = 1.0e-10
    @test B[3] != 0.0

    # Test arbitrary orientation
    normal = [1.0, 1.0, 1.0]
    up = [1.0, -1.0, 0.0]
    knot_rot = TorusKnot(R, r, 2, 3, current, [0, 0, 0], normal, up)
    @test dot(knot_rot.normal, knot_rot.up) ≈ 0.0 atol = 1.0e-10

    wire_rot = Wire(knot_rot, n_segments)
    B_rot = solver(wire_rot, [0.0, 0.0, 0.0])
    # The B field at the center should be parallel to the normal
    @test normalize(B_rot) ≈ normalize(SVector{3}(normal)) || normalize(B_rot) ≈ -normalize(SVector{3}(normal))
end

end # module test_knot
