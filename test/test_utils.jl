using Magnetostatics
using StaticArrays
using LinearAlgebra
using Test

@testset "Utils" begin
    # Define a function to verify the normal of the loop
    function verify_loop_normal(wire::Wire, expected_normal)
        points = wire.points
        # Approximate normal from the first 3 points (assuming they form a plane)
        # n ~ (p2 - p1) x (p3 - p1)
        p1 = points[1]
        p2 = points[2]
        p3 = points[3]

        v1 = p2 - p1
        v2 = p3 - p1

        computed_normal = normalize(cross(v1, v2))

        # Check if parallel
        return isapprox(computed_normal, normalize(expected_normal); atol = 1.0e-2)
    end

    @testset "Discretize Loop Rotation" begin
        # Test 1: Standard z-axis (should work as before)
        loop_z = discretize_loop(1.0, 100, 1.0; normal = SVector(0.0, 0.0, 1.0))
        @test verify_loop_normal(loop_z, SVector(0.0, 0.0, 1.0))

        # Test 2: x-axis
        loop_x = discretize_loop(1.0, 100, 1.0; normal = SVector(1.0, 0.0, 0.0))
        @test verify_loop_normal(loop_x, SVector(1.0, 0.0, 0.0))

        # Test 3: Arbitrary axis
        n_arb = SVector(1.0, 1.0, 1.0)
        loop_arb = discretize_loop(1.0, 100, 1.0; normal = n_arb)
        @test verify_loop_normal(loop_arb, n_arb)

        # Test 4: -z axis
        loop_neg_z = discretize_loop(1.0, 100, 1.0; normal = SVector(0.0, 0.0, -1.0))
        @test verify_loop_normal(loop_neg_z, SVector(0.0, 0.0, -1.0))
    end
end
