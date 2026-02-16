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
        # z-axis
        loop_z = discretize_loop(1.0, 100, 1.0; normal = SVector(0.0, 0.0, 1.0))
        @test verify_loop_normal(loop_z, SVector(0.0, 0.0, 1.0))

        # x-axis
        loop_x = discretize_loop(1.0, 100, 1.0; normal = SVector(1.0, 0.0, 0.0))
        @test verify_loop_normal(loop_x, SVector(1.0, 0.0, 0.0))

        # Arbitrary axis
        n_arb = SVector(1.0, 1.0, 1.0)
        loop_arb = discretize_loop(1.0, 100, 1.0; normal = n_arb)
        @test verify_loop_normal(loop_arb, n_arb)

        # -z axis
        loop_neg_z = discretize_loop(1.0, 100, 1.0; normal = SVector(0.0, 0.0, -1.0))
        @test verify_loop_normal(loop_neg_z, SVector(0.0, 0.0, -1.0))

        # CurrentLoop object dispatch
        cl = CurrentLoop(1.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0))
        loop_obj = discretize_loop(cl, 100)
        @test verify_loop_normal(loop_obj, SVector(0.0, 1.0, 0.0))
    end

    @testset "Wire Current Setting" begin
        Nx, Ny, Nz = 10, 10, 10
        x = range(-1.0, 1.0, length = Nx)
        y = range(-1.0, 1.0, length = Ny)
        z = range(-1.0, 1.0, length = Nz)
        current = 1.0
        width = 0.5

        # Test in-place
        J_inplace = zeros(Float64, 3, Nx, Ny, Nz)
        set_current_wire!(
            J_inplace, x, y, z,
            SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0),
            current, width
        )

        # Test out-of-place
        J_outofplace = set_current_wire(
            x, y, z,
            SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0),
            current, width
        )

        @test J_inplace ≈ J_outofplace

        # Verify direction (Should be purely z-directed for z-wire)
        @test all(isapprox.(J_inplace[1, :, :, :], 0.0, atol = 1.0e-10)) &&
            all(isapprox.(J_inplace[2, :, :, :], 0.0, atol = 1.0e-10)) &&
            any(abs.(J_inplace[3, :, :, :]) .> 0.0)
    end
end
