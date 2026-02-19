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

    @testset "Dipole Fieldline" begin
        # Test default parameters
        ϕ = 0.0
        x, y, z = Magnetostatics.dipole_fieldline(ϕ)

        # Verify dimensions
        @test all(length.((x, y, z)) .== 100)

        # Check start and end points (should be at origin)
        @test all(isapprox.(getindex.((x, y, z), 1), 0.0; atol = 1.0e-10)) &&
            all(isapprox.(getindex.((x, y, z), 100), 0.0; atol = 1.0e-10))

        # Check maximum radius (L-shell)
        r_vals = sqrt.(x .^ 2 + y .^ 2 + z .^ 2)
        @test maximum(r_vals) ≈ 2.5 atol = 1.0e-2

        # Test with custom arguments
        L, nP = 4.0, 51
        x2, y2, z2 = Magnetostatics.dipole_fieldline(π / 2, L, nP)
        r_vals2 = sqrt.(x2 .^ 2 + y2 .^ 2 + z2 .^ 2)
        @test maximum(r_vals2) ≈ 4.0 atol = 1.0e-2

        # At phi=pi/2, point with max r should be at (0, L, 0)
        max_idx = argmax(r_vals2)
        @test isapprox(x2[max_idx], 0.0; atol = 1.0e-10) &&
            isapprox(y2[max_idx], 4.0; atol = 1.0e-2) &&
            isapprox(z2[max_idx], 0.0; atol = 1.0e-2)
    end
end
