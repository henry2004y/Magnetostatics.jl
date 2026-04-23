module test_potential
using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

@testset "Vector Potential Solver" begin
    # Numerical Curl helper
    function curl_A(solver, source, r, h = 1.0e-5)
        # Central difference curl
        # (dA_z/dy - dA_y/dz, dA_x/dz - dA_z/dx, dA_y/dx - dA_x/dy)

        # Use central difference manually

        x, y, z = r[1], r[2], r[3]

        Ax(p) = solve(solver, source, SVector(p[1], p[2], p[3]))[1]
        Ay(p) = solve(solver, source, SVector(p[1], p[2], p[3]))[2]
        Az(p) = solve(solver, source, SVector(p[1], p[2], p[3]))[3]

        # dAz_dy
        dAz_dy = (Az([x, y + h, z]) - Az([x, y - h, z])) / (2 * h)
        dAy_dz = (Ay([x, y, z + h]) - Ay([x, y, z - h])) / (2 * h)
        Bx = dAz_dy - dAy_dz

        dAx_dz = (Ax([x, y, z + h]) - Ax([x, y, z - h])) / (2 * h)
        dAz_dx = (Az([x + h, y, z]) - Az([x - h, y, z])) / (2 * h)
        By = dAx_dz - dAz_dx

        dAy_dx = (Ay([x + h, y, z]) - Ay([x - h, y, z])) / (2 * h)
        dAx_dy = (Ax([x, y + h, z]) - Ax([x, y - h, z])) / (2 * h)
        Bz = dAy_dx - dAx_dy

        return SVector(Bx, By, Bz)
    end

    solverA = VectorPotential()
    solverB = BiotSavart()

    @testset "Wire" begin
        # Straight wire
        points = [SVector(0.0, 0.0, z) for z in range(-10.0, 10.0, length = 100)]
        wire = Wire(points, 1.0)

        r = SVector(1.0, 0.0, 0.0)

        B_numeric = curl_A(solverA, wire, r)
        B_exact = solve(solverB, wire, r)

        @test isapprox(B_numeric, B_exact, rtol = 1.0e-3)

        A_calc = solve(solverA, wire, r)
        @test solve(solverA, wire, Tuple(r)) ≈ A_calc
        @test solve(solverA, wire, Vector(r)) ≈ A_calc
        @test solverA(wire, r) ≈ A_calc
        @test solverA(wire, Tuple(r)) ≈ A_calc
        @test solverA(wire, Vector(r)) ≈ A_calc
    end

    @testset "Dipole" begin
        M = SVector(0.0, 0.0, 1.0)
        dipole = Dipole(M)

        r = SVector(1.0, 0.0, 0.0) # On x-axis

        # Analytical A for dipole at x-axis (theta=pi/2):
        # A = (mu0/4pi) * (m x r) / r^3
        #   = (mu0/4pi) * (z_hat x x_hat) / r^2
        #   = (mu0/4pi) * y_hat / r^2
        A_exact = SVector(0.0, Magnetostatics.μ0_4π / 1.0^2, 0.0)

        A_calc = solve(solverA, dipole, r)
        @test isapprox(A_calc, A_exact, rtol = 1.0e-5)

        @test solve(solverA, dipole, Tuple(r)) ≈ A_calc
        @test solve(solverA, dipole, Vector(r)) ≈ A_calc
        @test solverA(dipole, r) ≈ A_calc
        @test solverA(dipole, Tuple(r)) ≈ A_calc
        @test solverA(dipole, Vector(r)) ≈ A_calc

        # Numeric curl check
        r_check = SVector(0.5, 0.5, 0.5)
        B_numeric = curl_A(solverA, dipole, r_check)
        B_exact = dipole(r_check)

        @test isapprox(B_numeric, B_exact, rtol = 1.0e-3)
    end

    @testset "CurrentLoop" begin
        R = 1.0
        I = 1.0
        loop = CurrentLoop(R, I, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))

        # Off-axis point
        r = SVector(2.0, 0.0, 1.0)

        B_numeric = curl_A(solverA, loop, r)
        B_exact = getB_loop(r, loop)

        @test isapprox(B_numeric, B_exact, rtol = 1.0e-3)

        A_calc = solve(solverA, loop, r)
        @test solve(solverA, loop, Tuple(r)) ≈ A_calc
        @test solve(solverA, loop, Vector(r)) ≈ A_calc
        @test solverA(loop, r) ≈ A_calc
        @test solverA(loop, Tuple(r)) ≈ A_calc
        @test solverA(loop, Vector(r)) ≈ A_calc
    end

    @testset "HarrisSheet" begin
        B0_tail, L_tail = 20.0, 2.0
        tail = HarrisSheet(B0_tail, L_tail)
        r = SVector(1.0, 2.0, 3.0)

        A_calc = solve(solverA, tail, r)
        Ay_expected = -B0_tail * L_tail * log(cosh(r[3] / L_tail))
        @test A_calc ≈ SVector(0.0, Ay_expected, 0.0)

        @test solve(solverA, tail, Tuple(r)) ≈ A_calc
        @test solve(solverA, tail, Vector(r)) ≈ A_calc
        @test solverA(tail, r) ≈ A_calc
        @test solverA(tail, Tuple(r)) ≈ A_calc
        @test solverA(tail, Vector(r)) ≈ A_calc

        # Numeric curl check
        B_numeric = curl_A(solverA, tail, r)
        B_exact = tail(r)
        @test isapprox(B_numeric, B_exact, rtol = 1.0e-3)

        # Numerical stability check for large z
        r_large = SVector(0.0, 0.0, 1000.0 * L_tail) # z/L = 1000, would overflow cosh
        A_large = solve(solverA, tail, r_large)
        @test all(isfinite, A_large)
        @test A_large[2] ≈ -B0_tail * L_tail * (1000.0 - log(2.0)) rtol = 1.0e-5
    end

    @testset "UniformField" begin
        B0 = SVector(0.0, 0.0, -10.0)
        field = UniformField(B0)
        r = SVector(1.0, 2.0, 3.0)

        A_calc = solve(solverA, field, r)
        A_expected = 0.5 * cross(B0, r)
        @test A_calc ≈ A_expected

        @test solve(solverA, field, Tuple(r)) ≈ A_calc
        @test solve(solverA, field, Vector(r)) ≈ A_calc
        @test solverA(field, r) ≈ A_calc
        @test solverA(field, Tuple(r)) ≈ A_calc
        @test solverA(field, Vector(r)) ≈ A_calc

        # Numeric curl check
        B_numeric = curl_A(solverA, field, r)
        B_exact = field(r)
        @test isapprox(B_numeric, B_exact, rtol = 1.0e-3)
    end

    @testset "vector_potential convenience" begin
        B0 = SVector(1.0, 2.0, 3.0)
        field = UniformField(B0)
        r = SVector(0.1, 0.2, 0.3)
        @test vector_potential(field, r) ≈ solve(solverA, field, r)
        @test vector_potential(field, Tuple(r)) ≈ solve(solverA, field, r)
        @test vector_potential(field, Vector(r)) ≈ solve(solverA, field, r)
    end
end

end # module test_potential
