module test_poisson

using Test
using Magnetostatics
using LinearSolve: KLUFactorization, KrylovJL_CG
using StaticArrays

# Small uniform-Jz box: A_z should be non-negative (source is positive,
# Dirichlet = 0 on boundary), symmetric in x-y, and consistent between
# the direct factorisation and the Krylov CG solver.

const Nx, Ny, Nz = 16, 16, 16
const L = 1.0

const xs = range(0.0, L; length = Nx + 1)[1:(end - 1)]
const ys = range(0.0, L; length = Ny + 1)[1:(end - 1)]
const zs = range(0.0, L; length = Nz + 1)[1:(end - 1)]

const J_grid = let
    J = zeros(Float64, 3, Nx, Ny, Nz)
    J[3, :, :, :] .= 1.0  # uniform J_z
    J
end

@testset "PoissonSolver – direct (KLUFactorization)" begin
    solver = PoissonSolver(xs, ys, zs; alg = KLUFactorization())
    A = solve(solver, J_grid)

    @test size(A) == (3, Nx, Ny, Nz)

    # No Jx / Jy source → Ax, Ay should be zero everywhere
    @test all(iszero, A[1, :, :, :])
    @test all(iszero, A[2, :, :, :])

    Az = A[3, :, :, :]

    # All interior values non-negative (positive source, zero BC)
    @test all(>=(0.0), Az)

    # x-y symmetry: A_z(i,j,k) ≈ A_z(j,i,k) because the box is square
    @test isapprox(Az, permutedims(Az, (2, 1, 3)); atol = 1.0e-10)
end

@testset "PoissonSolver – Krylov CG matches direct" begin
    solver_direct = PoissonSolver(xs, ys, zs; alg = KLUFactorization())
    solver_cg = PoissonSolver(xs, ys, zs; alg = KrylovJL_CG())

    A_direct = solve(solver_direct, J_grid)
    A_cg = solve(solver_cg, J_grid)

    @test isapprox(A_cg, A_direct; rtol = 1.0e-4)
end

end # module test_poisson
