module test_poisson

using Test
using Magnetostatics
using LinearAlgebra: issymmetric, isposdef, norm
using LinearSolve: KLUFactorization, KrylovJL_CG, CholeskyFactorization
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

@testset "PoissonSolver – SPD, Cholesky, DSTFactorization & inhomogeneous BCs" begin
    xs_aniso = range(-0.5, 0.5; length = 10)
    ys_aniso = range(-0.6, 0.6; length = 12)
    zs_aniso = range(-0.4, 0.4; length = 8)
    nx, ny, nz = length(xs_aniso), length(ys_aniso), length(zs_aniso)

    solver_klu = PoissonSolver(xs_aniso, ys_aniso, zs_aniso; alg = KLUFactorization())
    solver_chol = PoissonSolver(xs_aniso, ys_aniso, zs_aniso; alg = CholeskyFactorization())
    solver_cg = PoissonSolver(xs_aniso, ys_aniso, zs_aniso; alg = KrylovJL_CG())
    solver_dst = PoissonSolver(xs_aniso, ys_aniso, zs_aniso; alg = DSTFactorization())

    @test solver_dst.L === nothing
    @test issymmetric(solver_klu.L)

    let xs_s = range(0.0, 1.0; length = 6)
        L_small = PoissonSolver(xs_s, xs_s, xs_s; alg = KLUFactorization()).L
        @test issymmetric(L_small)
        @test isposdef(Matrix(L_small))
    end

    J3 = zeros(Float64, 3, nx, ny, nz)
    for k in 2:(nz - 1), j in 2:(ny - 1), i in 2:(nx - 1)
        x, y, z = xs_aniso[i], ys_aniso[j], zs_aniso[k]
        J3[1, i, j, k] = exp(-4 * (x^2 + y^2 + z^2))
        J3[2, i, j, k] = (x - 0.1) * exp(-3 * (x^2 + y^2 + z^2))
        J3[3, i, j, k] = -0.7 * cos(y) * exp(-2 * z^2)
    end

    A_klu = Magnetostatics.solve(solver_klu, J3)
    A_chol = Magnetostatics.solve(solver_chol, J3)
    A_cg = Magnetostatics.solve(solver_cg, J3)
    A_dst = Magnetostatics.solve(solver_dst, J3)

    @test isapprox(A_chol, A_klu; atol = 1.0e-13, rtol = 1.0e-12)
    @test isapprox(A_dst, A_klu; atol = 1.0e-13, rtol = 1.0e-12)
    @test isapprox(A_cg, A_klu; atol = 1.0e-12, rtol = 1.0e-3)
    @test norm(A_cg - A_klu) <= 5.0e-4 * norm(A_klu)

    # Exact recovery of discrete harmonic function under inhomogeneous Dirichlet BC
    A_harm = zeros(Float64, 3, nx, ny, nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        x, y, z = xs_aniso[i], ys_aniso[j], zs_aniso[k]
        A_harm[1, i, j, k] = 2x - 3y + 0.5z
        A_harm[2, i, j, k] = x * y - 1.5z
        A_harm[3, i, j, k] = x - 0.25y + 4z + x * z
    end
    J_zero = zeros(Float64, 3, nx, ny, nz)
    A_rec_klu = Magnetostatics.solve(solver_klu, J_zero; boundary = A_harm)
    A_rec_dst = Magnetostatics.solve(solver_dst, J_zero; boundary = A_harm)
    @test isapprox(A_rec_klu, A_harm; atol = 1.0e-12, rtol = 1.0e-12)
    @test isapprox(A_rec_dst, A_harm; atol = 1.0e-12, rtol = 1.0e-12)

    # Multipole boundary accuracy check for a compact ring current
    N_m = 20
    xm = range(-1.0, 1.0; length = N_m)
    J_ring = zeros(Float64, 3, N_m, N_m, N_m)
    width = 2 * step(xm)
    R0 = 0.25
    for k in 1:N_m, j in 1:N_m, i in 1:N_m
        x, y, z = xm[i], xm[j], xm[k]
        rho = hypot(x, y)
        if rho > 1.0e-6
            w = exp(-((rho - R0)^2 + z^2) / width^2) / (π * width^2)
            J_ring[1, i, j, k] = (-y / rho) * w
            J_ring[2, i, j, k] = (x / rho) * w
        end
    end
    bc_mp = multipole_boundary(J_ring, xm, xm, xm; order = 2)
    @test size(bc_mp) == (3, N_m, N_m, N_m)
    @test maximum(abs, bc_mp[:, 2:(N_m - 1), 2:(N_m - 1), 2:(N_m - 1)]) == 0.0
    dst_solver = PoissonSolver(xm, xm, xm; alg = DSTFactorization())
    A_zero_bc = Magnetostatics.solve(dst_solver, J_ring)
    A_mp_bc = Magnetostatics.solve(dst_solver, J_ring; boundary = bc_mp)
    A_open = vector_potential(FFTSolver(OpenBoundary()), J_ring, xm, xm, xm)
    bmask = falses(N_m, N_m, N_m)
    bmask[1, :, :] .= true; bmask[N_m, :, :] .= true
    bmask[:, 1, :] .= true; bmask[:, N_m, :] .= true
    bmask[:, :, 1] .= true; bmask[:, :, N_m] .= true
    bnd_err_zero = maximum(abs, A_open[:, bmask])
    bnd_err_mp = maximum(abs, bc_mp[:, bmask] - A_open[:, bmask])
    @test bnd_err_mp < 0.05 * bnd_err_zero
    @test maximum(abs, A_mp_bc - A_open) < 0.35 * maximum(abs, A_zero_bc - A_open)

    # Off-center ring confirms order = 2 quadrupole beats order = 1 dipole
    J_off = zeros(Float64, 3, N_m, N_m, N_m)
    for k in 1:N_m, j in 1:N_m, i in 1:N_m
        x, y, z = xm[i] - 0.15, xm[j] + 0.1, xm[k] - 0.12
        rho = hypot(x, y)
        if rho > 1.0e-6
            w = exp(-((rho - 0.2)^2 + z^2) / width^2) / (π * width^2)
            J_off[1, i, j, k] = (-y / rho) * w
            J_off[2, i, j, k] = (x / rho) * w
        end
    end
    bc_o1 = multipole_boundary(J_off, xm, xm, xm; order = 1)
    bc_o2 = multipole_boundary(J_off, xm, xm, xm; order = 2)
    A_off_open = vector_potential(FFTSolver(OpenBoundary()), J_off, xm, xm, xm)
    @test maximum(abs, bc_o2[:, bmask] - A_off_open[:, bmask]) <
        0.5 * maximum(abs, bc_o1[:, bmask] - A_off_open[:, bmask])
end

end # module test_poisson
