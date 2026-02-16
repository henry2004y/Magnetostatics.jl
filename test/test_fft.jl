using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra
using FFTW

@testset "FFT Solver" begin
    # Test Setup: Current Loop
    # We will discretize a loop onto a grid and compare the field at the center
    # with the analytical result.

    Nx, Ny, Nz = 64, 64, 64
    dim = 2.0 # Physical size [-1, 1]
    dx = dim / Nx

    x = range(-dim / 2, dim / 2 - dx, length = Nx)
    y = range(-dim / 2, dim / 2 - dx, length = Ny)
    z = range(-dim / 2, dim / 2 - dx, length = Nz)

    J_grid = zeros(Float64, 3, Nx, Ny, Nz)

    # Define a loop
    R = 0.5
    I = 1.0
    width = 2 * dx # Gaussian width for soft current distribution

    # Deposit current onto grid (Gaussian limits aliasing)
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        xx, yy, zz = x[i], y[j], z[k]

        # Distance from loop wire
        # Loop in xy plane, radius R
        rho = sqrt(xx^2 + yy^2)
        dist_sq = (rho - R)^2 + zz^2

        # Profile factor
        val = exp(-dist_sq / width^2)

        # Direction of current: azimuthal phi_hat = (-y, x, 0) / rho
        if rho > 1.0e-10
            Jx = -yy / rho * I
            Jy = xx / rho * I
            Jz = 0.0

            # Normalize so that integrated current is I?
            # Actually, standard Gaussian deposition: I_density ~ I * val
            # Current density J magnitude should integrate over cross section to I.
            # Cross section integral of exp(-(r-R)^2/w^2 - z^2/w^2) dA
            # is pi * w^2.
            # So J0 = I / (pi * w^2)

            J0 = I / (π * width^2)

            J_grid[1, i, j, k] = Jx * J0 * val
            J_grid[2, i, j, k] = Jy * J0 * val
            J_grid[3, i, j, k] = Jz * J0 * val
        end
    end

    solver = FFTSolver()
    B_grid = solve(solver, J_grid, dx)

    # Check field at center (0,0,0)
    # Index for 0 is Nx/2 + 1
    c_i, c_j, c_k = Nx ÷ 2 + 1, Ny ÷ 2 + 1, Nz ÷ 2 + 1

    Bx_c = B_grid[1, c_i, c_j, c_k]
    By_c = B_grid[2, c_i, c_j, c_k]
    Bz_c = B_grid[3, c_i, c_j, c_k]

    # Analytical verification
    Bz_analytical = Magnetostatics.μ₀ * I / (2 * R)

    # FFT result will differ slightly due to:
    # 1. Finite grid resolution
    # 2. Gaussian spreading of current (diffuses the source)
    # 3. Periodic boundary conditions (images of the loop affect the center)

    # Being within 5-10% is usually expected for a naive test setup like this without padding
    @test isapprox(Bz_c, Bz_analytical, rtol = 0.1) &&
        abs(Bx_c) < 1.0e-8 && abs(By_c) < 1.0e-8
end
