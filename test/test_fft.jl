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

    # Test OpenBoundary Hockney solver
    open_solver = FFTSolver(OpenBoundary())
    B_open = solve(open_solver, J_grid, dx)
    Bz_open_c = B_open[3, c_i, c_j, c_k]
    @test isapprox(Bz_open_c, Bz_analytical, rtol = 0.05)
    @test abs(B_open[1, c_i, c_j, c_k]) < 1.0e-8
    @test abs(B_open[2, c_i, c_j, c_k]) < 1.0e-8

    # Point-current discrete Biot-Savart exactness (open boundary)
    let N = 16, h = (0.1, 0.12, 0.15)
        J_pt = zeros(Float64, 3, N, N, N)
        src_ijk = (2, 4, 3)
        dV = prod(h)
        J_pt[1, src_ijk...] = 0.5 / dV
        J_pt[2, src_ijk...] = -0.3 / dV
        J_pt[3, src_ijk...] = 1.2 / dV
        I_vec = SVector(0.5, -0.3, 1.2)
        s_open = FFTSolver(OpenBoundary())
        B_pt = solve(s_open, J_pt, h)
        s_cached = FFTSolver(OpenBoundary(), (N, N, N), h)
        B_pt_cached = solve(s_cached, J_pt, h)
        @test isapprox(B_pt, B_pt_cached; atol = 1.0e-15, rtol = 1.0e-13)
        for k in (1, 14), j in (2, 15), i in (1, 13)
            dr = SVector(
                (i - src_ijk[1]) * h[1],
                (j - src_ijk[2]) * h[2],
                (k - src_ijk[3]) * h[3]
            )
            r = norm(dr)
            B_exact = (Magnetostatics.μ₀ / (4π)) * cross(I_vec, dr) / r^3
            @test isapprox(SVector{3}(B_pt[:, i, j, k]), B_exact; rtol = 1.0e-12)
        end
        A_pt = vector_potential(s_open, J_pt, h)
        A_pt_cached = vector_potential(s_cached, J_pt, h)
        @test isapprox(A_pt, A_pt_cached; atol = 1.0e-15, rtol = 1.0e-13)
        probe = (14, 15, 13)
        dr = SVector(
            (probe[1] - src_ijk[1]) * h[1],
            (probe[2] - src_ijk[2]) * h[2],
            (probe[3] - src_ijk[3]) * h[3]
        )
        A_exact = (Magnetostatics.μ₀ / (4π)) * I_vec / norm(dr)
        @test isapprox(SVector{3}(A_pt[:, probe...]), A_exact; rtol = 1.0e-12)
    end

    # Cached periodic solver, range dispatch, and Coulomb-gauge vector potential
    let Nx = 24, Ny = 20, Nz = 16
        xs = range(-1.0, 1.0; length = Nx + 1)[1:Nx]
        ys = range(-0.8, 0.8; length = Ny + 1)[1:Ny]
        zs = range(-0.6, 0.6; length = Nz + 1)[1:Nz]
        J_smooth = zeros(Float64, 3, Nx, Ny, Nz)
        kx0 = 2π / (Nx * step(xs))
        ky0 = 2π / (Ny * step(ys))
        kz0 = 2π / (Nz * step(zs))
        for k in 1:Nz, j in 1:Ny, i in 1:Nx
            s = sin(kx0 * xs[i] + ky0 * ys[j])
            c = cos(kz0 * zs[k])
            J_smooth[1, i, j, k] = -ky0 * cos(kx0 * xs[i] + ky0 * ys[j]) * c
            J_smooth[2, i, j, k] = kx0 * cos(kx0 * xs[i] + ky0 * ys[j]) * c
            J_smooth[3, i, j, k] = s * sin(kz0 * zs[k])
        end
        p_uncached = FFTSolver()
        p_cached = FFTSolver(
            PeriodicBoundary(SVector(0.0, 0.0, 0.0)),
            (Nx, Ny, Nz), (step(xs), step(ys), step(zs))
        )
        B_u = solve(p_uncached, J_smooth, xs, ys, zs)
        B_c = solve(p_cached, J_smooth, xs, ys, zs)
        @test isapprox(B_u, B_c; atol = 1.0e-14, rtol = 1.0e-13)
        A_u = vector_potential(p_uncached, J_smooth, xs, ys, zs)
        A_c = vector_potential(p_cached, J_smooth, xs, ys, zs)
        @test isapprox(A_u, A_c; atol = 1.0e-14, rtol = 1.0e-13)
        B_curl = compute_curl(A_u, xs, ys, zs)
        @test isapprox(
            B_curl[:, 3:(Nx - 2), 3:(Ny - 2), 3:(Nz - 2)],
            B_u[:, 3:(Nx - 2), 3:(Ny - 2), 3:(Nz - 2)];
            rtol = 5.0e-2
        )
    end
end
