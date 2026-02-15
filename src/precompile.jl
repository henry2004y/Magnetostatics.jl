using PrecompileTools

@setup_workload begin
    let wire = Wire([SVector(0.0, 0.0, -1.0), SVector(0.0, 0.0, 1.0)], 1.0),
            loop = CurrentLoop(1.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0)),
            harris = HarrisSheet(1.0, 1.0),
            dipole = Dipole(SVector(0.0, 1.0, 0.0)),
            loop_analytic = CurrentLoopAnalytic(loop),
            r = SVector(0.5, 0.0, 0.0),
            points = [SVector(0.1, 0.1, 0.1)],
            bs_solver = BiotSavart(),
            vp_solver = VectorPotential(),
            fft_solver = FFTSolver(),
            Nx = 8, Ny = 8, Nz = 8,
            dx = 0.1,
            x_grid = collect(range(-1, 1, length = 10)),
            y_grid = collect(range(-1, 1, length = 10)),
            z_grid = collect(range(-1, 1, length = 10))

        J_fft = zeros(3, Nx, Ny, Nz)
        J_fft[3, :, :, :] .= 1.0

        @compile_workload begin
            # BiotSavart
            solve(bs_solver, wire, r)
            for p in points
                solve(bs_solver, wire, p)
            end

            # VectorPotential
            solve(vp_solver, wire, r)

            # Analytical Fields
            harris(r)
            dipole(r)
            loop_analytic(r)
            getB_loop(r, loop)
            getB_loop(r, SVector(0.0, 0.0, 0.0), 1.0, 1.0, SVector(0.0, 0.0, 1.0))

            # Utilities
            getB_mirror(0.5, 0.0, 0.0, 1.0, 1.0, 1.0)
            getB_bottle(0.5, 0.0, 0.0, 1.0, 1.0, 0.5, 1.0, 1.0)
            getB_tokamak_coil(0.5, 0.0, 0.0, 1.0, 2.0, 1.0, 1.0)
            getB_zpinch(0.5, 0.0, 0.0, 1.0, 0.1)

            # FFT Solver
            solve(fft_solver, J_fft, dx)

            # Helper set_current_wire!
            J_wire = zeros(3, 10, 10, 10)
            set_current_wire!(
                J_wire, x_grid, y_grid, z_grid,
                SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0), 1.0, 0.1
            )
            set_current_wire(
                x_grid, y_grid, z_grid,
                SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0), 1.0, 0.1
            )

            # Discretize loop
            discretize_loop(1.0, 10, 1.0)
            discretize_loop(loop, 10)
        end
    end
end
