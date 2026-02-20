using BenchmarkTools
using Magnetostatics
using StaticArrays

const SUITE = BenchmarkGroup()

SUITE["solvers"] = BenchmarkGroup()

# Biot-Savart
SUITE["solvers"]["biot_savart"] = BenchmarkGroup()
solver, wire, r = let
    r = SVector(0.5, 0.5, 0.0)
    loop = CurrentLoop(1.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))
    # creates a 100-segment wire
    wire = discretize_loop(loop, 100)
    solver = BiotSavart()

    solver, wire, r
end

SUITE["solvers"]["biot_savart"]["discretized_loop"] = @benchmarkable solve($solver, $wire, $r)

# FFT Solver
SUITE["solvers"]["fft"] = BenchmarkGroup()
fft_solver, J, dx = let
    solver = FFTSolver()
    dx = 0.1
    # 3x32x32x32 array
    J = rand(3, 32, 32, 32)
    solver, J, dx
end

SUITE["solvers"]["fft"]["small_grid"] = @benchmarkable solve($fft_solver, $J, $dx)

# Vector Potential
SUITE["solvers"]["potential"] = BenchmarkGroup()
pot_solver, pot_r, pot_loop = let
    solver = VectorPotential()
    r = SVector(0.5, 0.0, 0.1)
    loop = CurrentLoop(1.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))
    solver, r, loop
end

SUITE["solvers"]["potential"]["current_loop"] = @benchmarkable solve($pot_solver, $pot_loop, $pot_r)

# Analytical
SUITE["solvers"]["analytical"] = BenchmarkGroup()
r, loop = let
    r = SVector(0.5, 0.0, 0.1)
    # radius=3, current=1
    loop = CurrentLoop(3.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))

    r, loop
end

SUITE["solvers"]["analytical"]["current_loop"] = @benchmarkable getB_loop($r, $loop)

r_dipole, dipole = let
    r = SVector(0.5, 0.0, 0.1)
    dipole = Dipole(SVector(0.0, 0.0, 1.0))
    r, dipole
end

SUITE["solvers"]["analytical"]["dipole"] = @benchmarkable ($dipole)($r_dipole)
