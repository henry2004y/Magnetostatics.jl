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

# Analytical
SUITE["solvers"]["analytical"] = BenchmarkGroup()
r, loop = let
    r = SVector(0.5, 0.0, 0.1)
    # radius=3, current=1
    loop = CurrentLoop(3.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))

    r, loop
end

SUITE["solvers"]["analytical"]["current_loop"] = @benchmarkable getB_loop($r, $loop)
