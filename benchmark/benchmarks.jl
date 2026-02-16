using BenchmarkTools
using Magnetostatics
using StaticArrays

const SUITE = BenchmarkGroup()

SUITE["solvers"] = BenchmarkGroup()

# Biot-Savart
SUITE["solvers"]["biot_savart"] = BenchmarkGroup()
let
    r = SVector(0.5, 0.5, 0.0)
    # Wire from (0,0,-1) to (0,0,1) with current 1.0
    wire = Wire([SVector(0.0, 0.0, -1.0), SVector(0.0, 0.0, 1.0)], 1.0)
    solver = BiotSavart()

    SUITE["solvers"]["biot_savart"]["simple_wire"] = @benchmarkable solve($solver, $wire, $r)
end

# Analytical
SUITE["solvers"]["analytical"] = BenchmarkGroup()
let
    r = SVector(0.5, 0.0, 0.1)
    loop = CurrentLoop(3.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0)) # radius=3, current=1

    SUITE["solvers"]["analytical"]["current_loop"] = @benchmarkable getB_loop($r, $loop)
end
