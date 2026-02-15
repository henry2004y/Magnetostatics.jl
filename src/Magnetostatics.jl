module Magnetostatics

using LinearAlgebra
using StaticArrays
using SpecialFunctions

const μ₀ = 4π * 1.0e-7
const μ0_4π = 1.0e-7 # T*m/A (approximation for μ0/4π)

include("types.jl")
include("sources.jl")
include("solvers/biot_savart.jl")
include("solvers/analytical.jl")
include("solvers/fft.jl")
include("solvers/potential.jl")
include("utils.jl")

export AbstractMagneticField, AbstractCurrentSource, AbstractSolver
export BiotSavart, FFTSolver, VectorPotential, solve
export Wire, CurrentLoop, HarrisSheet, Dipole, CurrentLoopAnalytic
export discretize_loop, getB_loop
export getB_mirror, getB_bottle, getB_tokamak_coil, getB_tokamak_profile, getB_zpinch

end
