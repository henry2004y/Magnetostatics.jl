module Magnetostatics

using LinearAlgebra
using SparseArrays
using StaticArrays
using SpecialFunctions
using FFTW
using LinearSolve
using PrecompileTools

const μ₀ = 4π * 1.0e-7
const μ0_4π = 1.0e-7 # T*m/A (approximation for μ0/4π)

include("types.jl")
include("boundaries.jl")
include("sources.jl")
include("solvers/biot_savart.jl")
include("solvers/analytical.jl")
include("solvers/fft.jl")
include("solvers/potential.jl")
include("utils.jl")

export AbstractMagneticField, AbstractCurrentSource, AbstractSolver, BiotSavart, FFTSolver,
    VectorPotential, PoissonSolver, solve, vector_potential, Wire, CurrentLoop, TorusKnot, TrefoilKnot,
    HarrisSheet, Dipole, CurrentLoopAnalytic, InfiniteWire, UniformField,
    AsymmetricHarrisSheet, ForceFreeHarrisSheet, BifurcatedHarrisSheet, FadeevIsland,
    SurfaceCurrentMesh, compute_surface_current, current_density,
    AbstractBoundary, OpenBoundary, PeriodicBoundary, ConductingWall, ConductingSphere,
    PrescribedBoundary,
    discretize_loop, discretize_knot, getB_loop, set_current_wire!, set_current_wire,
    compute_curl!, compute_curl, getB_mirror, getB_bottle, getB_tokamak_coil,
    getB_tokamak_profile, getB_zpinch, image_source, dipole_fieldline, sph2cart, sphere_mesh

include("precompile.jl")

end
