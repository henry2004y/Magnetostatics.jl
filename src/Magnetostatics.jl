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

# Core types
export AbstractMagneticField, AbstractCurrentSource, AbstractSolver, AbstractBoundary
# Solvers
export BiotSavart, FFTSolver, VectorPotential, PoissonSolver, solve
# Core methods
export vector_potential, current_density, compute_curl!, compute_curl
# Source geometries
export Wire, CurrentLoop, TorusKnot, TrefoilKnot, InfiniteWire, SurfaceCurrentMesh
# Analytical fields
export HarrisSheet, AsymmetricHarrisSheet, ForceFreeHarrisSheet, BifurcatedHarrisSheet,
    FadeevIsland, Dipole, CurrentLoopAnalytic, UniformField, NullField
# Boundary conditions
export OpenBoundary, PeriodicBoundary, ConductingWall, ConductingSphere, PrescribedBoundary
# Specialized models and configurations
export AnalyticalMagnetosphere, draped_imf_field, getB_mirror, getB_bottle,
    getB_tokamak_coil, getB_tokamak_profile, getB_zpinch
# Utilities
export discretize_loop, discretize_knot, getB_loop, set_current_wire!, set_current_wire,
    compute_surface_current, image_source, dipole_fieldline, sph2cart, sphere_mesh

include("precompile.jl")

end
