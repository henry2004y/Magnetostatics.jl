module Magnetostatics

using LinearAlgebra
using SparseArrays
using StaticArrays
using SpecialFunctions
using FFTW
using LinearSolve
using PrecompileTools

using MagnetostaticsCore
import MagnetostaticsCore:
    μ₀, μ0_4π, _pos_eltype, _promote_float, _ellipke, _cel,
    _current_loop_unchecked, current_density, vector_potential,
    field_and_potential, evaluate!, poloidal_flux

include("boundaries.jl")
include("sources.jl")
include("solvers/biot_savart.jl")
include("solvers/analytical.jl")
include("solvers/fft.jl")
include("solvers/potential.jl")
include("solvers/reconstruction_axisymmetric.jl")
include("solvers/reconstruction3d.jl")
include("solvers/grad_shafranov.jl")
include("geometry.jl")
include("utils.jl")
include("grid_field.jl")
include("topology.jl")

# Core types
export AbstractMagneticField, AbstractCurrentSource, AbstractSolver, AbstractBoundary
export CompositeField, ScaledField, TransformedField, GridMagneticField
export translate_field, rotate_field
export field_jacobian, classify_magnetic_null, locate_magnetic_null,
    locate_poloidal_critical_point, field_covariance
# Solvers
export BiotSavart, FFTSolver, VectorPotential, PoissonSolver, DSTFactorization, solve
export multipole_boundary
export AxisymmetricFieldReconstruction, reconstruct_axisymmetric_field
export FieldReconstruction3D, reconstruct_field_3d
export GradShafranovReconstruction, reconstruct_grad_shafranov
export flux_function, thermal_pressure
# Core methods
export vector_potential, field_and_potential, evaluate!, current_density,
    compute_curl!, compute_curl
# Magnetic geometry and plasma drifts
export b_unit, field_magnitude, grad_b, magnetic_curvature, drift_velocities,
    mirror_ratio, loss_cone_angle, magnetic_helicity
# Source geometries
export Wire, CurrentLoop, TorusKnot, TrefoilKnot, InfiniteWire, WireSegment,
    SurfaceCurrentMesh, MultiSource
# Analytical fields
export HarrisSheet, AsymmetricHarrisSheet, ForceFreeHarrisSheet, BifurcatedHarrisSheet,
    FadeevIsland, Dipole, LineDipole2D, CurrentLoopAnalytic, ToroidalField,
    SphericalHarmonicField, SphericalHarmonicWorkspace, UniformField, NullField
# Boundary conditions
export OpenBoundary, PeriodicBoundary, ConductingWall, ConductingSphere, PrescribedBoundary
# Specialized models and configurations
export AnalyticalMagnetosphere, draped_imf_field,
    MagneticMirror, MagneticBottle, getB_mirror, getB_bottle,
    TokamakCoils, TokamakProfile, getB_tokamak_coil, getB_tokamak_profile,
    ZPinch, getB_zpinch,
    SolovevEquilibrium, poloidal_flux, TaylorState, ScrewPinch,
    ParkerSpiral, GuideFieldHarrisSheet, GoldHoyleFluxTube,
    FLARECoils, Solenoid, HelmholtzCoil, AntiHelmholtzCoil,
    SphereMagnet, CylinderMagnet, CuboidMagnet


# Utilities
export discretize_loop, discretize_knot, getB_loop, set_current_wire!, set_current_wire,
    compute_surface_current, image_source, dipole_fieldline, sph2cart, sphere_mesh


include("precompile.jl")

end
