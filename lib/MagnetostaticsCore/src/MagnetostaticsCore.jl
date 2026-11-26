module MagnetostaticsCore

using LinearAlgebra: cross, det, dot, I, normalize
using PrecompileTools: @compile_workload, @setup_workload
using StaticArrays: @SMatrix, SMatrix, SVector

const μ₀ = 4π * 1.0e-7
const μ0_4π = 1.0e-7

include("types.jl")
include("elliptic.jl")
include("sources.jl")
include("composition.jl")
include("fields/uniform.jl")
include("fields/dipole.jl")
include("fields/wire.jl")
include("fields/spherical_harmonic.jl")
include("precompile.jl")

export AbstractMagneticField, AbstractCurrentSource, AbstractSolver
export CompositeField, ScaledField, TransformedField, translate_field, rotate_field
export CurrentLoop
export UniformField, NullField, Dipole, LineDipole2D
export InfiniteWire, WireSegment, CurrentLoopAnalytic, getB_loop
export ToroidalField, SphericalHarmonicField, SphericalHarmonicWorkspace
export current_density, vector_potential, field_and_potential, evaluate!, poloidal_flux

end
