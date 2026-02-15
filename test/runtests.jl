using Test
using Magnetostatics
using LinearAlgebra
using StaticArrays

@testset "Magnetostatics.jl" begin
    include("test_biot_savart.jl")
    include("test_analytical.jl")
    include("test_fft.jl")
    include("test_potential.jl")
    include("test_utils.jl")
end
