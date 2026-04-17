using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

@testset "Code quality (Aqua.jl)" begin
    using Aqua
    Aqua.test_all(Magnetostatics)
end

@testset "Magnetostatics.jl" begin
    include("test_biot_savart.jl")
    include("test_analytical.jl")
    include("test_fft.jl")
    include("test_potential.jl")
    include("test_utils.jl")
    @testset "knot" include("test_knot.jl")
    @testset "boundaries" include("test_boundaries.jl")
    @testset "Poisson" include("test_poisson.jl")
end
