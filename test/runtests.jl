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
    include("test_current_density.jl")
    @testset "geometry" include("test_geometry.jl")
    @testset "equilibria" include("test_equilibria.jl")
    @testset "astrophysical" include("test_astrophysical.jl")
    @testset "flare" include("test_flare.jl")
    @testset "knot" include("test_knot.jl")
    @testset "boundaries" include("test_boundaries.jl")
    @testset "Poisson" include("test_poisson.jl")
    @testset "reconstruction" include("test_reconstruction.jl")
    @testset "3D reconstruction" include("test_reconstruction3d.jl")
    @testset "Physical reconstruction cases" include("test_reconstruction_cases.jl")
    @testset "Reconstruction uncertainty" include("test_reconstruction_uncertainty.jl")
    @testset "Spacecraft Grad-Shafranov reconstruction" include("test_grad_shafranov.jl")
    @testset "FieldTracer integration" include("test_fieldtracer.jl")
end
