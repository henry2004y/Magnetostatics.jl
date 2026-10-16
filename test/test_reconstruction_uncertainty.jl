module test_reconstruction_uncertainty

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

include(joinpath(
    @__DIR__, "..", "docs", "src", "examples", "reconstruction_uncertainty.jl"
))
using .ReconstructionUncertaintyDemo
const Demo = ReconstructionUncertaintyDemo

@testset "Conditional noise propagation and boundary information" begin
    experiment = Demo.run_ensemble(; pairs = 64)
    data = experiment.data
    truth_values = Demo.samples(data.truth, data.heldout)
    for result in (experiment.free, experiment.bounded)
        @test length(result.fields) == 128
        @test size(result.summary.average) == size(truth_values)
        @test all(isfinite, result.summary.variance)
        @test all(result.summary.variance .>= 0)
        @test result.summary.average ≈
            Demo.samples(result.nominal, data.heldout) atol = 1.0e-12
        @test norm(result.summary.average - truth_values) / norm(truth_values) < 0.1
        @test all(n.status == :converged for n in result.nulls)
        @test all(n.residual < 1.0e-11 && n.positive == 1 for n in result.nulls)
        @test maximum(norm(n.point) for n in result.nulls) < 0.15
        refined = Demo.locate_null(result.nominal; h = 5.0e-6)
        coarse = Demo.locate_null(result.nominal)
        @test refined.status == coarse.status == :converged
        @test norm(refined.point - coarse.point) < 1.0e-8
    end
    @test experiment.bounded.summary.rms_spread < experiment.free.summary.rms_spread
    null_variance(nulls) = let
        points = Demo.point_matrix(n.point for n in nulls)
        center = sum(points; dims = 1) / length(nulls)
        sum(abs2, points .- center) / (length(nulls) - 1)
    end
    @test null_variance(experiment.bounded.nulls) < null_variance(experiment.free.nulls)

    point = SVector(0.13, 0.17, -0.23)
    covariances = map((false, true)) do constrained
        result = constrained ? experiment.bounded : experiment.free
        covariance = Demo.prediction_covariance(data, point; experiment.sigma, constrained)
        draws = hcat((field(point) for field in result.fields)...)
        centered = draws .- sum(draws; dims = 2) / size(draws, 2)
        empirical = centered * centered' / size(draws, 2)
        @test issymmetric(covariance)
        @test minimum(eigvals(Symmetric(covariance))) > 0
        @test norm(empirical - covariance) / norm(covariance) < 0.3
        covariance
    end
    @test tr(covariances[2]) < tr(covariances[1])
    @test maximum(eigvals(Symmetric(covariances[2] - covariances[1]))) < 1.0e-18

    doubled = Demo.run_ensemble(; pairs = 64, sigma = 2experiment.sigma)
    for name in (:free, :bounded)
        original = getproperty(experiment, name)
        scaled = getproperty(doubled, name)
        @test scaled.summary.variance ≈ 4original.summary.variance rtol = 1.0e-7
    end
    @test Demo.locate_null(NullField()).status == :degenerate
    @test_throws ArgumentError Demo.run_ensemble(; pairs = 1)
    @test_throws ArgumentError Demo.run_ensemble(; sigma = 0.0)
end

@testset "Invisible current can change null topology" begin
    case = Demo.ambiguity_case()
    data = case.data
    @test Demo.samples(case.alternative, data.positions) ==
        Demo.samples(data.truth, data.positions)
    @test all(
        dot(n, case.alternative(p)) == value
            for (p, n, value) in zip(
                eachrow(data.surface), eachrow(data.normals), data.boundary
            )
    )
    @test maximum(abs, case.hidden.centers) + case.hidden.length_scale < 1
    @test minimum(
        norm(p - c) for p in eachrow(data.positions), c in eachrow(case.hidden.centers)
    ) > case.hidden.length_scale
    baseline = Demo.locate_null(data.truth)
    changed = Demo.locate_null(case.alternative)
    @test baseline.status == changed.status == :converged
    @test baseline.positive == 1
    @test changed.positive == 2
    for root in (case.positive_root, case.negative_root)
        @test root.status == :converged
        @test root.residual < 1.0e-11
        @test 0.15 < norm(root.point) < 0.3
        @test norm(data.truth(root.point)) > 1.0e-5
        refined = Demo.locate_null(case.alternative; start = root.point, h = 5.0e-6)
        @test refined.status == :converged
        @test refined.point ≈ root.point atol = 1.0e-8
    end
    @test norm(case.positive_root.point - case.negative_root.point) > 0.3
    point = SVector(0.1, 0.02, 0.01)
    conditional = Demo.prediction_covariance(
        data, point; sigma = 1.0e-5, constrained = true
    )
    @test norm(case.alternative(point) - data.truth(point)) > 10sqrt(tr(conditional))
    @test norm(current_density(case.hidden, point)) > 1.0
    for field in (data.truth, case.alternative), p in (
            SVector(0.07, -0.03, 0.02), SVector(0.18, 0.05, -0.06),
        )
        @test abs(tr(Demo.jacobian(field, p))) < 1.0e-8
    end
end

end
