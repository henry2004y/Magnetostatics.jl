module test_grad_shafranov

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

point_matrix(points) = [p[k] for p in vec(collect(points)), k in 1:3]

function manufactured_problem()
    ex, ez = normalize(SVector(1.0, 2.0, 1.0)), normalize(SVector(-1.0, 0.0, 1.0))
    frame = hcat(ex, cross(ez, ex), ez)
    origin = SVector(0.3, -0.2, 0.1)
    A(x, y) = 0.001x^2 - 0.0005y^2 + 0.0003x * y + 0.001x + 0.003y
    Bz(a) = 0.01 + 0.2a
    Pt(a) = 1000.0 - 0.001a / Magnetostatics.μ₀
    p(x, y) = Pt(A(x, y)) - Bz(A(x, y))^2 / (2Magnetostatics.μ₀)
    function local_field(x, y)
        return SVector(0.003 - 0.001y + 0.0003x, -0.002x - 0.001 - 0.0003y, Bz(A(x, y)))
    end
    times = collect(range(0.0, 1.0; length = 31))
    xs = 2 .* (times .- 0.5)
    observed = point_matrix(frame * local_field(x, 0.0) for x in xs)
    pressures = [p(x, 0.0) for x in xs]
    kwargs = (;
        crossing_velocity = 2ex + 3ez, invariant_axis = ez, origin,
        y = collect(range(-0.15, 0.15; length = 31)),
        spatial_degree = 2, pressure_degree = 1, axial_degree = 1,
    )
    return (; times, observed, pressures, kwargs, frame, origin, A, p, local_field)
end

@testset "Chebyshev derivatives" begin
    for u in (-1.0, -0.31, 0.0, 0.47, 1.0)
        @test Magnetostatics._chebyshev_sum([2.0], u) == (2.0, 0.0, 0.0)
        value, derivative, second = Magnetostatics._chebyshev_sum([0.0, 0.0, 0.0, 1.0], u)
        @test value ≈ 4u^3 - 3u atol = 1.0e-14
        @test derivative ≈ 12u^2 - 3 atol = 1.0e-14
        @test second ≈ 24u atol = 1.0e-14
    end
end

@testset "Planar equilibrium, frame and analytic currents" begin
    problem = manufactured_problem()
    field = @test_logs(
        (:warn, r"GSR source functions extrapolate"),
        reconstruct_grad_shafranov(
            problem.times, problem.observed, problem.pressures; problem.kwargs...
        )
    )
    @test field isa GradShafranovReconstruction
    @test field.diagnostics.track_field_rms < 1.0e-12
    @test field.diagnostics.transverse_pressure_rms < 1.0e-10
    @test field.diagnostics.axial_field_rms < 1.0e-12
    @test field.diagnostics.gs_residual_rms < 1.0e-9
    for x in (-1.0, -0.37, 0.13, 1.0), y in (-0.15, -0.023, 0.0, 0.071, 0.15)
        point = problem.origin + problem.frame * SVector(x, y, 0.0)
        expected = problem.frame * problem.local_field(x, y)
        @test field(point) ≈ expected atol = 1.0e-10
        @test field(point + 0.7problem.frame[:, 3]) ≈ expected atol = 1.0e-10
        expected_flux = problem.A(x, y) - problem.A(-1.0, 0.0)
        @test flux_function(field, point) ≈ expected_flux atol = 1.0e-10
        @test thermal_pressure(field, point) ≈ problem.p(x, y) atol = 1.0e-7
        Ax = 0.002x + 0.0003y + 0.001
        Ay = -0.001y + 0.0003x + 0.003
        expected_J = problem.frame * SVector(0.2Ay, -0.2Ax, -0.001) / Magnetostatics.μ₀
        @test current_density(field, point) ≈ expected_J atol = 1.0e-6
    end
    point = problem.origin + problem.frame * SVector(0.13, 0.071, 0.0)
    h = 1.0e-5
    gradient = hcat(
        [
            (
                field(point + h * SVector{3}(ntuple(i -> i == k, 3))) -
                    field(point - h * SVector{3}(ntuple(i -> i == k, 3)))
            ) / (2h)
                for k in 1:3
        ]...
    )
    curl = SVector(
        gradient[3, 2] - gradient[2, 3], gradient[1, 3] - gradient[3, 1],
        gradient[2, 1] - gradient[1, 2]
    )
    @test abs(tr(gradient)) < 1.0e-9
    @test current_density(field, point) ≈ curl / Magnetostatics.μ₀ atol = 1.0e-5
    pressure_gradient = SVector{3}(
        ntuple(3) do k
            step = h * SVector{3}(ntuple(i -> i == k, 3))
            plus, minus = thermal_pressure(field, point + step),
                thermal_pressure(field, point - step)
            (plus - minus) / (2h)
        end
    )
    force = cross(current_density(field, point), field(point))
    @test force ≈ pressure_gradient atol = 1.0e-6
    @test_throws DomainError field(problem.origin + 1.01problem.frame[:, 1])
    @test_throws DomainError thermal_pressure(field, problem.origin + problem.frame[:, 2])
    @test_throws ArgumentError field(SVector(NaN, 0.0, 0.0))
    reversed = @test_logs(
        (:warn, r"GSR source functions extrapolate"),
        reconstruct_grad_shafranov(
            problem.times, reverse(problem.observed; dims = 1), reverse(problem.pressures);
            merge(
                problem.kwargs, (;
                    crossing_velocity = -2problem.frame[:, 1] + 3problem.frame[:, 3],
                )
            )...
        )
    )
    @test reversed(point) ≈ field(point) atol = 1.0e-10
    shifted = @test_logs(
        (:warn, r"GSR source functions extrapolate"),
        reconstruct_grad_shafranov(
            problem.times .+ 1.0e6, problem.observed, problem.pressures; problem.kwargs...
        )
    )
    @test shifted(point) ≈ field(point) atol = 1.0e-10
    times = collect(range(0.0, 1.0; length = 31)) .^ 1.3
    xs = 2 .* (times .- 0.5)
    observed = point_matrix(problem.frame * problem.local_field(x, 0.0) for x in xs)
    pressures = [problem.p(x, 0.0) for x in xs]
    nonuniform = @test_logs(
        (:warn, r"GSR source functions extrapolate"),
        reconstruct_grad_shafranov(times, observed, pressures; problem.kwargs...)
    )
    @test nonuniform(point) ≈ field(point) atol = 1.0e-10
end

function island_problem(; noise = 0.0)
    B0, L, epsilon, guide, offset = 20.0e-9, 1.0e6, 0.3, 5.0e-9, 0.2e6
    island = FadeevIsland(B0, L, L, epsilon)
    function truth(point)
        field = island(SVector(point[1] + π * L, 0.0, point[2] + offset))
        return SVector(field[1], field[3], guide)
    end
    function pressure(point)
        f = cosh((point[2] + offset) / L) - epsilon * cos(point[1] / L)
        return 0.1e-9 + B0^2 * (1 - epsilon^2) / (2Magnetostatics.μ₀ * f^2)
    end
    times = collect(range(-40.0, 40.0; length = 121))
    positions = [SVector(50.0e3 * t, 0.0, 0.0) for t in times]
    observed = point_matrix(truth(p) for p in positions)
    observed .+= noise .* reshape(sin.(1:length(observed)), size(observed))
    pressures = [pressure(p) for p in positions]
    kwargs = (;
        crossing_velocity = SVector(50.0e3, 0.0, 0.0),
        invariant_axis = SVector(0.0, 0.0, 1.0),
        y = collect(range(-0.3e6, 0.3e6; length = 61)),
        spatial_degree = 10, pressure_degree = 6, axial_degree = 0,
    )
    return (; times, observed, pressures, kwargs, truth, pressure)
end

@testset "Mock satellite crossing a Fadeev island" begin
    for noise in (0.0, 0.02e-9)
        problem = island_problem(; noise)
        field = @test_logs(
            (:warn, r"GSR source functions extrapolate"),
            reconstruct_grad_shafranov(
                problem.times, problem.observed, problem.pressures; problem.kwargs...
            )
        )
        heldout = [
            SVector(x, y, 0.0) for x in range(-1.7e6, 1.7e6; length = 17),
                y in (-0.23e6, 0.17e6, 0.27e6)
        ]
        expected = point_matrix(problem.truth(p) for p in heldout)
        predicted = point_matrix(field(p) for p in heldout)
        @test norm(predicted - expected) / norm(expected) < 0.03
        pressures = [problem.pressure(p) for p in heldout]
        @test norm([thermal_pressure(field, p) for p in heldout] - pressures) /
            norm(pressures) < 0.03
        expected_J = point_matrix(
            SVector(
                0.0, 0.0, 20.0e-9 * (1 - 0.3^2) /
                    (
                    Magnetostatics.μ₀ * 1.0e6 *
                        (cosh((p[2] + 0.2e6) / 1.0e6) - 0.3cos(p[1] / 1.0e6))^2
                )
            ) for p in heldout
        )
        predicted_J = point_matrix(current_density(field, p) for p in heldout)
        @test norm(predicted_J - expected_J) / norm(expected_J) < 0.03
        refined = @test_logs(
            (:warn, r"GSR source functions extrapolate"),
            reconstruct_grad_shafranov(
                problem.times, problem.observed, problem.pressures;
                problem.kwargs..., max_step = 2000.0
            )
        )
        @test norm(point_matrix(refined(p) for p in heldout) - predicted) / norm(expected) <
            0.005
    end
end

@testset "Input validation and inconsistent pressure data" begin
    problem = manufactured_problem()
    fit(; kwargs...) = reconstruct_grad_shafranov(
        problem.times, problem.observed, problem.pressures;
        merge(problem.kwargs, (; kwargs...))...
    )
    for kwargs in (
            (; crossing_velocity = zeros(3)), (; crossing_velocity = problem.frame[:, 3]),
            (; invariant_axis = zeros(3)), (; invariant_axis = [1.0, 0.0]),
            (; origin = [0.0, Inf, 0.0]), (; y = [-0.1, 0.1]),
            (; y = (-0.1, 0.0, 0.1)),
            (; y = [0.0, 0.1, 0.2]), (; y = [-0.1, 0.01, 0.1]),
            (; y = [-0.1, 0.0, NaN]), (; spatial_degree = 0),
            (; spatial_degree = 30), (; pressure_degree = 0), (; axial_degree = -1),
            (; max_step = 0.0), (; max_step = NaN),
        )
        @test_throws ArgumentError fit(; kwargs...)
    end
    @test_throws ArgumentError reconstruct_grad_shafranov(
        reverse(problem.times), problem.observed, problem.pressures; problem.kwargs...
    )
    @test_throws ArgumentError reconstruct_grad_shafranov(
        problem.times, problem.observed[:, 1:2], problem.pressures; problem.kwargs...
    )
    @test_throws ArgumentError reconstruct_grad_shafranov(
        problem.times, fill(NaN, size(problem.observed)), problem.pressures;
        problem.kwargs...
    )
    @test_throws ArgumentError reconstruct_grad_shafranov(
        problem.times, problem.observed, -problem.pressures; problem.kwargs...
    )
    @test_throws ArgumentError reconstruct_grad_shafranov(
        problem.times, repeat([0.0 0.0 1.0], length(problem.times)), problem.pressures;
        merge(
            problem.kwargs, (;
                crossing_velocity = SVector(1.0, 0.0, 0.0),
                invariant_axis = SVector(0.0, 0.0, 1.0),
            )
        )...
    )
    noisy_pressure = problem.pressures + sin.(1:length(problem.pressures))
    field = @test_logs(
        (:warn, r"GSR source functions extrapolate"),
        reconstruct_grad_shafranov(
            problem.times, problem.observed, noisy_pressure; problem.kwargs...
        )
    )
    @test field.diagnostics.transverse_pressure_rms > 0.1
    @test_logs fit(; y = [-0.001, 0.0, 0.001])
    @test_logs(
        (:warn, r"GSR strip exceeds"), (:warn, r"GSR source functions extrapolate"),
        fit(; y = [-0.5, 0.0, 0.5])
    )
end

end # module test_grad_shafranov
