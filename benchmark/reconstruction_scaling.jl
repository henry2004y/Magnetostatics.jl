using BenchmarkTools
using Magnetostatics
using StaticArrays
using LinearAlgebra

point_matrix(points) = [p[k] for p in vec(collect(points)), k in 1:3]
lattice(coordinates) = point_matrix(Iterators.product(ntuple(_ -> coordinates, 3)...))
sample_field(field, positions) =
    point_matrix(field(SVector{3}(p)) for p in eachrow(positions))

function reconstruction_scaling()
    domain = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
    truth = CurrentLoopAnalytic(
        CurrentLoop(
            1.0, 1.0e4, SVector(0.2, -0.1, 2.0), normalize(SVector(0.3, 0.2, 1.0))
        )
    )
    heldout = lattice((-0.37, -0.13, 0.11, 0.39))
    expected = sample_field(truth, heldout)
    point = SVector(0.11, -0.17, 0.23)
    println("# Julia $(VERSION), threads=$(Threads.nthreads()), CPU=$(Sys.CPU_NAME)")
    println(
        "kernel,centers,probes,length_scale,regularization,relative_B_error," *
            "relative_curl_error,support_fraction,iterations,fit_ms,fit_bytes," *
            "B_ns,B_bytes,J_ns,J_bytes"
    )
    for kernel in (:gaussian, :wendland_c2), n in (3, 4, 5), probes in (3, 5),
            (scale, penalty) in ((0.6, 1.0e-8), (1.2, 1.0e-8), (1.2, 1.0e8))
        positions = lattice(range(-0.45, 0.45; length = probes))
        observed = sample_field(truth, positions)
        observed .+= 1.0e-5 .* reshape(sin.(1:length(observed)), size(observed))
        kwargs = (;
            domain, length_scale = scale, kernel,
            basis_centers = lattice(range(-0.5, 0.5; length = n)),
            uncertainties = fill(1.0e-5, size(observed)), regularization = penalty,
            regularization_grid = (5, 5, 5), abstol = 1.0e-12, reltol = 1.0e-12,
        )
        field = reconstruct_field_3d(positions, observed; kwargs...)
        error = norm(sample_field(field, heldout) - expected) / norm(expected)
        currents = point_matrix(current_density(field, p) for p in eachrow(heldout))
        curl_error = Magnetostatics.μ₀ * scale * norm(currents) / norm(expected)
        fit = minimum(
            @benchmark(
                reconstruct_field_3d($positions, $observed; $kwargs...),
                samples = 3, evals = 1, seconds = 0.2
            )
        )
        B = minimum(@benchmark $field($point) seconds = 0.2)
        J = minimum(@benchmark current_density($field, $point) seconds = 0.2)
        println(
            join(
                (
                    kernel, n^3, probes^3, scale, penalty, error, curl_error,
                    field.diagnostics.sampled_support_fraction, field.diagnostics.iterations,
                    fit.time / 1.0e6, fit.memory, B.time, B.memory, J.time, J.memory,
                ), ','
            )
        )
    end
    return
end

reconstruction_scaling()
