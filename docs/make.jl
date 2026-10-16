using Documenter, DocumenterVitepress
using Magnetostatics

makedocs(;
    sitename = "Magnetostatics.jl",
    modules = [Magnetostatics],
    authors = "Hongyang Zhou <hyzhou@umich.edu> and contributors",
    repo = "https://github.com/henry2004y/Magnetostatics.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/henry2004y/Magnetostatics.jl",
        devbranch = "main",
        devurl = "dev"
    ),
    pages = [
        "Home" => "index.md",
        "Forward Methods" => "forward.md",
        "Inverse Methods" => [
            "Choosing a Method" => "inverse/index.md",
            "Axisymmetric Reconstruction" => "inverse/axisymmetric.md",
            "3D Reconstruction" => "inverse/three_dimensional.md",
            "Spacecraft Grad-Shafranov" => "inverse/grad_shafranov.md",
        ],
        "Forward Examples" => [
            "Sources and Geometry" => [
                "Magnetic Dipole" => "examples/dipole.md",
                "Current Loop" => "examples/current_loop.md",
                "Trefoil Knot" => "examples/trefoil_knot.md",
                "Magnetic Mirror" => "examples/magnetic_mirror.md",
            ],
            "Plasma and Laboratory Fields" => [
                "Tokamak" => "examples/tokamak.md",
                "Solov'ev Equilibrium" => "examples/solovev.md",
                "Taylor State" => "examples/taylor_state.md",
                "Current Sheet" => "examples/current_sheet.md",
                "FLARE Experiment and Probe Fitting" => "examples/flare.md",
                "Z-Pinch" => "examples/z_pinch.md",
            ],
            "Space and Astrophysical Fields" => [
                "Parker Spiral" => "examples/parker_spiral.md",
                "Magnetosphere" => "examples/magnetosphere.md",
                "Empirical Magnetosphere" => "examples/empirical_magnetosphere.md",
            ],
            "Numerical Solvers" => [
                "Biot-Savart Solver" => "examples/biot_savart.md",
                "FFT Solver" => "examples/fft.md",
                "Poisson Solver" => "examples/poisson.md",
            ],
        ],
        "Inverse Examples" => [
            "3D Field Reconstruction" => "examples/reconstruction3d.md",
            "Reconstruction Cases" => "examples/reconstruction_cases.md",
            "Uncertainty and Hidden Topology" => "examples/reconstruction_uncertainty.md",
            "Spacecraft GSR" => "examples/grad_shafranov.md",
        ],
        "Field Geometry and Diagnostics" => "examples/magnetic_geometry.md",
        "API Reference" => "api.md",
    ],
    warnonly = [:missing_docs, :linkcheck],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/henry2004y/Magnetostatics.jl",
    target = "build",
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true
)
