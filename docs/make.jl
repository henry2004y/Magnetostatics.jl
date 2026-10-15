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
        "Examples" => [
            "Magnetic Dipole" => "examples/dipole.md",
            "Current Loop" => "examples/current_loop.md",
            "Trefoil Knot" => "examples/trefoil_knot.md",
            "Magnetic Mirror" => "examples/magnetic_mirror.md",
            "Magnetic Geometry" => "examples/magnetic_geometry.md",
            "Tokamak" => "examples/tokamak.md",
            "Solov'ev Equilibrium" => "examples/solovev.md",
            "Taylor State" => "examples/taylor_state.md",
            "Current Sheet" => "examples/current_sheet.md",
            "FLARE Experiment" => "examples/flare.md",
            "3D Field Reconstruction" => "examples/reconstruction3d.md",
            "Parker Spiral" => "examples/parker_spiral.md",
            "Magnetosphere" => "examples/magnetosphere.md",

            "Empirical Magnetosphere" => "examples/empirical_magnetosphere.md",
            "Z-Pinch" => "examples/z_pinch.md",
            "Biot-Savart Solver" => "examples/biot_savart.md",
            "FFT Solver" => "examples/fft.md",
            "Poisson Solver" => "examples/poisson.md",
        ],
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
