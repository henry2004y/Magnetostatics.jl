using Documenter
using Magnetostatics

makedocs(;
    sitename = "Magnetostatics.jl",
    modules = [Magnetostatics],
    repo = Remotes.GitHub("henry2004y", "Magnetostatics.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)
