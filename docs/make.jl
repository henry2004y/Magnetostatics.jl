using Documenter
using Magnetostatics

makedocs(;
    sitename = "Magnetostatics.jl",
    modules = [Magnetostatics],
    remotes = nothing,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)
