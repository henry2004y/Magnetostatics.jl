using Documenter, DocumenterVitepress
using Magnetostatics

makedocs(;
    sitename = "Magnetostatics.jl",
    modules = [Magnetostatics],
    authors = "Hongyang Zhou <hyzhou@umich.edu> and contributors",
    repo = "https://github.com/henry2004y/Magnetostatics.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/henry2004y/Magnetostatics.jl",
        devbranch = "master",
        devurl = "dev"
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    warnonly = [:missing_docs, :linkcheck],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/henry2004y/Magnetostatics.jl",
    target = "build",
    devbranch = "master",
    branch = "gh-pages",
    push_preview = true
)
