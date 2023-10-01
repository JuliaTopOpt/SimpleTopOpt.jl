using SimpleTopOpt
using Documenter

DocMeta.setdocmeta!(SimpleTopOpt, :DocTestSetup, :(using SimpleTopOpt); recursive = true)

makedocs(;
    modules = [SimpleTopOpt],
    authors = "Matthew Meeker <mmmeeeker@gmail.com> and contributors",
    repo = "https://github.com/juliatopopt/SimpleTopOpt.jl/blob/{commit}{path}#{line}",
    sitename = "SimpleTopOpt.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliatopopt.github.io/SimpleTopOpt.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/juliatopopt/SimpleTopOpt.jl", devbranch = "main")
