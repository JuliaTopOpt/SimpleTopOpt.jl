using SimpleTopOpt
using Documenter

DocMeta.setdocmeta!(SimpleTopOpt, :DocTestSetup, :(using SimpleTopOpt); recursive=true)

makedocs(;
    modules=[SimpleTopOpt],
    authors="mjachi <mmmeeeker@gmail.com> and contributors",
    repo="https://github.com/mjachi/SimpleTopOpt.jl/blob/{commit}{path}#{line}",
    sitename="SimpleTopOpt.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mjachi.github.io/SimpleTopOpt.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mjachi/SimpleTopOpt.jl",
    devbranch="main",
)
