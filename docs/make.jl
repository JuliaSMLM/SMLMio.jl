using SMLMio
using Documenter

makedocs(;
    modules=[SMLMio],
    sitename="SMLMio.jl",
    checkdocs=:exports,
    doctest=false, # QA runs every docstring's jldoctest, decision 0018.
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaSMLM.github.io/SMLMio.jl",
    ),
    pages=[
        "Home" => "index.md",
        "SMITE files" => "smite.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSMLM/SMLMio.jl.git",
    devbranch="main",
)
