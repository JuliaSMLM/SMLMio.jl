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
        edit_link="main", # skip asking the remote for its default branch
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Guide" => ["SMITE files" => "guide/smite.md"],
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSMLM/SMLMio.jl.git",
    devbranch="main",
)
