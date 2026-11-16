using SMLMio, Test

# Docs.hasdoc is Julia 1.11+; a docstring is an entry for the name's binding in its module's docs.
_hasdoc(m, s) = (b = Base.Docs.Binding(m, s); haskey(Base.Docs.meta(b.mod), b))

@testset "every exported name has a docstring" begin
    # ?name and the API page need each exported name's docstring attached to it.
    undocumented = [n for n in names(SMLMio) if n !== :SMLMio && !_hasdoc(SMLMio, n)]
    @test undocumented == Symbol[]
end

# The README is the docs Home's front door: the same quick start, and the same "Where it fits".
@testset "README matches the docs Home" begin
    root = dirname(@__DIR__)
    readme = read(joinpath(root, "README.md"), String)
    home = replace(read(joinpath(root, "docs", "src", "index.md"), String),
                   "```@example quickstart" => "```julia")
    section(text, name) = match(Regex("^## $name\\n(.*?)(?=^## |\\z)", "ms"), text).captures[1]
    # Home's quick start is followed by links to the other pages
    through_code(s) = s[1:last(findlast("```\n", s))]
    @test through_code(section(readme, "Quick start")) == through_code(section(home, "Quick start"))
    @test section(readme, "Where it fits") == section(home, "Where it fits")
    # The examples use SMLMCore's types, so Installation adds SMLMCore
    for text in (readme, home)
        @test occursin("Pkg.add(\"SMLMCore\")", section(text, "Installation"))
    end
end
