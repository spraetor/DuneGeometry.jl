using Documenter
using DuneGeometry

push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "DuneGeometry.jl Documentation",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Index" => "index.md",
    ],
    modules = [DuneGeometry],
    doctest = false,
)

deploydocs(
    repo = "github.com/spraetor/DuneGeometry.jl.git",
    devbranch = "main"
)
