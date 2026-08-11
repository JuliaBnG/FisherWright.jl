using Pkg

Pkg.develop(path = joinpath(@__DIR__, ".."))

using Documenter
using FisherWright

makedocs(
    sitename = "FisherWright.jl",
    authors = "Xijiang Yu",
    modules = [FisherWright],
    checkdocs = :exports,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://xijiang.org/JuliaBnG/FisherWright",
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Simulation" => "manual/simulation.md",
            "Results and export" => "manual/results.md",
            "Utilities" => "manual/utilities.md",
        ],
        "API reference" => "api.md",
    ],
)
