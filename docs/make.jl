using Documenter
using FisherWright

makedocs(
    sitename = "FisherWright.jl",
    authors = "Xijiang Yu",
    modules = [FisherWright],
    checkdocs = :exports,
    format = Documenter.HTML(),
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

deploydocs(
    repo = "github.com/JuliaBnG/FisherWright.jl.git",
    deploy_repo = "github.com/JuliaBnG/juliabng.github.io.git",
    dirname = "FisherWright",
)
