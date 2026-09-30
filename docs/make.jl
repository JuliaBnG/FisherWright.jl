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
            "Population genetics validation" => "manual/validation.md",
        ],
        "API reference" => "api.md",
    ],
)

if !isempty(get(ENV, "DOCUMENTER_KEY", ""))
    deploydocs(
        repo = "github.com/JuliaBnG/FisherWright.jl.git",
        deploy_repo = "github.com/JuliaBnG/juliabng.github.io.git",
        dirname = "FisherWright",
        forcepush = true,
    )
elseif get(ENV, "GITHUB_ACTIONS", "") == "true"
    @warn "Skipping documentation deployment because DOCUMENTER_KEY is not configured."
end
