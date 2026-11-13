using Documenter
using FisherWright

# JuliaBnG Sepia and Black themes, added to Documenter's theme picker.
# The files in src/assets/juliabng/ are built by docs-themes/build.jl in the
# juliabng.github.io repository; rebuild there rather than editing them here.
const JULIABNG_ASSETS = [
    [asset("assets/juliabng/$t.css"; class = :css, islocal = true,
           attributes = Dict(Symbol("data-theme-name") => t))
     for t in ("juliabng-sepia", "juliabng-black")];
    "assets/juliabng/themes.js"
]

makedocs(
    sitename = "FisherWright.jl",
    authors = "Xijiang Yu",
    modules = [FisherWright],
    checkdocs = :exports,
    format = Documenter.HTML(assets = JULIABNG_ASSETS),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Simulation" => "manual/simulation.md",
            "Results and export" => "manual/results.md",
            "Utilities" => "manual/utilities.md",
            "Population genetics validation" => "manual/validation.md",
            "Validation: changing population size" => "manual/validation-demography.md",
            "Validation: selection" => "manual/validation-selection.md",
            "Validation: mating schemes" => "manual/validation-mating.md",
        ],
        "API reference" => "api.md",
        "API reference: population models" => "api-models.md",
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
