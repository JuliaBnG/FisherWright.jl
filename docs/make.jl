using Documenter, FisherWright
makedocs(sitename = "FisherWright.jl", modules = [FisherWright], format = Documenter.HTML())
deploydocs(repo = "github.com/JuliaBnG/FisherWright.jl.git")
