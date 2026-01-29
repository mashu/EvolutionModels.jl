using EvolutionModels
using Documenter

DocMeta.setdocmeta!(EvolutionModels, :DocTestSetup, :(using EvolutionModels); recursive=true)

makedocs(;
    modules=[EvolutionModels],
    authors="Mateusz Kaduk <mateusz.kaduk@gmail.com> and contributors",
    sitename="EvolutionModels.jl",
    format=Documenter.HTML(;
        canonical="https://mashu.github.io/EvolutionModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theoretical Guide" => "guide.md",
    ],
)

deploydocs(;
    repo="github.com/mashu/EvolutionModels.jl",
    devbranch="main",
)
