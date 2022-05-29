using Documenter, MadDiff

const _PAGES = [
    "Home" => "index.md",
    "Quick Start"=>"guide.md",
    "How it works?" => "algorithms.md",
    "API Manual" => [
        "MadDiffCore" => "core.md",
        "MadDiffModels" => "models.md",
        "MadDiffMOI" => "moi.md",
    ]
]

makedocs(
    sitename = "MadDiff",
    modules = [MadDiff],
    authors = "Sungho Shin",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        sidebar_sitename = true,
        collapselevel = 1,
    ),
    pages = _PAGES
)

makedocs(
    sitename = "MadDiff",
    authors = "Sungho Shin",
    format = Documenter.LaTeX()
    pages = _PAGES
)

deploydocs(
    repo = "github.com/sshin23/MadDiff.jl.git"
)

