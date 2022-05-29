using Documenter, MadDiff

makedocs(
    sitename = "MadDiff",
    modules = [MadDiff],
    authors = "Sungho Shin, Francois Pacaud, and contributors.",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        sidebar_sitename = true,
        collapselevel = 1,
    ),
    pages = [
        "Home" => "index.md",
        "Quick Start"=>"guide.md",
        "How it works?" => "algorithms.md",
        "API Manual" => [
            "MadDiffCore" => "core.md",
            "MadDiffModels" => "models.md",
            "MadDiffMOI" => "moi.md",
        ]
    ],
)

deploydocs(
    repo = "github.com/sshin23/MadDiff.jl.git"
)
