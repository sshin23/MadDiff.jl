using Documenter, MadDiff

makedocs(
    sitename = "MadDiff.jl",
    modules = [MadDiff],
    authors = "Sungho Shin, Francois Pacaud, and contributors.",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        sidebar_sitename = true,
    ),
    pages = [
        "Home" => "index.md",
        "User Guide" => [
            "Getting Started"=>"guide.md",
            "Examples" => "examples.md",
        ],
        "How it works?" => "algorithms.md",
        "Citing MadDiff" => "citation.md",
        "Contributing" => "contrib.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/sshin23/MadDiff.jl.git"
)
