include("build.jl")

Pkg.test([
    "MadDiffCore", "MadDiffModels", "MadDiffMOI", "MadDiff"
]; coverage = true)
