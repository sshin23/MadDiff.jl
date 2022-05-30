include("build.jl")

Pkg.test("MadDiff"; coverage = true)
