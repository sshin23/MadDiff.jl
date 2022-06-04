include("build.jl")

const _IS_COVERAGE = Base.JLOptions().code_coverage == 1
const _SUBMODULES = [
    "MadDiff", "MadDiffCore", "MadDiffSpecialFunctions", "MadDiffModels", "MadDiffMOI"
]

for submod in _SUBMODULES
    Pkg.test(submod; coverage = _IS_COVERAGE)
end
