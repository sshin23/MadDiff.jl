using Test, Pkg, MadDiff

const _IS_COVERAGE = Base.JLOptions().code_coverage == 1
const _SUBMODULES = [
    "MadDiffCore", "MadDiffModels", "MadDiffMOI"
]

for submod in _SUBMODULES
    @testset "$submod test" begin
        Pkg.test(submod; coverage = _IS_COVERAGE)
    end
end
