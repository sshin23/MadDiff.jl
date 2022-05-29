using Test, Pkg, MadDiff

const _SUBMODULES = [
    "MadDiffCore", "MadDiffModels", "MadDiffMOI"
]

for submod in _SUBMODULES
    @testset "$submod test" begin
        Pkg.test(submod)
    end
end
