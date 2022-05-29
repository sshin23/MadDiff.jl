"""
    MadDiffMOI

`MadDiffMOI` is a submodule of `MadDiff`. `MadDifMOI` allows solving nonlinear optimization problems specified by `MathOptInterface` (https://github.com/jump-dev/JuMP.jl/tree/od/moi-nonlinear).
"""
module MadDiffMOI

import MadDiffCore: MadDiffCore, ifelse
using SpecialFunctions
using MathOptInterface: MathOptInterface, eval_objective

include("moi_interface.jl")

end # module
