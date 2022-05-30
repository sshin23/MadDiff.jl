"""
    MadDiffMOI

`MadDiffMOI` is a submodule of `MadDiff`. `MadDifMOI` allows solving nonlinear optimization problems specified by `MathOptInterface` (https://github.com/jump-dev/JuMP.jl/tree/od/moi-nonlinear).
"""
module MadDiffMOI

import MadDiffCore: MadDiffCore, ifelse
using SpecialFunctions
using MathOptInterface: MathOptInterface, eval_objective

const MOI = MathOptInterface
const X = MadDiffCore.Variable()
const P = MadDiffCore.Parameter()
const _UNIDICT = [eval(f) for f in MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS]
const _MULTIDICT = [eval(f) for f in MOI.Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS]

include("moi_interface.jl")

export MadDiffAD

end # module
