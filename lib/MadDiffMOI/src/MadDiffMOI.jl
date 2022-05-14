module MadDiffMOI

using MadDiffCore: MadDiffCore, NLPCore, Expression, Constant
import MadDiffModels: MadDiffModels, MadDiffModel, variable, constraint
using SpecialFunctions
using MathOptInterface

export MadDiffAD

include("moi_interface.jl")

end # module
