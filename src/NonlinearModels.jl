module NonlinearModels

export Expression, PrintSymbol

import Base: string,show,print, getindex, add_sum ,+,-,*,^,/,inv, one
import DiffRules: diffrule, diffrules

import MadNLP: NonlinearProgram, INITIAL
import Ipopt: IpoptProblem, createProblem

include("expression.jl")
include("nlp.jl")
include("Interfaces/common.jl")
include("Interfaces/madnlp.jl")
include("Interfaces/ipopt.jl")

end # module
