module NonlinearModels

export Expression, PrintSymbol

import Base: getindex,+,-,*,^,/,inv,add_sum,
    exp,sin,cos,log,tan,asin,acos,atan,sinh,cosh,tanh,
    string,show,print

import MadNLP: NonlinearProgram, INITIAL
import Ipopt: IpoptProblem, createProblem

include("expression.jl")
include("Interfaces/common.jl")
include("Interfaces/madnlp.jl")
include("Interfaces/ipopt.jl")

end # module
