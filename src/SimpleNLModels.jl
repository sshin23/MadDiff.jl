module SimpleNLModels

import Base: string,show,print,getindex,add_sum,+,-,*,^,/,==,<=,>=
import DiffRules: diffrule
import MadNLP, Ipopt
import JuMP: optimize!, value, dual, getobjectivevalue, setvalue, set_lower_bound, set_upper_bound # to avoid conflict

const diffrules = [
    (:Base,:inv,1),(:Base,:abs2,1),
    (:Base,:exp,1),(:Base,:exp2,1),(:Base,:exp10,1),
    (:Base,:log,1),(:Base,:log2,1),(:Base,:log10,1),
    (:Base,:sin,1),(:Base,:cos,1),(:Base,:tan,1),(:Base,:csc,1),(:Base,:sec,1),(:Base,:cot,1),
    (:Base,:asin,1),(:Base,:acos,1),(:Base,:atan,1),(:Base,:acsc,1),(:Base,:asec,1),(:Base,:acot,1),
    (:Base,:sind,1),(:Base,:cosd,1),(:Base,:tand,1),(:Base,:cscd,1),(:Base,:secd,1),(:Base,:cotd,1),
    (:Base,:asind,1),(:Base,:acosd,1),(:Base,:atand,1),(:Base,:acscd,1),(:Base,:asecd,1),(:Base,:acotd,1),
    (:Base,:sinh,1),(:Base,:cosh,1),(:Base,:tanh,1),(:Base,:csch,1),(:Base,:sech,1),(:Base,:coth,1),
    (:Base,:asinh,1),(:Base,:acosh,1),(:Base,:atanh,1),(:Base,:acsch,1),(:Base,:asech,1),(:Base,:acoth,1)
]

for (M,f,nargs) in diffrules
    @eval import $M.$f
end

const DEFAULT_VAR_STRING = "x"
const DEFAULT_PAR_STRING = "p"

export Source, Variable, Parameter, Term, func, deriv,
    variable, parameter, constraint, objective, instantiate!, optimize!,
    value, dual, getobjectivevalue, setvalue, set_lower_bound, set_upper_bound

const Reals = [Int,Float64]

include("print.jl")
include("expression.jl")
include("nlp.jl")
include("Interfaces/common.jl")
include("Interfaces/madnlp.jl")
include("Interfaces/ipopt.jl")

end # module
