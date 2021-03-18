module SimpleNLModels

import Base: string,show,print,getindex,setindex!,add_sum,+,-,*,^,/,==,<=,>=,one,zero,ImmutableDict
import DiffRules: diffrule
import JuMP: optimize!, value, dual, objective_value, setvalue, set_lower_bound, set_upper_bound, lower_bound, upper_bound, num_variables, num_constraints # to avoid conflict
import SpecialFunctions
import Requires: @require

const diffrules = [
    (:Base,:inv,1),(:Base,:abs2,1),
    (:Base,:exp,1),(:Base,:exp2,1),(:Base,:exp10,1),
    (:Base,:log,1),(:Base,:log2,1),(:Base,:log10,1),
    (:Base,:sin,1),(:Base,:cos,1),(:Base,:tan,1),(:Base,:csc,1),(:Base,:sec,1),(:Base,:cot,1),
    (:Base,:asin,1),(:Base,:acos,1),(:Base,:atan,1),(:Base,:acsc,1),(:Base,:asec,1),(:Base,:acot,1),
    (:Base,:sind,1),(:Base,:cosd,1),(:Base,:tand,1),(:Base,:cscd,1),(:Base,:secd,1),(:Base,:cotd,1),
    (:Base,:asind,1),(:Base,:acosd,1),(:Base,:atand,1),(:Base,:acscd,1),(:Base,:asecd,1),(:Base,:acotd,1),
    (:Base,:sinh,1),(:Base,:cosh,1),(:Base,:tanh,1),(:Base,:csch,1),(:Base,:sech,1),(:Base,:coth,1),
    (:Base,:asinh,1),(:Base,:acosh,1),(:Base,:atanh,1),(:Base,:acsch,1),(:Base,:asech,1),(:Base,:acoth,1),
    (:SpecialFunctions,:erfi,1),(:SpecialFunctions,:loggamma,1),(:SpecialFunctions,:erfcinv,1),
    (:SpecialFunctions,:erfcx,1),(:SpecialFunctions,:invdigamma,1),(:SpecialFunctions,:bessely1,1),
    (:SpecialFunctions,:besselj1,1),(:SpecialFunctions,:dawson,1),(:SpecialFunctions,:airyaiprime,1),
    (:SpecialFunctions,:erf,1),(:SpecialFunctions,:digamma,1),(:SpecialFunctions,:gamma,1),
    (:SpecialFunctions,:airyai,1),(:SpecialFunctions,:airybi,1),(:SpecialFunctions,:erfinv,1),
    (:SpecialFunctions,:bessely0,1),(:SpecialFunctions,:erfc,1),(:SpecialFunctions,:trigamma,1),
    (:SpecialFunctions,:airybiprime,1),(:SpecialFunctions,:besselj0,1),
    (:SpecialFunctions,:bessely,2),(:SpecialFunctions,:besselj,2),(:SpecialFunctions,:hankelh2,2),
    (:SpecialFunctions,:hankelh1,2),(:SpecialFunctions,:besselk,2),(:SpecialFunctions,:beta,2),
    (:SpecialFunctions,:besseli,2),(:SpecialFunctions,:polygamma,2),(:SpecialFunctions,:logbeta,2)
]

for (M,f,nargs) in diffrules
    @eval import $M.$f
end

const DEFAULT_VAR_STRING = "x"
const DEFAULT_PAR_STRING = "p"

export Source, Variable, Parameter, Term, func, deriv, 
    SimpleModel, variable, parameter, constraint, objective, instantiate!, optimize!,
    value, dual, objective_value, setvalue, set_lower_bound, set_upper_bound, lower_bound, upper_bound, set_optimizer,
    num_variables, num_constraints


include("expression.jl")
include("printexpression.jl")
include("nlp.jl")
include("print.jl")
include("Interfaces/common.jl")

function __init__()
    @require Ipopt="b6b21f68-93f8-5de0-b562-5493be1d77c9" @eval begin
        include("Interfaces/ipopt.jl")
    end
    @require MadNLP = "2621e9c9-9eb4-46b1-8089-e8c72242dfb6" @eval begin
        include("Interfaces/madnlp.jl")
    end
end

end # module
