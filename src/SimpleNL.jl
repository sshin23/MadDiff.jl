module SimpleNL

import Base: ==, >=, <=, getindex, setindex!, add_sum, copyto!, ndims, size, copyto!, convert, string, print, show, log2, max, csch, acoth, tanh, asecd, ^, *, cospi, csc, abs, hypot, /, rem2pi, -, acos, min, exp10, tan, acsc, log10, log1p, \, acotd, sin, -, exp2, acot, acosh, atan, mod, cotd, atanh, sind, rem, mod2pi, secd, atan, sinh, +, cot, cosh, acsch, sinpi, atand, cos, acosd, rad2deg, +, tand, cscd, asinh, cbrt, asin, transpose, inv, acscd, expm1, log, deg2rad, asec, asind, cosd, abs2, sech, sqrt, asech, sec, exp, coth
import SpecialFunctions: erfi, bessely, besselj, loggamma, erfcinv, hankelh2, hankelh1, erfcx, besselk, beta, invdigamma, bessely1, besselj1, dawson, airyaiprime, erf, digamma, gamma, airyai, airybi, erfinv, bessely0, erfc, trigamma, besseli, polygamma, logbeta, airybiprime, besselj0
import Requires: @require

export Constant, Variable, Parameter, Field, Gradient, Jacobian, Hessian, SparseHessian, SparseJacobian, SparseGradient, variable, parameter, objective, constraint,optimize!,instantiate!, num_variables, num_constraints, value, dual, objective_value, setvalue, set_lower_bound, set_upper_bound, lower_bound, upper_bound, set_optimizer, function_evaluator, gradient_evaluator, sparse_gradient_evaluator, hessian_evaluator, sparse_hessian_evaluator, field_evaluator, jacobian_evaluator, sparse_jacobian_evaluator, nlp_evaluator

abstract type Expression end
abstract type Gradient end
abstract type Hessian end
abstract type Field end
abstract type Entry end

mutable struct MyRef{T}
    x::T
end
@inline getindex(ref::MyRef) = ref.x
@inline setindex!(ref::MyRef,val) = ref.x = val

include("functionlist.jl")
include("expression.jl")
include("gradient.jl")
include("hessian.jl")
include("field.jl")
include("sourcesink.jl")
include("jacobian.jl")
include("lagrangianhessian.jl")
include("sparse.jl")
include("rules.jl")
include("template.jl")
include("utils.jl")
include("evaluator.jl")
include("model.jl")
include("print.jl")

function __init__()
    @require Ipopt="b6b21f68-93f8-5de0-b562-5493be1d77c9" @eval begin
        include("Interfaces/ipopt.jl")
    end
    @require MadNLP = "2621e9c9-9eb4-46b1-8089-e8c72242dfb6" @eval begin
        include("Interfaces/madnlp.jl")
    end
end

end
