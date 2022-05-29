module MadDiffCore

import NaNMath: sin, cos, tan, asin, acos, acosh, atanh, log, log2, log10, lgamma, log1p, pow, sqrt

import Base: zero, one, ==, >=, <=, >, <, getindex, setindex!, add_sum, copyto!, ndims, size, copyto!, convert, string, print, show, max, csch, acoth, tanh, asecd, ^, *, cospi, csc, abs, hypot, /, rem2pi, -, min, exp10, acsc, \, acotd, -, exp2, acot, atan, mod, cotd, sind, rem, mod2pi, secd, atan, sinh, +, cot, cosh, acsch, sinpi, atand, acosd, rad2deg, +, tand, cscd, asinh, cbrt, transpose, inv, acscd, expm1, deg2rad, asec, asind, cosd, abs2, sech, asech, sec, exp, coth, RefValue

import SpecialFunctions: SpecialFunctions, erfi, bessely, besselj, loggamma, erfcinv, hankelh2, hankelh1, erfcx, besselk, beta, invdigamma, bessely1, besselj1, dawson, airyaiprime, erf, digamma, gamma, airyai, airybi, erfinv, bessely0, erfc, trigamma, besseli, polygamma, logbeta, airybiprime, besselj0

export Source, Sink, Constant, Expression, Variable, Parameter, Field, Gradient, Jacobian, Hessian, SparseHessian, SparseJacobian, SparseGradient, function_evaluator, gradient_evaluator, sparse_gradient_evaluator, hessian_evaluator, sparse_hessian_evaluator, field_evaluator, jacobian_evaluator, sparse_jacobian_evaluator, obj, cons!, grad!, jac_coord!, hess_coord!

abstract type AbstractExpression{T <: AbstractFloat} end 
abstract type Expression{T <: AbstractFloat} <: AbstractExpression{T} end
abstract type Gradient{T <: AbstractFloat} <: AbstractExpression{T} end
abstract type Hessian{T <: AbstractFloat} <: AbstractExpression{T} end
abstract type Field{T <: AbstractFloat} <: AbstractExpression{T} end
abstract type Entry{T <: AbstractFloat} <: AbstractExpression{T} end

include("functionlist.jl")
include("expression.jl")
include("gradient.jl")
include("hessian.jl")
include("field.jl")
include("sourcesink.jl")
include("nlpcore.jl")
include("sparse.jl")
include("rules.jl")
include("utils.jl")
include("evaluator.jl")
include("print.jl")

end
