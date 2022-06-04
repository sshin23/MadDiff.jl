"""
    MadDiffCore
Core algorithm for MadDiff.
"""
module MadDiffCore

import Base: RefValue, ==, >=, <=, >, <, +, -, -, ^, *, /, getindex, setindex!, add_sum, show, csch, tanh, csc, abs, exp10, acsc, acotd, exp2, acot, atan, cotd, sind, secd, atan, sinh, cot, cosh, acsch, atand, acosd, tand, cscd, asinh, cbrt, inv, acscd, asec, asind, cosd, abs2, sech, asech, sec, exp, coth

import NaNMath: NaNMath, sin, cos, tan, asin, acos, acosh, atanh, log, log2, log10, lgamma, log1p, pow, sqrt

"""
    AbstractExpression{T <: AbstractFloat}
Abstract type for expression, gradient, hessian, entry, and field evaluators.
"""
abstract type AbstractExpression{T <: AbstractFloat} end

"""
    Expression{T <: AbstractFloat}
Abstract type for expression evaluators.
"""
abstract type Expression{T <: AbstractFloat} <: AbstractExpression{T} end

"""
    Gradient{T <: AbstractFloat}
Abstract type for gradient evaluators.
"""
abstract type Gradient{T <: AbstractFloat} <: AbstractExpression{T} end

"""
    Hessian{T <: AbstractFloat}
Abstract type for hessian evaluators.
"""
abstract type Hessian{T <: AbstractFloat} <: AbstractExpression{T} end

"""
    Entry{T <: AbstractFloat}
Abstract type for entry evaluators.
"""
abstract type Entry{T <: AbstractFloat} <: AbstractExpression{T} end

"""
    Field{T <: AbstractFloat}
Abstract type for field evaluators.
"""
abstract type Field{T <: AbstractFloat} <: AbstractExpression{T} end

include("utils.jl")
include("expression.jl")
include("gradient.jl")
include("hessian.jl")
include("field.jl")
include("sourcesink.jl")
include("sparse.jl")
include("special.jl")
include("nlpcore.jl")
include("evaluator.jl")
include("register.jl")
include("functionlist.jl")
include("rules.jl")
include("print.jl")

export Source, Sink, Constant, Expression, Variable, Parameter, Field, Gradient, Jacobian, Hessian, SparseHessian, SparseJacobian, SparseGradient, Evaluator, FieldEvaluator, GradientEvaluator, SparseGradientEvaluator, HessianEvaluator, SparseHessianEvaluator, JacobianEvaluator, SparseJacobianEvaluator, function_evaluator, gradient_evaluator, sparse_gradient_evaluator, hessian_evaluator, sparse_hessian_evaluator, field_evaluator, jacobian_evaluator, sparse_jacobian_evaluator, obj, cons!, grad!, jac_coord!, hess_coord!, non_caching_eval, default_eval, ifelse, @register_univariate, @register_bivariate

end
