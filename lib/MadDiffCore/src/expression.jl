# abstract type AbstractConstant{T } <: Expression{T,RT} end

# """
#     Constant{T} <: Expression{T,RT}
# `Expression` for constants.

#     Constant(x::T) where T
# Returns a `Constant` with value `x`.
# # Example

# ```julia-repl
# julia> e = Constant(1.)
# 1.0
# julia> non_caching_eval(e, [1.,2.,3.])
# 1.0
# ```
# """
# struct Constant{T,R <: Real} <:  AbstractConstant{T}
#     x::R
# end

# """
#     Constant{T}(x::R) where {T, R <: Real}
# Returns a `Constant{T,R}` whose value is `x`.
# """
# Constant{T}(x::R) where {T, R <: Real} = Constant{T,R}(x)

struct ExpressionNull{T,RT} <: Expression{T,RT} end

abstract type AbstractVariable{T,RT}  <: Expression{T,RT}  end
"""
    Variable{T} <: Expression{T,RT}
`Expression` for variables.
"""
struct Variable{T,RT <: Ref{T}} <:  AbstractVariable{T,RT}
    index::Int
    ref::RT
end

"""
    Variable{T}(n::Int) where T
Returns a `Variable{T}` whose index is `n`.
"""
Variable{T,RT}(n::Int) where {T,RT} = Variable{T,RT}(n,RT())

"""
    Variable(n::Int) 
Returns a `Variable{Float64}` whose index is `n`
# Example

```julia-repl
julia> e = Variable(2)
x[2]
julia> non_caching_eval(e, [1.,2.,3.])
2.0
```
"""
Variable(n::Int) = Variable{Float64,RefValue{Float64}}(n)

abstract type AbstractParameter{T,RT} <: Expression{T,RT}  end
"""
    Parameter{T} <: Expression{T,RT}
`Expression` for parameters.
"""
struct Parameter{T,RT <: Ref{T}} <: AbstractParameter{T,RT}
    index::Int
    ref::RT
end


"""
    Parameter{T}(n::Int) where T
Returns a `Parameter{T}` whose index is `n`.
"""
Parameter{T,RT}(n::Int) where {T,RT} = Parameter{T,RT}(n,RT())

"""
    Parameter(n::Int) 
Returns a `Parameter{Float64}` whose index is `n`
# Example

```julia-repl
julia> e = Parameter(3)
p[3]
julia> non_caching_eval(e, [1.,2.,3.], [4.,5.,6.])
6.0
```
"""
Parameter(n::Int) = Parameter{Float64,RefValue{Float64}}(n)

"""
    Expression1{T, F <: Function ,E <: Expression{T,RT}}  <: Expression{T,RT}
`Expression` for univariate function
"""
struct Expression1{T, RT <: Ref{T}, F <: Function ,E <: Expression{T,RT}}  <: Expression{T,RT}
    e1::E
    ref::RT
end
Expression1(f::F,e1::E) where {T, RT, F <: Function, E <: Expression{T,RT}} =  Expression1{T,RT,F,E}(e1,RT())

"""
    Expression2{T, F <: Function,E1, E2} <: Expression{T,RT}
`Expression` for bivariate function
"""
struct Expression2{T, RT <: Ref{T}, F <: Function,E1, E2} <: Expression{T,RT}
    e1::E1
    e2::E2
    ref::RT
end
Expression2(f::F,e1::E1, e2::E2) where {T, RT, F, E1 <: Union{Expression{T,RT},Real}, E2 <: Union{Expression{T,RT},Real}} = Expression2{T,RT,F,E1,E2}(e1,e2,RT())


"""
    ExpressionIfElse{T,E0 <: Expression{T,RT}, E1, E2} <: Expression{T,RT}
`Expression` for `ifelse`
"""
struct ExpressionIfElse{T, RT <: Ref{T}, E0 <: Expression{T,RT}, E1, E2} <: Expression{T,RT}
    e0::E0
    e1::E1
    e2::E2
    ref::RT
    bref::RefValue{Bool}
end
ExpressionIfElse(e0::E0,e1::E1,e2::E2) where {T, RT, E0 <: Expression{T,RT}, E1, E2} = ExpressionIfElse{T,E0,E1,E2}(e0,e1,e2,RT(),RefValue{Bool}(true))

"""
    ExpressionSum{T, E <: Expression{T,RT}, I} <: Expression{T,RT}
`Expression` for a summation of `Expression`s
"""
struct ExpressionSum{T,RT <: Ref{T}, E <: Expression{T,RT}, I} <: Expression{T,RT}
    inner::I
    es::Vector{E}
    ref::RT
end
ExpressionSum(es::Vector{E}) where {T, RT, E <: Expression{T,RT}} = ExpressionSum{T,RT,eltype(es),Nothing}(nothing,es,RefValue{T}())
ExpressionSum(inner::E,es) where {T, RT, E <: ExpressionSum{T,RT}} = ExpressionSum{T,RT,eltype(es),typeof(inner)}(inner,es,ref(inner))

