# abstract type AbstractConstant{T } <: Expression{T} end

# """
#     Constant{T} <: Expression{T}
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

struct ExpressionNull{T} <: Expression{T} end
ExpressionNull() = ExpressionNull{Float64}()

abstract type AbstractVariable{T } <: Expression{T}  end
"""
    Variable{T} <: Expression{T}
`Expression` for variables.
"""
struct Variable{T} <:  AbstractVariable{T}
    index::Int
    ref::RefValue{T}
end

"""
    Variable{T}(n::Int) where T
Returns a `Variable{T}` whose index is `n`.
"""
Variable{T}(n::Int) where T = Variable{T}(n,RefValue{T}())

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
Variable(n::Int) = Variable{Float64}(n)

abstract type AbstractParameter{T } <: Expression{T}  end
"""
    Parameter{T} <: Expression{T}
`Expression` for parameters.
"""
struct Parameter{T} <: AbstractParameter{T}
    index::Int
    ref::RefValue{T}
end


"""
    Parameter{T}(n::Int) where T
Returns a `Parameter{T}` whose index is `n`.
"""
Parameter{T}(n::Int) where T= Parameter{T}(n,RefValue{T}())

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
Parameter(n::Int) = Parameter{Float64}(n)

"""
    Expression1{T, F <: Function ,E <: Expression{T}}  <: Expression{T}
`Expression` for univariate function
"""
struct Expression1{T, F <: Function ,E <: Expression{T}}  <: Expression{T}
    e1::E
    ref::RefValue{T}
end
Expression1(f::F,e1::E) where {T, F <: Function, E <: Expression{T}} =  Expression1{T,F,E}(e1,RefValue{T}())

"""
    Expression2{T, F <: Function,E1, E2} <: Expression{T}
`Expression` for bivariate function
"""
struct Expression2{T, F <: Function,E1, E2} <: Expression{T}
    e1::E1
    e2::E2
    ref::RefValue{T}
end
Expression2(f::F,e1::E1, e2::E2) where {T, F, E1 <: Union{Expression{T},Real}, E2 <: Union{Expression{T},Real}} = Expression2{T,F,E1,E2}(e1,e2,RefValue{T}())


"""
    ExpressionIfElse{T,E0 <: Expression{T}, E1, E2} <: Expression{T}
`Expression` for `ifelse`
"""
struct ExpressionIfElse{T,E0 <: Expression{T}, E1, E2} <: Expression{T}
    e0::E0
    e1::E1
    e2::E2
    ref::RefValue{T}
    bref::RefValue{Bool}
end
ExpressionIfElse(e0::E0,e1::E1,e2::E2) where {T, E0 <: Expression{T}, E1, E2} = ExpressionIfElse{T,E0,E1,E2}(e0,e1,e2,RefValue{T}(),RefValue{Bool}(true))

"""
    ExpressionSum{T, E <: Expression{T}, I} <: Expression{T}
`Expression` for a summation of `Expression`s
"""
struct ExpressionSum{T, E <: Expression{T}, I} <: Expression{T}
    inner::I
    es::Vector{E}
    ref::RefValue{T}
end
ExpressionSum(es::Vector{E}) where {T, E <: Expression{T}} = ExpressionSum{T,eltype(es),Nothing}(nothing,es,RefValue{T}())
ExpressionSum(inner::E,es) where {T, E <: ExpressionSum{T}} = ExpressionSum{T,eltype(es),typeof(inner)}(inner,es,ref(inner))

