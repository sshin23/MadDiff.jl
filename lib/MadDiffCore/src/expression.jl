"""
    Constant{T <: AbstractFloat} <: Expression{T}
`Expression` for constants.

    Constant(x::T) where T <: AbstractFloat
Returns a `Constant` with value `x`.
# Example

```julia-repl
julia> e = Constant(1.)
1.0
julia> non_caching_eval(e, [1.,2.,3.])
1.0
```
"""
struct Constant{T <: AbstractFloat} <: Expression{T}
    x::T
end

"""
    Variable{T <: AbstractFloat} <: Expression{T}
`Expression` for variables.
"""
struct Variable{T <: AbstractFloat} <: Expression{T}
    index::Int
    ref::RefValue{T}
end

"""
    Variable{T}(n::Int) where T <: AbstractFloat
Returns a `Variable{T}` whose index is `n`.
"""
Variable{T}(n::Int) where T <: AbstractFloat = Variable{T}(n,RefValue{T}(0.))

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

"""
    Parameter{T <: AbstractFloat} <: Expression{T}
`Expression` for parameters.
"""
struct Parameter{T <: AbstractFloat} <: Expression{T}
    index::Int
    ref::RefValue{T}
end


"""
    Parameter{T}(n::Int) where T <: AbstractFloat
Returns a `Parameter{T}` whose index is `n`.
"""
Parameter{T}(n::Int) where T <: AbstractFloat= Parameter{T}(n,RefValue{T}(0.))

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
    Expression1{T <: AbstractFloat, F <: Function ,E <: Expression{T}}  <: Expression{T}
`Expression` for univariate function
"""
struct Expression1{T <: AbstractFloat, F <: Function ,E <: Expression{T}}  <: Expression{T}
    e1::E
    ref::RefValue{T}
    Expression1(f::F,e1::E) where {T <: AbstractFloat, F <: Function, E <: Expression{T}} =
        new{T,F,E}(e1,RefValue{T}(0.))
end

"""
    Expression2{T <: AbstractFloat, F <: Function,E1, E2} <: Expression{T}
`Expression` for bivariate function
"""
struct Expression2{T <: AbstractFloat, F <: Function,E1, E2} <: Expression{T}
    e1::E1
    e2::E2
    ref::RefValue{T}
    Expression2(f::F,e1::E1,e2::E2) where {T <: AbstractFloat, F, E1 <: Expression{T}, E2 <: Expression{T}} =
        new{T,F,E1,E2}(e1,e2,RefValue{T}(0.))
    Expression2(f::F,e1::E1,e2::E2) where {
        T <: AbstractFloat,
        F,
        E1,
        E2 <: Expression{T}
    } = new{T,F,E1,E2}(e1,e2,RefValue{T}(0.))
    Expression2(f::F,e1::E1,e2::E2) where {
        T <: AbstractFloat,
        F,
        E1 <: Expression{T},
        E2
    } = new{T,F,E1,E2}(e1,e2,RefValue{T}(0.))
end

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
    ExpressionIfElse(e0::E0,e1::E1,e2::E2) where {T, E0 <: Expression{T}, E1, E2} = new{T,E0,E1,E2}(e0,e1,e2,RefValue{T}(0.),RefValue{Bool}(true))
end

"""
    ExpressionSum{T <: AbstractFloat, E <: Expression{T}, I} <: Expression{T}
`Expression` for a summation of `Expression`s
"""
struct ExpressionSum{T <: AbstractFloat, E <: Expression{T}, I} <: Expression{T}
    inner::I
    es::Vector{E}
    ref::RefValue{T}
    ExpressionSum(es::Vector{E}) where {T <: AbstractFloat, E <: Expression{T}} = new{T,eltype(es),Nothing}(nothing,es,RefValue{T}(0.))
    ExpressionSum(inner::E,es) where {T <: AbstractFloat, E <: ExpressionSum{T}} = new{T,eltype(es),typeof(inner)}(inner,es,ref(inner))
end

