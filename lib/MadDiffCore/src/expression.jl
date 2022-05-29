struct Constant{T <: AbstractFloat} <: Expression{T}
    ref::RefValue{T}
    Constant(x::T) where T <: AbstractFloat  = new{T}(RefValue{T}(x))
end
struct Variable{T <: AbstractFloat} <: Expression{T}
    index::Int
    ref::RefValue{T}
    Variable{T}(n::Int) where T <: AbstractFloat = new{T}(n,RefValue{T}(0.))
end
struct Parameter{T <: AbstractFloat} <: Expression{T}
    index::Int
    ref::RefValue{T}
    Parameter{T}(n::Int) where T <: AbstractFloat= new{T}(n,RefValue{T}(0.))
end

struct Expression1{T <: AbstractFloat, F <: Function ,E} <: Expression{T}
    e1::E
    ref::RefValue{T}
    Expression1(f::F,e1::E) where {T, F, E <: Expression{T}} = new{T,F,E}(e1,RefValue{T}(0.))
end
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
struct ExpressionIfElse{T,E0 <: Expression{T}, E1, E2} <: Expression{T}
    e0::E0
    e1::E1
    e2::E2
    ref::RefValue{Bool}
    ExpressionIfElse(e0::E0,e1::E1,e2::E2) where {T, E0 <: Expression{T}, E1, E2} = new{T,E0,E1,E2}(e0,e1,e2,RefValue{Bool}(true))
end

struct ExpressionSum{T <: AbstractFloat, E <: Expression{T}, I} <: Expression{T}
    inner::I
    es::Vector{E}
    ref::RefValue{T}
    ExpressionSum(es::Vector{E}) where {T <: AbstractFloat, E <: Expression{T}} = new{T,eltype(es),Nothing}(nothing,es,RefValue{T}(0.))
    ExpressionSum(inner::E,es) where {T <: AbstractFloat, E <: ExpressionSum{T}} = new{T,eltype(es),typeof(inner)}(inner,es,ref(inner))
end

