struct Constant{T <: AbstractFloat} <: Expression{T}
    ref::MyRef{T}
    Constant(x::T) where T <: AbstractFloat  = new{T}(MyRef{T}(x))
end
struct Variable{T <: AbstractFloat} <: Expression{T}
    index::Int
    ref::MyRef{T}
    Variable{T}(n::Int) where T <: AbstractFloat = new{T}(n,MyRef{T}(0.))
end
struct Parameter{T <: AbstractFloat} <: Expression{T}
    index::Int
    ref::MyRef{T}
    Parameter{T}(n::Int) where T <: AbstractFloat= new{T}(n,MyRef{T}(0.))
end

struct Expression1{T <: AbstractFloat, F <: Function ,E} <: Expression{T}
    e1::E
    ref::MyRef{T}
    Expression1{T}(f::F,e1::E) where {T, F, E} = new{T,F,E}(e1,MyRef{T}(0.))
end
struct Expression2{T <: AbstractFloat, F <: Function,E1, E2} <: Expression{T}
    e1::E1
    e2::E2
    ref::MyRef{T}
    Expression2{T}(f::F,e1::E1,e2::E2) where {T, F, E1, E2} =
        new{T,F,E1,E2}(e1,e2,MyRef{T}(0.))
end
struct ExpressionSum{T <: AbstractFloat, E <: Expression{T}, I} <: Expression{T}
    inner::I
    es::Vector{E}
    ref::MyRef{T}
    ExpressionSum{T}(es::Vector{E}) where {T, E} =
        new{T,eltype(es),Nothing}(nothing,es,MyRef(0.))
    ExpressionSum{T}(inner::Vector{E},es) where {T, E} =
        new{T,eltype(es),typeof(inner)}(inner,es,ref(inner))
end

@inline function (e::ExpressionSum{T,E,I})(x,p=nothing) where {T,E,I}
    inner(e)(x,p)
    @simd for i in eachindex(e.es)
        @inbounds addrefval(e,e.es[i](x,p))
    end
    return refval(e)
end
@inline function (e::ExpressionSum{T,E,Nothing})(x,p=nothing) where {T,E}
    setrefval(e,.0)
    @simd for i in eachindex(e.es)
        @inbounds addrefval(e,e.es[i](x,p))
    end
    return refval(e)
end
@inline (e::Variable{T})(x,p=nothing) where {T <: Real} = setrefval(e,getindex(x,index(e)))
@inline (e::Parameter{T})(x,p=nothing) where {T <: AbstractFloat} = setrefval(e,getindex(p,index(e)))
@inline (e::Constant{T})(x,p=nothing) where {T <: AbstractFloat} = refval(e)

