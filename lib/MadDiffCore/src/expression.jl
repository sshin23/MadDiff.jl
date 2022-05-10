struct Constant <: Expression
    ref::MyRef{Float64}
    Constant(x::R) where R <: Real  = new(MyRef(x))
end
struct Variable <: Expression
    index::Int
    ref::MyRef{Float64}
    Variable(n::Int) = new(n,MyRef(0.))
end
struct Parameter <: Expression
    index::Int
    ref::MyRef{Float64}
    Parameter(n::Int) = new(n,MyRef(0.))
end

struct Expression1{F <: Function ,E <: Expression} <: Expression
    e1::E
    ref::MyRef{Float64}
    Expression1(f::F,e1::E) where {F,E} = new{F,E}(e1,MyRef(0.))
end
struct Expression2{F <: Function,E1,E2} <: Expression
    e1::E1
    e2::E2
    ref::MyRef{Float64}
    Expression2(f::F,e1::E1,e2::E2) where {F,E1,E2} = new{F,E1,E2}(e1,e2,MyRef(0.))
end
struct ExpressionSum{E <: Expression,I} <: Expression
    inner::I
    es::Vector{E}
    ref::MyRef{Float64}
    ExpressionSum(es) = new{eltype(es),Nothing}(nothing,es,MyRef(0.))
    ExpressionSum(inner,es) = new{eltype(es),typeof(inner)}(inner,es,ref(inner))
end

@inline function (e::ExpressionSum{E,I})(x,p=nothing) where {E,I}
    inner(e)(x,p)
    @simd for i in eachindex(e.es)
        @inbounds addrefval(e,e.es[i](x,p))
    end
    return refval(e)
end
@inline function (e::ExpressionSum{E,Nothing})(x,p=nothing) where E
    setrefval(e,.0)
    @simd for i in eachindex(e.es)
        @inbounds addrefval(e,e.es[i](x,p))
    end
    return refval(e)
end
@inline (e::Variable)(x,p=nothing)  = setrefval(e,getindex(x,index(e)))
@inline (e::Parameter)(x,p=nothing)  = setrefval(e,getindex(p,index(e)))
@inline (e::Constant)(x,p=nothing)  = refval(e)

