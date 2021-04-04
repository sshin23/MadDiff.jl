struct GradientNull <: Gradient end
struct Gradient0 <: Gradient
    index::Int
    offset::Int
end
struct Gradient1{F,D1 <: Gradient} <: Gradient
    d1::D1
    fref1::MyRef{Float64}
    ref::MyRef{Float64}
    function Gradient1(e::Expression1{F,E1}, indexer = nothing) where {F,E1}
        d1 = Gradient(e.e1,indexer)
        return new{F,typeof(d1)}(d1,ref(e.e1),MyRef(1.))
    end
end
struct Gradient2F1{F,D1 <: Gradient, R<: Real} <: Gradient
    a::R
    d1::D1
    fref1::MyRef{Float64}
    ref::MyRef{Float64}
    function Gradient2F1(e::Expression2{F,E1,E2}, indexer = nothing) where {F,E1<:Real,E2}
        g1 = Gradient(e.e2,indexer)
        return new{F,typeof(g1),typeof(e.e1)}(e.e1,g1,ref(e.e2),MyRef(0.))
    end
end
struct Gradient2F2{F,D1 <: Gradient, R<: Real} <: Gradient
    a::R
    d1::D1
    fref1::MyRef{Float64}
    ref::MyRef{Float64}
    function Gradient2F2(e::Expression2{F,E1,E2}, indexer = nothing) where {F,E1,E2<:Real}
        g1 = Gradient(e.e1,indexer)
        return new{F,typeof(g1),typeof(e.e2)}(e.e2,g1,ref(e.e1),MyRef(0.))
    end
end
struct Gradient2{F,D1 <: Gradient,D2 <: Gradient} <: Gradient
    d1::D1
    d2::D2
    fref1::MyRef{Float64}
    fref2::MyRef{Float64}
    ref1::MyRef{Float64}
    ref2::MyRef{Float64}
    function Gradient2(e::Expression2{F,E1,E2}, indexer = nothing) where {F,E1,E2}
        g1 = Gradient(e.e1,indexer)
        g2 = Gradient(e.e2,indexer)
        return new{F,typeof(g1),typeof(g2)}(g1,g2,ref(e.e1),ref(e.e2),MyRef(0.),MyRef(0.))
    end
end
struct GradientSum{D,I} <: Gradient
    inner::I
    ds::Vector{D}
end

const GRADIENT_NULL = GradientNull()

@inline (::GradientNull)(z,x,p=nothing,d0=1) = nothing
@inline (d::Gradient0)(y,x,p=nothing,d0=1) = (@inbounds y[d.offset] += d0; nothing)
@inline (d::Gradient0)((j,y)::Tuple{Int,M},x,p=nothing,d0=1) where M <: AbstractMatrix = (@inbounds y[j,d.offset] += d0; nothing)
@inline (d::Gradient0)((j,y)::Tuple{Int,M},x,p=nothing,d0=1) where M <: AbstractVector = (@inbounds y[d.offset] += d0; nothing)

@inline function (d::GradientSum{D,I})(y,x,p=nothing,d0=1) where {D,I}
    inner(d)(y,x,p,d0)
    @simd for i in eachindex(d.ds)
        @inbounds d.ds[i](y,x,p,d0)
    end
end
@inline function (d::GradientSum{D,Nothing})(y,x,p=nothing,d0=1) where D
    @simd for i in eachindex(d.ds)
        @inbounds d.ds[i](y,x,p,d0)
    end
end

Gradient(e::V,(row,indexer)::Tuple{Int,Dict{Tuple{Int,Int},Int}}) where V <: Variable = Gradient0(index(e),set_indexer!(indexer,row,index(e)))
Gradient(e::V,indexer=nothing) where V <: Variable = Gradient0(index(e),set_indexer!(indexer,index(e)) )
Gradient(e::V,::Nothing) where V <: Variable = Gradient0(index(e),index(e))
Gradient(e::V,::Tuple{Int,Nothing}) where V <: Variable = Gradient0(index(e),index(e))
Gradient(e::P,::T) where {P<:Parameter,T} = GRADIENT_NULL
Gradient(e::C, indexer = nothing ) where {C<:Constant,T} = GRADIENT_NULL
Gradient(e::ExpressionSum{E,I},indexer = nothing) where {E,I} = GradientSum(Gradient(inner(e)),[Gradient(ee,indexer) for ee in e.es])
Gradient(e::ExpressionSum{E,Nothing},indexer = nothing) where {E,I} = GradientSum(nothing,[Gradient(ee,indexer) for ee in e.es])
Gradient(e::Expression1{F,E}, indexer = nothing) where {F,E} = Gradient1(e,indexer)
Gradient(e::Expression2{F,E1,E2}, indexer = nothing) where {F,E1,E2} = Gradient2(e,indexer)
Gradient(e::Expression2{F,E1,E2}, indexer = nothing) where {F,E1<:Real,E2} = Gradient2F1(e,indexer)
Gradient(e::Expression2{F,E1,E2}, indexer = nothing) where {F,E1,E2<:Real} = Gradient2F2(e,indexer)
