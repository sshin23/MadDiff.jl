"""
    GradientNull{T,RT} <: Gradient{T,RT}
`Gradient` of `Parameter` or `Real`.
"""
struct GradientNull{T,RT} <: Gradient{T,RT} end

"""
    Gradient0{T,RT} <: Gradient{T,RT}
`Gradient` of `Variable`.
"""
struct Gradient0{T,RT} <: Gradient{T,RT}
    index::Int
    offset::Int
end

"""
    Gradient1{T, F, D1 <: Gradient} <: Gradient{T,RT}
`Gradient` of `Expression1`.
"""
struct Gradient1{T, RT <: Ref{T}, F, D1 <: Gradient} <: Gradient{T,RT}
    d1::D1
    fref1::RT
    ref::RT
    function Gradient1(e::Expression1{T,RT,F,E1}, indexer = nothing) where {T,RT,F,E1}
        d1 = Gradient(e.e1,indexer)
        return new{T,RT,F,typeof(d1)}(d1,ref(e.e1),RT())
    end
end

"""
    Gradient2F1{T, F, D1 <: Gradient, R<: Real} <: Gradient{T,RT}
`Gradient` of `Expression2` whose first argument is `<: Real`.
"""
struct Gradient2F1{T, RT <: Ref{T}, F, D1 <: Gradient, R<: Real} <: Gradient{T,RT}
    a::R
    d1::D1
    fref1::RT
    ref::RT
    function Gradient2F1(e::Expression2{T,RT,F,E1,E2}, indexer = nothing) where {T,RT,F,E1<:Real,E2 <: Expression{T,RT}}
        g1 = Gradient(e.e2,indexer)
        return new{T,RT,F,typeof(g1),typeof(e.e1)}(e.e1,g1,ref(e.e2),RT())
    end
end

"""
    Gradient2F2{T, F,D1 <: Gradient, R<: Real} <: Gradient{T,RT}
`Gradient` of `Expression2` whose second argument is `<: Real`.
"""
struct Gradient2F2{T, RT <: Ref{T}, F,D1 <: Gradient, R<: Real} <: Gradient{T,RT}
    a::R
    d1::D1
    fref1::RT
    ref::RT
    function Gradient2F2(e::Expression2{T,RT,F,E1,E2}, indexer = nothing) where {T,RT,F,E1 <: Expression{T,RT}, E2<:Real}
        g1 = Gradient(e.e1,indexer)
        return new{T,RT,F,typeof(g1),typeof(e.e2)}(e.e2,g1,ref(e.e1),RT())
    end
end

"""
    Gradient2{T, F,D1 <: Gradient,D2 <: Gradient} <: Gradient{T,RT}
`Gradient` of `Expression2`.
"""
struct Gradient2{T, RT <: Ref{T}, F,D1 <: Gradient,D2 <: Gradient} <: Gradient{T,RT}
    d1::D1
    d2::D2
    fref1::RT
    fref2::RT
    ref1::RT
    ref2::RT
end
for (T1,T2) in [(:Expression,:Expression),(:Expression,:Real),(:Real,:Expression)]
    @eval function Gradient2(e::Expression2{T,RT,F,E1,E2}, indexer = nothing) where {T,RT,F,E1 <: $T1,E2 <: $T2}
        g1 = Gradient(e.e1,indexer)
        g2 = Gradient(e.e2,indexer)
        return Gradient2{T,RT,F,typeof(g1),typeof(g2)}(g1,g2,ref(e.e1),ref(e.e2),RT(),RT())
    end
end

"""
    GradientSum{T,D <: Gradient{T,RT},I} <: Gradient{T,RT}
`Gradient` of `ExpressionSum`.
"""
struct GradientSum{T, RT, D <: Gradient{T,RT}, I} <: Gradient{T,RT}
    inner::I
    ds::Vector{D}
end

"""
    GradientIfElse{T, G1, G2} <: Gradient{T,RT}
`Gradient` of `ExpressionIfElse`
"""
struct GradientIfElse{T, RT <: Ref{T}, G1, G2} <: Gradient{T,RT}
    d1::G1
    d2::G2
    bref::RefValue{Bool}
    function GradientIfElse(e::E, indexer = nothing) where {T, RT, E <: ExpressionIfElse{T,RT}}
        d1 = Gradient(e.e1,indexer)
        d2 = Gradient(e.e2,indexer)
        return new{T,RT,typeof(d1),typeof(d2)}(
            d1,
            d2,
            e.bref
        )
    end
end

"""
    Gradient(e :: Expression{T,RT}) where T
Returns the `Gradient` of an absraction `e`.
"""
Gradient(e::Real,indexer=nothing) = GradientNull{typeof(float(e)),RefValue{typeof(float(e))}}()
Gradient(e::V,(row,indexer)::Tuple{Int,Dict{Tuple{Int,Int},Int}}) where {T, RT, V <: AbstractVariable{T,RT}} = Gradient0{T,RT}(index(e),set_indexer!(indexer,row,index(e)))
Gradient(e::V,indexer=nothing) where {T, RT, V <: AbstractVariable{T,RT}} = Gradient0{T,RT}(index(e),set_indexer!(indexer,index(e)) )
Gradient(e::V,::Nothing) where {T, RT, V <: AbstractVariable{T,RT}} = Gradient0{T,RT}(index(e),index(e))
Gradient(e::V,::Tuple{Int,Nothing}) where {T, RT, V <: AbstractVariable{T,RT}} = Gradient0{T,RT}(index(e),index(e))
Gradient(e::P,::N) where {T, RT, P<:AbstractParameter{T,RT}, N} = GradientNull{T,RT}()
Gradient(e::ExpressionNull{T,RT}, indexer = nothing ) where {T,RT} = GradientNull{T,RT}()
Gradient(e::ExpressionSum{T,RT, E,I},indexer = nothing) where {T,RT,E,I} = GradientSum(Gradient(inner(e),indexer),[Gradient(ee,indexer) for ee in e.es])
Gradient(e::ExpressionSum{T,RT, E,Nothing},indexer = nothing) where {T,RT,E,I} = GradientSum(nothing,[Gradient(ee,indexer) for ee in e.es])
Gradient(e::Expression1{T,RT,F,E}, indexer = nothing) where {T,RT,F,E} = Gradient1(e,indexer)
Gradient(e::Expression2{T,RT,F,E1,E2}, indexer = nothing) where {T,RT,F,E1,E2} = Gradient2(e,indexer)
Gradient(e::Expression2{T,RT,F,E1,E2}, indexer = nothing) where {T,RT,F,E1<:Real,E2} = Gradient2F1(e,indexer)
Gradient(e::Expression2{T,RT,F,E1,E2}, indexer = nothing) where {T,RT,F,E1,E2<:Real} = Gradient2F2(e,indexer)
Gradient(e::E, indexer = nothing) where {E <: ExpressionIfElse} = GradientIfElse(e)
