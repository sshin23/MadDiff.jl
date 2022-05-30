"""
    GradientNull{T <: AbstractFloat} <: Gradient{T}
`Gradient` of `Parameter` or `Constant`.
"""
struct GradientNull{T <: AbstractFloat} <: Gradient{T} end

"""
    Gradient0{T <: AbstractFloat} <: Gradient{T}
`Gradient` of `Variable`.
"""
struct Gradient0{T <: AbstractFloat} <: Gradient{T}
    index::Int
    offset::Int
end

"""
    Gradient1{T <: AbstractFloat, F, D1 <: Gradient} <: Gradient{T}
`Gradient` of `Expression1`.
"""
struct Gradient1{T <: AbstractFloat, F, D1 <: Gradient} <: Gradient{T}
    d1::D1
    fref1::RefValue{T}
    ref::RefValue{T}
    function Gradient1(e::Expression1{T,F,E1}, indexer = nothing) where {T,F,E1 <: Expression{T}}
        d1 = Gradient(e.e1,indexer)
        return new{T,F,typeof(d1)}(d1,ref(e.e1),RefValue{T}(1.))
    end
end

"""
    Gradient2F1{T <: AbstractFloat, F, D1 <: Gradient, R<: Real} <: Gradient{T}
`Gradient` of `Expression2` whose first argument is `<: Real`.
"""
struct Gradient2F1{T <: AbstractFloat, F, D1 <: Gradient, R<: Real} <: Gradient{T}
    a::R
    d1::D1
    fref1::RefValue{T}
    ref::RefValue{T}
    function Gradient2F1(e::Expression2{T,F,E1,E2}, indexer = nothing) where {T,F,E1<:Real,E2 <: Expression{T}}
        g1 = Gradient(e.e2,indexer)
        return new{T,F,typeof(g1),typeof(e.e1)}(e.e1,g1,ref(e.e2),RefValue{T}(0.))
    end
end

"""
    Gradient2F2{T <: AbstractFloat, F,D1 <: Gradient, R<: Real} <: Gradient{T}
`Gradient` of `Expression2` whose second argument is `<: Real`.
"""
struct Gradient2F2{T <: AbstractFloat, F,D1 <: Gradient, R<: Real} <: Gradient{T}
    a::R
    d1::D1
    fref1::RefValue{T}
    ref::RefValue{T}
    function Gradient2F2(e::Expression2{T,F,E1,E2}, indexer = nothing) where {T,F,E1 <: Expression{T},E2<:Real}
        g1 = Gradient(e.e1,indexer)
        return new{T,F,typeof(g1),typeof(e.e2)}(e.e2,g1,ref(e.e1),RefValue{T}(0.))
    end
end

"""
    Gradient2{T <: AbstractFloat, F,D1 <: Gradient,D2 <: Gradient} <: Gradient{T}
`Gradient` of `Expression2`.
"""
struct Gradient2{T <: AbstractFloat, F,D1 <: Gradient,D2 <: Gradient} <: Gradient{T}
    d1::D1
    d2::D2
    fref1::RefValue{T}
    fref2::RefValue{T}
    ref1::RefValue{T}
    ref2::RefValue{T}
    function Gradient2(e::Expression2{T,F,E1,E2}, indexer = nothing) where {T,F,E1 <: Expression{T},E2 <: Expression{T}}
        g1 = Gradient(e.e1,indexer)
        g2 = Gradient(e.e2,indexer)
        return new{T,F,typeof(g1),typeof(g2)}(g1,g2,ref(e.e1),ref(e.e2),RefValue{T}(0.),RefValue{T}(0.))
    end
    function Gradient2(e::Expression2{T,F,E1,E2}, indexer = nothing) where {T,F,E1 <: Expression{T},E2}
        g1 = Gradient(e.e1,indexer)
        g2 = Gradient(e.e2,indexer)
        return new{T,F,typeof(g1),typeof(g2)}(g1,g2,ref(e.e1),ref(e.e2),RefValue{T}(0.),RefValue{T}(0.))
    end
    function Gradient2(e::Expression2{T,F,E1,E2}, indexer = nothing) where {T,F,E1,E2 <: Expression{T}}
        g1 = Gradient(e.e1,indexer)
        g2 = Gradient(e.e2,indexer)
        return new{T,F,typeof(g1),typeof(g2)}(g1,g2,ref(e.e1),ref(e.e2),RefValue{T}(0.),RefValue{T}(0.))
    end
end

"""
    GradientSum{T <: AbstractFloat,D <: Gradient{T},I} <: Gradient{T}
`Gradient` of `ExpressionSum`.
"""
struct GradientSum{T <: AbstractFloat,D <: Gradient{T},I} <: Gradient{T}
    inner::I
    ds::Vector{D}
end

"""
    GradientIfElse{T, G1, G2} <: Gradient{T}
`Gradient` of `ExpressionIfElse`
"""
struct GradientIfElse{T, G1, G2} <: Gradient{T}
    d1::G1
    d2::G2
    bref::RefValue{Bool}
    function GradientIfElse(e::E, indexer = nothing) where {T, E <: ExpressionIfElse{T}}
        d1 = Gradient(e.e1,indexer)
        d2 = Gradient(e.e2,indexer)
        return new{T,typeof(d1),typeof(d2)}(
            d1,
            d2,
            e.bref
        )
    end
end

"""
    Gradient(e :: Expression{T}) where T
Returns the `Gradient` of an absraction `e`.
"""
Gradient(e::V,(row,indexer)::Tuple{Int,Dict{Tuple{Int,Int},Int}}) where {T, V <: Variable{T}} = Gradient0{T}(index(e),set_indexer!(indexer,row,index(e)))
Gradient(e::V,indexer=nothing) where {T, V <: Variable{T}} = Gradient0{T}(index(e),set_indexer!(indexer,index(e)) )
Gradient(e::V,::Nothing) where {T, V <: Variable{T}} = Gradient0{T}(index(e),index(e))
Gradient(e::V,::Tuple{Int,Nothing}) where {T, V <: Variable{T}} = Gradient0{T}(index(e),index(e))
Gradient(e::P,::N) where {T, P<:Parameter{T}, N} = GradientNull{T}()
Gradient(e::C, indexer = nothing ) where {T, C<:Constant{T}} = GradientNull{T}()
Gradient(e::ExpressionSum{T,E,I},indexer = nothing) where {T,E,I} = GradientSum(Gradient(inner(e),indexer),[Gradient(ee,indexer) for ee in e.es])
Gradient(e::ExpressionSum{T,E,Nothing},indexer = nothing) where {T,E,I} = GradientSum(nothing,[Gradient(ee,indexer) for ee in e.es])
Gradient(e::Expression1{T,F,E}, indexer = nothing) where {T,F,E} = Gradient1(e,indexer)
Gradient(e::Expression2{T,F,E1,E2}, indexer = nothing) where {T,F,E1,E2} = Gradient2(e,indexer)
Gradient(e::Expression2{T,F,E1,E2}, indexer = nothing) where {T,F,E1<:Real,E2} = Gradient2F1(e,indexer)
Gradient(e::Expression2{T,F,E1,E2}, indexer = nothing) where {T,F,E1,E2<:Real} = Gradient2F2(e,indexer)
Gradient(e::E, indexer = nothing) where {E <: ExpressionIfElse} = GradientIfElse(e)
