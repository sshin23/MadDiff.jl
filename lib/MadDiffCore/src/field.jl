struct FieldNull{T} <: Field{T} end
struct Field1{T,E <: Term{T},I}  <: Field{T}
    inner::I
    es::Vector{E}
    Field1(es::Vector{E}) where {T, E <: Term{T}} = new{T,E,Nothing}(nothing,es)
    Field1(inner::F,es::Vector{E}) where {T,E <: Term{T},F} = new{T,E,typeof(inner)}(inner,es)
    Field1(::Nothing,es::Vector{E}) where {T,E <: Term{T}} = new{T,E,Nothing}(nothing,es)
end
struct IndexedExpression{T <: AbstractFloat, E <: Expression{T}} <: Entry{T}
    index::Int
    e::E
end

@inline (F::FieldNull{T})(y,x,p=nothing) where T = nothing
@inline function (F::Field1{T,E,I})(y,x,p=nothing) where {T,E,I}
    inner(F)(y,x,p)
    @simd for i in eachindex(F.es)
        F.es[i](y,x,p)
    end
end
@inline function (F::Field1{T,E,Nothing})(y,x,p=nothing) where {T,E}
    @simd for i in eachindex(F.es)
        F.es[i](y,x,p)
    end
end
@inline (e::IndexedExpression{T,E})(y,x,p=nothing) where {T,E} = (y[e.index] = e.e(x,p))

# Field(e::E) where E <: IndexedExpression= Field1([e])
Field(f::FieldNull{T},e::E) where {T, E <: IndexedExpression} = Field1([e])
Field(f::F,e::E) where {F <: Field, E <: IndexedExpression} = _field(f,e) ? f : Field1(f,[e])

function _field(F::Field1{T,E,I},e) where {T,E,I} 
    if e isa eltype(F.es)
        push!(F.es,e)
        return true
    else
        return _field(inner(F),e)
    end
end
function _field(F::Field1{T,E,Nothing},e) where {T,E}
    if e isa eltype(F.es)
        push!(F.es,e)
        return true
    else
        return false
    end    
end
