struct FieldNull <: Field end
struct Field1{E,I}  <: Field
    inner::I
    es::Vector{E}
    Field1(es::Vector{E}) where E = new{E,Nothing}(nothing,es)
    Field1(inner::F,es::Vector{E}) where {E,F} = new{E,typeof(inner)}(inner,es)
    Field1(::Nothing,es::Vector{E}) where E = new{E,Nothing}(nothing,es)
end
struct IndexedExpression{T <: AbstractFloat, E <: Expression{T}} <: Entry{T}
    index::Int
    e::E
end

FIELD_NULL = FieldNull()

@inline (F::FieldNull)(y,x,p=nothing) = nothing
@inline function (F::Field1{E,I})(y,x,p=nothing) where {E,I}
    inner(F)(y,x,p)
    @simd for i in eachindex(F.es)
        F.es[i](y,x,p)
    end
end
@inline function (F::Field1{E,Nothing})(y,x,p=nothing) where E
    @simd for i in eachindex(F.es)
        F.es[i](y,x,p)
    end
end
@inline (e::IndexedExpression{T,E})(y,x,p=nothing) where {T,E} = (y[e.index] = e.e(x,p))

Field() = Sink{Field}(FIELD_NULL)
# Field(e::E) where E <: IndexedExpression= Field1([e])
Field(f::FieldNull,e::E) where E <: IndexedExpression = Field1([e])
Field(f::F,e::E) where {F <: Field, E <: IndexedExpression} = _field(f,e) ? f : Field1(f,[e])


function _field(F::Field1{E,I},e) where {E,I} 
    if e isa eltype(F.es)
        push!(F.es,e)
        return true
    else
        return _field(inner(F),e)
    end
end
function _field(F::Field1{E,Nothing},e) where E
    if e isa eltype(F.es)
        push!(F.es,e)
        return true
    else
        return false
    end    
end
