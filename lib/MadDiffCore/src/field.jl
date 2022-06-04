struct FieldNull{T,RT} <: Field{T,RT} end
struct Field1{T,RT,E <: AbstractExpression{T,RT},I}  <: Field{T,RT}
    inner::I
    es::Vector{E}
    Field1(es::Vector{E}) where {T,RT, E <: AbstractExpression{T,RT}} = new{T,RT,E,Nothing}(nothing,es)
    Field1(inner::F,es::Vector{E}) where {T,RT,E <: AbstractExpression{T,RT},F} = new{T,RT,E,typeof(inner)}(inner,es)
    Field1(::Nothing,es::Vector{E}) where {T,RT,E <: AbstractExpression{T,RT}} = new{T,RT,E,Nothing}(nothing,es)
end
struct IndexedExpression{T <: AbstractFloat, RT <: Ref{T}, E <: Expression{T,RT}} <: Entry{T,RT}
    index::Int
    e::E
end

Field(f::FieldNull{T,RT},e::E) where {T, RT, E <: IndexedExpression} = Field1([e])
Field(f::F,e::E) where {F <: Field, E <: IndexedExpression} = _field(f,e) ? f : Field1(f,[e])

function _field(F::Field1{T,RT,E,I},e) where {T,RT,E,I} 
    if e isa eltype(F.es)
        push!(F.es,e)
        return true
    else
        return _field(inner(F),e)
    end
end
function _field(F::Field1{T,RT,E,Nothing},e) where {T,RT,E}
    if e isa eltype(F.es)
        push!(F.es,e)
        return true
    else
        return false
    end    
end

# Jacobian
struct JacobianEntry{T,RT,E <: Gradient{T,RT}} <: Entry{T,RT}
    index::Int
    e::E
end
Jacobian(f::Field1{T,RT,G,I},indexer = nothing) where {T,RT,G,I} = Field1(Jacobian(inner(f),indexer),[JacobianEntry(index(ie),Gradient(ie.e,(index(ie),indexer))) for ie in f.es])
Jacobian(f::Field1{T,RT,G,Nothing},indexer = nothing) where {T,RT,G} = Field1(nothing,[JacobianEntry(index(ie),Gradient(ie.e,(index(ie),indexer))) for ie in f.es])
Jacobian(::MadDiffCore.FieldNull{T,RT}) where {T,RT} = FieldNull{T,RT}()
