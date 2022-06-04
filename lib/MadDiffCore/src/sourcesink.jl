mutable struct Sink{I}
    inner::Union{Nothing,I}
end
struct Source{I <: Union{Variable,Parameter}} end
@inline (F::Sink{F1})(y,x,p=nothing) where F1 <: Field = inner(F)(y,x,p)

setindex!(s::Sink{Field{T,RT}},e::Expression{T,RT},i) where {T,RT} =
    (s.inner = Field(inner(s),IndexedExpression(i,e)))
getindex(::Source{Variable{T,RT}},n) where {T,RT} = Variable{T,RT}(n)
getindex(::Source{Parameter{T,RT}},n) where {T,RT} = Parameter{T,RT}(n)

Variable{T,RT}() where {T,RT} = Source{Variable{T,RT}}()
Parameter{T,RT}() where {T,RT} = Source{Parameter{T,RT}}()
Field{T,RT}() where {T,RT} = Sink{Field{T,RT}}(FieldNull{T,RT}())
Variable() = Variable{Float64,RefValue{Float64}}()
Parameter() = Parameter{Float64,RefValue{Float64}}()
Field() = Field{Float64,RefValue{Float64}}()
