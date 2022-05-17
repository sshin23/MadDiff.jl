mutable struct Sink{T}
    inner::Union{Nothing,T}
end
struct Source{T <: Union{Variable,Parameter}} end
@inline (F::Sink{F1})(y,x,p=nothing) where F1 <: Field = inner(F)(y,x,p)

setindex!(s::Sink{Field{T}},e::Expression{T},i) where T =
    (s.inner = Field(inner(s),IndexedExpression(i,e)))
getindex(::Source{Variable{T}},n) where T = Variable{T}(n)
getindex(::Source{Parameter{T}},n) where T = Parameter{T}(n)

Variable{T}() where T <: AbstractFloat = Source{Variable{T}}()
Parameter{T}() where T <: AbstractFloat = Source{Parameter{T}}()
Field{T}() where T = Sink{Field{T}}(FieldNull{T}())
Variable() = Variable{Float64}()
Parameter() = Parameter{Float64}()
Field() = Field{Float64}()
