mutable struct Sink{T}
    inner::Union{Nothing,T}
end
struct Source{T <: Union{Variable,Parameter}} end
@inline (F::Sink{Field})(y,x,p=nothing) = inner(F)(y,x,p)

setindex!(s::Sink{Field},e,i) = (s.inner = Field(inner(s),IndexedExpression{Float64}(i,e)))
getindex(::Source{Variable{T}},n) where {T <: Real} = Variable{T}(n)
getindex(::Source{Parameter},n) = Parameter(n)

Variable{T}() where T <: AbstractFloat = Source{Variable{T}}()
Variable{T}() where T <: AbstractFloat = Source{Variable{T}}()
Variable() = Variable{Float64}()
Parameter() = Source{Parameter{Float64}}()
