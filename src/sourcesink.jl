mutable struct Sink{T} inner::Union{Nothing,T} end
struct Source{T <: Union{Variable,Parameter}} end
@inline (F::Sink{Field})(y,x,p=nothing) = inner(F)(y,x,p)

setindex!(s::Sink{Field},e,i) = (s.inner = Field(inner(s),IndexedExpression(i,e)))
getindex(::Source{Variable},n) = Variable(n)
getindex(::Source{Parameter},n) = Parameter(n)

Variable() = Source{Variable}()
Parameter() = Source{Parameter}()
