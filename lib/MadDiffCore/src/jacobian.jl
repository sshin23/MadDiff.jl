struct JacobianEntry{E} <: Entry
    index::Int
    e::E
end
# struct Jacobian{G <: Gradient, I}
#     inner::I
#     ds::Vector{Tuple{Int,G}}
# end
# @inline function (J::Jacobian{G,I})(y,x,p=nothing) where {G,I}
#     inner(J)(y,x)
#     @simd for i in eachindex(J.ds)
#         (j,d) = J.ds[i]
#         d((j,y),x,p)
#     end
# end
# @inline function (J::Jacobian{G,Nothing})(y,x,p=nothing) where G
#     @simd for i in eachindex(J.ds)
#         (j,d) = J.ds[i]
#         d((j,y),x,p)
#     end
# end
@inline (e::JacobianEntry{E})(y,x,p=nothing) where {E} = e.e((index(e),y),x,p)
# Jacobian(s::Sink{Field}) = Jacobian(inner(s))
Jacobian(f::Field1{G,I},indexer = nothing) where {G,I} = Field1(Jacobian(inner(f),indexer),[JacobianEntry(index(ie),Gradient(ie.e,(index(ie),indexer))) for ie in f.es])
Jacobian(f::Field1{G,Nothing},indexer = nothing) where G = Field1(nothing,[JacobianEntry(index(ie),Gradient(ie.e,(index(ie),indexer))) for ie in f.es])


