struct Dummy end
const DUMMY = Dummy()
getindex(::Dummy,key::Int) = 0.
getindex(::Tuple{Int,Dummy},key::Int) = 0.
setindex!(::Dummy,val,key) = nothing
setindex!(::Tuple{Int,Dummy},val,key) = nothing

index(e) = e.index
index1(e) = e.index1
index2(e) = e.index2
inner(e) = e.inner
ref(e) = e.ref
ref1(e) = e.ref1
ref2(e) = e.ref2
refval(e) = e.ref.x
frefval(e) = e.fref.x
frefval1(e) = e.fref1.x
frefval2(e) = e.fref2.x
refval1(e) = e.ref1.x
refval2(e) = e.ref2.x
refval11(e) = e.ref11.x
refval12(e) = e.ref12.x
refval21(e) = e.ref21.x
refval22(e) = e.ref22.x
setrefval(e,val) = ref(e)[] = val
setrefval1(e,val) = ref1(e)[] = val
setrefval2(e,val) = ref2(e)[] = val
addrefval(e,val) = ref(e)[] += val

# non-caching variants
@inline non_caching_eval(f::Sink{Field},y,x,p=nothing) = inner(f)(y,x,p)
@inline function non_caching_eval(f::Field1{E,I},y,x,p=nothing) where {E,I}
    non_caching_eval(inner(f),y,x,p)
    @simd for i in eachindex(f.es)
        non_caching_eval(f.es[i],y,x,p)
    end
end
@inline function non_caching_eval(f::Field1{E,Nothing},y,x,p=nothing) where E
    @simd for i in eachindex(f.es)
        non_caching_eval(f.es[i],y,x,p)
    end
end

@inline non_caching_eval(e::Variable{T},x,p=nothing) where {T <: AbstractFloat}  = @inbounds getindex(x,index(e))
@inline non_caching_eval(e::Parameter{T},x,p=nothing) where {T <: AbstractFloat} = @inbounds getindex(p,index(e))
@inline non_caching_eval(e::Constant{T},x,p=nothing) where {T <: AbstractFloat}  = refval(e)
@inline function non_caching_eval(e::ExpressionSum{T,E,I},x,p=nothing) where {T,E,I}
    res = non_caching_eval(e.inner,x,p)
    @simd for i in eachindex(e.es)
        @inbounds res = add_sum(res,non_caching_eval(e.es[i],x,p))
    end
    return res
end
@inline function non_caching_eval(e::ExpressionSum{T,E,Nothing},x,p=nothing) where {T,E}
    res = non_caching_eval(e.es[1],x,p)
    @simd for i in 2:length(e.es)
        @inbounds res = add_sum(res,non_caching_eval(e.es[i],x,p))
    end
    return res
end

@inline non_caching_eval(::GradientNull{T},z,x,p=nothing,d0=1) where T = nothing
@inline non_caching_eval(d::Gradient0{T},y,x,p=nothing,d0=1) where T = (@inbounds y[d.offset] += d0; nothing)
@inline non_caching_eval(d::Gradient0{T},(j,y)::Tuple{Int,M},x,p=nothing,d0=1) where {T, M <: AbstractMatrix} = (@inbounds y[j,d.offset] += d0; nothing)
@inline non_caching_eval(d::Gradient0{T},(j,y)::Tuple{Int,M},x,p=nothing,d0=1) where {T, M <: AbstractVector} = (@inbounds y[d.offset] += d0; nothing)
@inline function non_caching_eval(d::GradientSum{T,D,I},y,x,p=nothing,d0=1) where {T,D,I}
    non_caching_eval(inner(d),y,x,p,d0)
    @simd for i in eachindex(d.ds)
        non_caching_eval(d.ds[i],y,x,p,d0)
    end
end
@inline function non_caching_eval(d::GradientSum{T,D,Nothing},y,x,p=nothing,d0=1)  where {T,D}
    @simd for i in eachindex(d.ds)
        non_caching_eval(d.ds[i],y,x,p,d0)
    end
end
@inline non_caching_eval(e::IndexedExpression{E},y,x,p=nothing) where {E} = (y[index(e)] = non_caching_eval(e.e,x,p))
@inline non_caching_eval(e::JacobianEntry{E},y,x,p=nothing) where {E} = non_caching_eval(e.e,(index(e),y),x,p)
@inline non_caching_eval(e::FieldNull,y,x,p=nothing) = nothing
