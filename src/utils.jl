index(e) = e.index
index1(e) = e.index1
index2(e) = e.index2
inner(e) = e.inner
ref(e) = e.ref
ref1(e) = e.ref1
ref2(e) = e.ref2
# ref11(e) = e.ref11
# ref12(e) = e.ref12
# ref21(e) = e.ref21
# ref22(e) = e.ref22
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
# addrefval1(e,val) = ref1(e)[] += val
# addrefval2(e,val) = ref2(e)[] += val

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

@inline non_caching_eval(e::Variable,x,p=nothing)  = @inbounds getindex(x,index(e))
@inline non_caching_eval(e::Parameter,x,p=nothing)  = @inbounds getindex(p,index(e))
@inline non_caching_eval(e::Constant,x,p=nothing)  = refval(e)
@inline function non_caching_eval(e::ExpressionSum{E,I},x,p=nothing) where {E,I}
    res = non_caching_eval(e.inner,x,p)
    @simd for i in eachindex(e.es)
        @inbounds res = add_sum(res,non_caching_eval(e.es[i],x,p))
    end
    return res
end
@inline function non_caching_eval(e::ExpressionSum{E,Nothing},x,p=nothing) where E
    res = non_caching_eval(e.es[1],x,p)
    @simd for i in 2:length(e.es)
        @inbounds res = add_sum(res,non_caching_eval(e.es[i],x,p))
    end
    return res
end

@inline non_caching_eval(::GradientNull,z,x,p=nothing,d0=1) = nothing
@inline non_caching_eval(d::Gradient0,y,x,p=nothing,d0=1) = (@inbounds y[d.offset] += d0; nothing)
@inline non_caching_eval(d::Gradient0,(j,y)::Tuple{Int,M},x,p=nothing,d0=1) where M <: AbstractMatrix = (@inbounds y[j,d.offset] += d0; nothing)
@inline non_caching_eval(d::Gradient0,(j,y)::Tuple{Int,M},x,p=nothing,d0=1) where M <: AbstractVector = (@inbounds y[d.offset] += d0; nothing)
@inline function non_caching_eval(d::GradientSum{D,I} where {D,I},y,x,p=nothing,d0=1)
    non_caching_eval(inner(d),y,x,p,d0)
    @simd for i in eachindex(d.ds)
        non_caching_eval(d.ds[i],y,x,p,d0)
    end
end
@inline function non_caching_eval(d::GradientSum{D,Nothing} where D,y,x,p=nothing,d0=1)
    @simd for i in eachindex(d.ds)
        non_caching_eval(d.ds[i],y,x,p,d0)
    end
end
@inline non_caching_eval(e::IndexedExpression{E},y,x,p=nothing) where {E} = (y[index(e)] = e.e(x,p))
@inline non_caching_eval(e::JacobianEntry{E},y,x,p=nothing) where {E} = e.e((index(e),y),x,p)
@inline non_caching_eval(e::FieldNull,y,x,p=nothing) = nothing
