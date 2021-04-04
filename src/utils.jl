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
    res = e.inner(x,p)
    @simd for i in eachindex(e.es)
        @inbounds res += e.es[i](x,p)
    end
    return res
end
@inline function non_caching_eval(e::ExpressionSum{E,Nothing},x,p=nothing) where E
    res = 0
    @simd for i in eachindex(e.es)
        @inbounds res += e.es[i](x,p)
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


get_terms(e::E) where E <: Expression = (terms=Expression[];_get_terms!(e,terms);terms)
_get_terms!(e::Expression2{typeof(+),E1,E2},terms) where {E1,E2} = (_get_terms!(e.e1,terms);_get_terms!(e.e2,terms))
function _get_terms!(e::ExpressionSum{E,I},terms) where {E,I}
    _get_terms!(e.inner,terms)
    for e in e.es
        _get_terms!(e,terms)
    end
end
function _get_terms!(e::ExpressionSum{E,Nothing},terms) where E
    for e in e.es
        _get_terms!(e,terms)
    end
end
_get_terms!(e::E,terms) where E = push!(terms,e)

get_entries(f::Sink{Field}) = get_entries(inner(f))
get_entries(f::Field1{E,I}) where {E,I} = (entries=[];_get_entries!(f,entries);entries)
function _get_entries!(f::Field1{E,I},entries) where {E,I}
    _get_entries!(inner(f),entries)
    for e in f.es
        _get_entries!(e,entries)
    end
end
function _get_entries!(f::Field1{E,Nothing},entries) where E
    for e in f.es
        _get_entries!(e,entries)
    end
end
_get_entries!(e::E,entries) where E = push!(entries,e)

get_entries_expr(f::Sink{Field}) = get_entries_expr(inner(f))
function get_entries_expr(f::Field1{E,I}) where {E,I}
    entries = get_entries(f)
    fs = Vector{Expression}(undef,length(entries))
    for c in entries
        fs[index(c)] = c.e
    end
    return fs
end

sparsity(e::E) where E = (indices=Int[];_sparsity!(e::E,indices);indices)
_sparsity!(e::Variable,indices) = @inbounds union!(indices,index(e))
_sparsity!(e::Parameter,indices) = nothing
_sparsity!(e::Constant,indices) = nothing
_sparsity!(e::Expression1{F,F1},indices) where {F,F1} = _sparsity!(e.e1,indices)
_sparsity!(e::Expression2{F,F1,F2},indices) where {F,F1,F2} = (_sparsity!(e.e1,indices);_sparsity!(e.e2,indices))
_sparsity!(e::Expression2{F,F1,F2},indices) where {F,F1<:Real,F2} = _sparsity!(e.e2,indices)
_sparsity!(e::Expression2{F,F1,F2},indices) where {F,F1,F2<:Real} = _sparsity!(e.e1,indices)
function _sparsity!(e::ExpressionSum{E,I},indices) where {E,I}
    _sparsity!(e.inner,indices)
    @simd for i in eachindex(e.es)
        @inbounds  _sparsity!(e.es[i],indices)
    end
end
function _sparsity!(e::ExpressionSum{E,Nothing},indices) where E
    @simd for i in eachindex(e.es)
        @inbounds _sparsity!(e.es[i],indices)
    end
    return res
end



@inline non_caching_grad(e::Variable,y,x,p=nothing,d0=1) = (@inbounds y[e.index] += d0; nothing)
