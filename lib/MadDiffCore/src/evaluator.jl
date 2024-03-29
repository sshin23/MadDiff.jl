for fname in [:default_eval]
    @eval begin
        @inline $fname(e::V,x,p=nothing) where {V <: AbstractVariable} = setrefval(e,getindex(x,index(e)))
        @inline $fname(e::P,x,p=nothing) where {P <: AbstractParameter} = setrefval(e,getindex(p,index(e)))
        @inline $fname(e::E, x, p = nothing) where {E <: ExpressionIfElse} =
            setrefval(e,setbrefval(e,non_caching_eval(e.e0,x,p)) ? $fname(e.e1,x,p) : $fname(e.e2,x,p))
    end
end


for fname in [:non_caching_eval]
    @eval begin
        @inline $fname(e::V,x,p=nothing) where {V <: AbstractVariable}  = @inbounds getindex(x,index(e))
        @inline $fname(e::P,x,p=nothing) where {P <: AbstractParameter} = @inbounds getindex(p,index(e))
        @inline $fname(laghess::LagrangianHessian,z,x,p,l,s) = $fname(laghess.f,z,x,p,l,s)
        @inline $fname(h::H,z,x,p,l,s) where H <: Hessian = $fname(h,z,x,p,s)
        @inline $fname(lagentry::LagrangianEntry{T,E},z,x,p,l) where {T,E} = $fname(lagentry.e,z,x,p,l[index(lagentry)])
        @inline $fname(e::E, x, p = nothing) where {E <: ExpressionIfElse} =
            non_caching_eval(e.e0,x,p) ? $fname(e.e1,x,p) : $fname(e.e2,x,p)


        @inline $fname(::HessianNull{T},z,x,p=nothing,h0=1) where T = nothing
        @inline function $fname(h::HessianD00{T},z,x,p=nothing,h0 = 1) where T
            islower(h) && @inbounds z[index1(h)::Int,index2(h)::Int] += h0
            return 
        end
        @inline function $fname(h::HessianD00S{T},z,x,p=nothing,h0 = 1) where T
            islower(h) && @inbounds z[index(h)::Int] += h0
            return 
        end
        @inline function $fname(h::HessianD10{T,H},z,x,p=nothing,h0 = 1) where {T,H}
            $fname(h.h,z,x,p,h0*refval(h))
            return
        end
        @inline function $fname(h::HessianD11{T,H},z,x,p=nothing,h0 = 1) where {T,H}
            $fname(h.h,z,x,p,h0*refval1(h)*refval2(h))
            return
        end
        @inline function $fname(h::HessianD21{T,H1,H2},z,x,p=nothing,h0=1) where {T,H1,H2}
            $fname(h.h1,z,x,p,h0*refval(h)*refval1(h))
            $fname(h.h2,z,x,p,h0*refval(h)*refval2(h))
            return
        end
        @inline function $fname(h::HessianD20{T,H1,H2},z,x,p=nothing,h0=1) where {T,H1,H2}
            $fname(h.h1,z,x,p,h0*refval1(h))
            $fname(h.h2,z,x,p,h0*refval2(h))
        end
        @inline function $fname(h::Hessian02{T,H11,H12,H21,H22},z,x,p=nothing,h0=1) where {T,H11,H12,H21,H22}
            $fname(h.h11,z,x,p,h0*refval11(h)*refval21(h))
            $fname(h.h12,z,x,p,h0*refval11(h)*refval22(h)) 
            $fname(h.h21,z,x,p,h0*refval12(h)*refval21(h)) 
            $fname(h.h22,z,x,p,h0*refval12(h)*refval22(h))
            return
        end
        @inline function $fname(h::Hessian11a{T,H},z,x,p=nothing,h0=1) where {T,H}
            $fname(h.h1,z,x,p,h0*refval(h))
            return
        end
        @inline function $fname(h::Hessian22a{T,H1,H2},z,x,p=nothing,h0=1) where {T,H1,H2}
            $fname(h.h1,z,x,p,h0*refval1(h))
            $fname(h.h2,z,x,p,h0*refval2(h))
            return
        end
        @inline function $fname(h::Hessian22m{T,H1,H2,H12},z,x,p=nothing,h0=1) where {T,H1,H2,H12}
            $fname(h.h1,z,x,p,h0*refval1(h))
            $fname(h.h2,z,x,p,h0*refval2(h))
            $fname(h.h12,z,x,p,h0)
            $fname(h.h21,z,x,p,h0)
            return
        end
        @inline $fname(h::H,z,x,p=nothing,h0=1.) where {H <: HessianIfElse } = brefval(h) ? $fname(h.h1,z,x,p,h0) : $fname(h.h2,z,x,p,h0)
    end
end


for fname in [:default_eval, :non_caching_eval]
    @eval begin
        @inline $fname(a::Real,x,p=nothing) = a
        @inline $fname(e::ExpressionNull{T},x,p=nothing) where T = 0.
        @inline $fname(::GradientNull{T},z,x,p=nothing,d0=1) where T = nothing
        @inline $fname(d::Gradient0{T},y,x,p=nothing,d0=1) where T = (@inbounds y[d.offset] += d0; nothing)
        @inline $fname(d::Gradient0{T},(j,y)::Tuple{Int,M},x,p=nothing,d0=1) where {T, M <: AbstractMatrix} = (@inbounds y[j,d.offset] += d0; nothing)
        @inline $fname(d::Gradient0{T},(j,y)::Tuple{Int,M},x,p=nothing,d0=1) where {T, M <: AbstractVector} = (@inbounds y[d.offset] += d0; nothing)
        @inline $fname(d::G,y,x,p=nothing,d0=1.) where {G <: GradientIfElse} = brefval(d) ? $fname(d.d1,y,x,p,d0) : $fname(d.d2,y,x,p,d0)
        
        @inline $fname(e::FieldNull{T},y,x,p=nothing) where T = nothing
        @inline $fname(e::JacobianEntry{T,E},y,x,p=nothing) where {T,E} = $fname(e.e,(index(e),y),x,p)
        @inline $fname(e::IndexedExpression{T,E},y,x,p=nothing) where {T,E} = (y[e.index] = $fname(e.e,x,p))
        @inline function $fname(d::GradientSum{T,D,I},y,x,p=nothing,d0=1) where {T,D,I}
            $fname(inner(d),y,x,p,d0)
            @simd for i in eachindex(d.ds)
                $fname(d.ds[i],y,x,p,d0)
            end
        end
        @inline function $fname(d::GradientSum{T,D,Nothing},y,x,p=nothing,d0=1)  where {T,D}
            @simd for i in eachindex(d.ds)
                $fname(d.ds[i],y,x,p,d0)
            end
        end
        @inline $fname(f::Sink{Field{T}},y,x,p=nothing) where T = inner(f)(y,x,p)
        @inline function $fname(f::Field1{T,E,I},y,x,p=nothing) where {T,E,I}
            $fname(inner(f),y,x,p)
            @simd for i in eachindex(f.es)
                $fname(f.es[i],y,x,p)
            end
        end
        @inline function $fname(f::Field1{T,E,Nothing},y,x,p=nothing) where {T,E}
            @simd for i in eachindex(f.es)
                $fname(f.es[i],y,x,p)
            end
        end
    end
end




@inline function non_caching_eval(f::Field1{T,E,I},z,x,p,l,s) where {T,E,I}
    non_caching_eval(inner(f),z,x,p,l,s)
    @simd for i in eachindex(f.es)
        non_caching_eval(f.es[i],z,x,p,l)
    end
end
@inline function non_caching_eval(h::HessianSum{T,I,H},z,x,p=nothing,h0=1) where {T,I,H} 
    non_caching_eval(inner(h),z,x,p,h0)
    @simd for i in eachindex(h.hs)
        @inbounds non_caching_eval(h.hs[i],z,x,p,h0)
    end
end
@inline function non_caching_eval(h::HessianSum{T,Nothing,H},z,x,p=nothing,h0=1) where {T,H}
    @simd for i in eachindex(h.hs)
        @inbounds non_caching_eval(h.hs[i],z,x,p,h0)
    end
end
@inline function default_eval(e::ExpressionSum{T,E,Nothing},x,p=nothing) where {T,E}
    setrefval(e,.0)
    @simd for i in eachindex(e.es)
        @inbounds addrefval(e,default_eval(e.es[i],x,p))
    end
    return refval(e)
end
@inline function non_caching_eval(e::ExpressionSum{T,E,Nothing},x,p=nothing) where {T,E}
    res = non_caching_eval(e.es[1],x,p)
    @simd for i in 2:length(e.es)
        @inbounds res = add_sum(res,non_caching_eval(e.es[i],x,p))
    end
    return res
end
@inline function default_eval(e::ExpressionSum{T,E,I},x,p=nothing) where {T,E,I}
    default_eval(inner(e),x,p)
    @simd for i in eachindex(e.es)
        @inbounds addrefval(e,default_eval(e.es[i],x,p))
    end
    return refval(e)
end
@inline function non_caching_eval(e::ExpressionSum{T,E,I},x,p=nothing) where {T,E,I}
    res = non_caching_eval(e.inner,x,p)
    @simd for i in eachindex(e.es)
        @inbounds res = add_sum(res,non_caching_eval(e.es[i],x,p))
    end
    return res
end

struct Evaluator{E}
    e::E
end
@inline (ev::Evaluator{E})(x,p=nothing) where E = non_caching_eval(ev.e,x,p)

struct FieldEvaluator{F}
    f::F
end
@inline (ev::FieldEvaluator{F})(y,x,p=nothing) where F = non_caching_eval(ev.f,y,x,p)

struct GradientEvaluator{E,D}
    e::E
    d::D
    function GradientEvaluator(e)
        d = Gradient(e)
        new{typeof(e),typeof(d)}(e,d)
    end
end
@inline function (ev::GradientEvaluator{E,D})(y,x,p=nothing) where {E, D}
    y.=0
    default_eval(ev.e,x,p)
    non_caching_eval(ev.d,y,x,p)
end

struct SparseGradientEvaluator{E,D}
    e::E
    d::D
    sparsity::Vector{Int}
    function SparseGradientEvaluator(e)
        d,sparsity = SparseGradient(e)
        new{typeof(e),typeof(d)}(e,d,sparsity)
    end
end
@inline function (ev::SparseGradientEvaluator{E,D})(y,x,p=nothing) where {E,D}
    y.=0
    default_eval(ev.e,x,p)
    non_caching_eval(ev.d,y,x,p)
end

struct HessianEvaluator{E,D,H}
    e::E
    d::D
    h::H
    function HessianEvaluator(e)
        d = Gradient(e)
        h = Hessian(e,d)
        new{typeof(e),typeof(d),typeof(h)}(e,d,h)
    end
end
@inline function (ev::HessianEvaluator{E,D,H})(z,x,p=nothing) where {E,D,H}
    z.=0
    default_eval(ev.e,x,p)
    default_eval(ev.d,DUMMY,x,p)
    non_caching_eval(ev.h,z,x,p)
end


struct SparseHessianEvaluator{E,D,H}
    e::E
    d::D
    h::H
    sparsity::Vector{Tuple{Int,Int}}
    function SparseHessianEvaluator(e)
        d = Gradient(e)
        h,sparsity = SparseHessian(e,d)
        new{typeof(e),typeof(d),typeof(h)}(e,d,h,sparsity)
    end
end
@inline function (ev::SparseHessianEvaluator{E,D,H})(z,x,p=nothing) where {E,D,H}
    z.=0
    default_eval(ev.e,x,p)
    default_eval(ev.d,DUMMY,x,p)
    non_caching_eval(ev.h,z,x,p)
end

struct JacobianEvaluator{F,J}
    f::F
    j::J
    function JacobianEvaluator(f)
        j = Jacobian(f)
        new{typeof(f),typeof(j)}(f,j)
    end
end
@inline function (ev::JacobianEvaluator{F,J})(y,x,p=nothing) where {F,J}
    y.=0
    default_eval(ev.f,DUMMY,x,p)
    non_caching_eval(ev.j,y,x,p)
end

struct SparseJacobianEvaluator{F,J}
    f::F
    j::J
    sparsity::Vector{Tuple{Int,Int}}
    function SparseJacobianEvaluator(f)
        j,sparsity = SparseJacobian(f)
        new{typeof(f),typeof(j)}(f,j,sparsity)
    end
end
@inline function (ev::SparseJacobianEvaluator{F,J})(y,x,p=nothing) where {F,J}
    y.=0
    default_eval(ev.f,DUMMY,x,p)
    non_caching_eval(ev.j,y,x,p)
end

FieldEvaluator(f::Sink) = FieldEvaluator(inner(f))
JacobianEvaluator(f::Sink) = JacobianEvaluator(inner(f))
SparseJacobianEvaluator(f::Sink) = SparseJacobianEvaluator(inner(f))
