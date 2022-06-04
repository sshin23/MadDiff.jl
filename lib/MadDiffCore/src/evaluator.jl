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
        @inline $fname(lagentry::LagrangianEntry{T,RT,E},z,x,p,l) where {T,RT,E} = $fname(lagentry.e,z,x,p,l[index(lagentry)])
        @inline $fname(e::E, x, p = nothing) where {E <: ExpressionIfElse} =
            non_caching_eval(e.e0,x,p) ? $fname(e.e1,x,p) : $fname(e.e2,x,p)


        @inline $fname(::HessianNull{T,RT},z,x,p=nothing,h0=1) where {T,RT} = nothing
        @inline function $fname(h::H,z,x,p=nothing,h0 = 1) where H <: HessianD00
            islower(h) && @inbounds z[index1(h)::Int,index2(h)::Int] += h0
            return 
        end
        @inline function $fname(h::H,z,x,p=nothing,h0 = 1) where H <: HessianD00S
            islower(h) && @inbounds z[index(h)::Int] += h0
            return 
        end
        @inline function $fname(h::H,z,x,p=nothing,h0 = 1) where H <: HessianD10
            $fname(h.h,z,x,p,h0*refval(h))
            return
        end
        @inline function $fname(h::H,z,x,p=nothing,h0 = 1) where H <: HessianD11
            $fname(h.h,z,x,p,h0*refval1(h)*refval2(h))
            return
        end
        @inline function $fname(h::H,z,x,p=nothing,h0=1) where H <: HessianD21
            $fname(h.h1,z,x,p,h0*refval(h)*refval1(h))
            $fname(h.h2,z,x,p,h0*refval(h)*refval2(h))
            return
        end
        @inline function $fname(h::H,z,x,p=nothing,h0=1) where H <: HessianD20
            $fname(h.h1,z,x,p,h0*refval1(h))
            $fname(h.h2,z,x,p,h0*refval2(h))
        end
        @inline function $fname(h::H,z,x,p=nothing,h0=1) where H <: Hessian02
            $fname(h.h11,z,x,p,h0*refval11(h)*refval21(h))
            $fname(h.h12,z,x,p,h0*refval11(h)*refval22(h)) 
            $fname(h.h21,z,x,p,h0*refval12(h)*refval21(h)) 
            $fname(h.h22,z,x,p,h0*refval12(h)*refval22(h))
            return
        end
        @inline function $fname(h::H,z,x,p=nothing,h0=1) where H <: Hessian11a
            $fname(h.h1,z,x,p,h0*refval(h))
            return
        end
        @inline function $fname(h::H,z,x,p=nothing,h0=1) where H <: Hessian22a
            $fname(h.h1,z,x,p,h0*refval1(h))
            $fname(h.h2,z,x,p,h0*refval2(h))
            return
        end
        @inline function $fname(h::H,z,x,p=nothing,h0=1) where H <: Hessian22m
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
        @inline $fname(e::ExpressionNull{T,RT},x,p=nothing) where {T, RT} = 0.
        @inline $fname(::GradientNull{T,RT},z,x,p=nothing,d0=1) where {T, RT} = nothing
        @inline $fname(d::Gradient0{T,RT},y,x,p=nothing,d0=1) where {T, RT} = (@inbounds y[d.offset] += d0; nothing)
        @inline $fname(d::Gradient0{T,RT},(j,y)::Tuple{Int,M},x,p=nothing,d0=1) where {T, RT, M <: AbstractMatrix} = (@inbounds y[j,d.offset] += d0; nothing)
        @inline $fname(d::Gradient0{T,RT},(j,y)::Tuple{Int,M},x,p=nothing,d0=1) where {T, RT, M <: AbstractVector} = (@inbounds y[d.offset] += d0; nothing)
        @inline $fname(d::G,y,x,p=nothing,d0=1.) where {G <: GradientIfElse} = brefval(d) ? $fname(d.d1,y,x,p,d0) : $fname(d.d2,y,x,p,d0)
        
        @inline $fname(e::FieldNull{T,RT},y,x,p=nothing) where {T,RT} = nothing
        @inline $fname(e::JacobianEntry{T,RT,E},y,x,p=nothing) where {T,RT,E} = $fname(e.e,(index(e),y),x,p)
        @inline $fname(e::IndexedExpression{T,RT,E},y,x,p=nothing) where {T,RT,E} = (y[e.index] = $fname(e.e,x,p))
        @inline function $fname(d::GradientSum{T,RT,D,I},y,x,p=nothing,d0=1) where {T,RT,D,I}
            $fname(inner(d),y,x,p,d0)
            @simd for i in eachindex(d.ds)
                $fname(d.ds[i],y,x,p,d0)
            end
        end
        @inline function $fname(d::GradientSum{T,RT,D,Nothing},y,x,p=nothing,d0=1)  where {T,RT,D}
            @simd for i in eachindex(d.ds)
                $fname(d.ds[i],y,x,p,d0)
            end
        end
        @inline $fname(f::Sink{Field{T,RT}},y,x,p=nothing) where {T,RT} = inner(f)(y,x,p)
        @inline function $fname(f::Field1{T,RT,E,I},y,x,p=nothing) where {T,RT,E,I}
            $fname(inner(f),y,x,p)
            @simd for i in eachindex(f.es)
                $fname(f.es[i],y,x,p)
            end
        end
        @inline function $fname(f::Field1{T,RT,E,Nothing},y,x,p=nothing) where {T,RT,E}
            @simd for i in eachindex(f.es)
                $fname(f.es[i],y,x,p)
            end
        end
    end
end




@inline function non_caching_eval(f::Field1{T,RT,E,I},z,x,p,l,s) where {T,RT,E,I}
    non_caching_eval(inner(f),z,x,p,l,s)
    @simd for i in eachindex(f.es)
        non_caching_eval(f.es[i],z,x,p,l)
    end
end
@inline function non_caching_eval(h::HessianSum{T,RT,I,H},z,x,p=nothing,h0=1) where {T,RT,I,H} 
    non_caching_eval(inner(h),z,x,p,h0)
    @simd for i in eachindex(h.hs)
        @inbounds non_caching_eval(h.hs[i],z,x,p,h0)
    end
end
@inline function non_caching_eval(h::HessianSum{T,RT,Nothing,H},z,x,p=nothing,h0=1) where {T,RT,H}
    @simd for i in eachindex(h.hs)
        @inbounds non_caching_eval(h.hs[i],z,x,p,h0)
    end
end
@inline function default_eval(e::ExpressionSum{T,RT,E,Nothing},x,p=nothing) where {T,RT,E}
    setrefval(e,.0)
    @simd for i in eachindex(e.es)
        @inbounds addrefval(e,default_eval(e.es[i],x,p))
    end
    return refval(e)
end
@inline function non_caching_eval(e::ExpressionSum{T,RT,E,Nothing},x,p=nothing) where {T,RT,E}
    res = non_caching_eval(e.es[1],x,p)
    @simd for i in 2:length(e.es)
        @inbounds res = add_sum(res,non_caching_eval(e.es[i],x,p))
    end
    return res
end
@inline function default_eval(e::ExpressionSum{T,RT,E,I},x,p=nothing) where {T,RT,E,I}
    default_eval(inner(e),x,p)
    @simd for i in eachindex(e.es)
        @inbounds addrefval(e,default_eval(e.es[i],x,p))
    end
    return refval(e)
end
@inline function non_caching_eval(e::ExpressionSum{T,RT,E,I},x,p=nothing) where {T,RT,E,I}
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
