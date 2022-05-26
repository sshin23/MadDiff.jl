for fname in [:default_eval, :default_eval_threaded]
    @eval begin
        @inline $fname(e::Variable{T},x,p=nothing) where {T <: Real} = setrefval(e,getindex(x,index(e)))
        @inline $fname(e::Parameter{T},x,p=nothing) where {T <: AbstractFloat} = setrefval(e,getindex(p,index(e)))
    end


    for (f,df,ddf) in f_nargs_1
        @eval begin
            @inline $fname(e::Expression1{T,typeof($f),E},x,p=nothing) where {T,E} = setrefval(e,$f($fname(e.e1,x,p)))
            
            @inline function $fname(d::Gradient1{T,typeof($f),F1},y,x,p=nothing,d0=1.) where {T,F1}
                setrefval(d,$df(frefval1(d)))
                $fname(d.d1,y,x,p,d0*refval(d))
                return 
            end
            
        end
    end

    for (f,df1,df2,ddf11,ddf12,ddf22) in f_nargs_2
        @eval begin
            @inline $fname(e::Expression2{T,typeof($f),F1,F2},x,p=nothing) where {T,F1,F2} = setrefval(e,$f($fname(e.e1,x,p),$fname(e.e2,x,p)))
            @inline $fname(e::Expression2{T,typeof($f),F1,F2},x,p=nothing) where {T,F1<:Real,F2} = setrefval(e,$f(e.e1,$fname(e.e2,x,p)))
            @inline $fname(e::Expression2{T,typeof($f),F1,F2},x,p=nothing) where {T,F1,F2<:Real} = setrefval(e,$f($fname(e.e1,x,p),e.e2))

            @inline function $fname(d::Gradient2F1{T,typeof($f),F1,R},y,x,p=nothing,d0=1.) where {T,F1,R}
                setrefval(d,$df2(d.a,frefval1(d)))
                $fname(d.d1,y,x,p,d0*refval(d))
                return 
            end
            @inline function $fname(d::Gradient2F2{T,typeof($f),F1,R},y,x,p=nothing,d0=1.) where {T,F1,R}
                setrefval(d,$df1(frefval1(d),d.a))
                $fname(d.d1,y,x,p,d0*refval(d))
                return 
            end
            @inline function $fname(d::Gradient2{T,typeof($f),D1,D2},y,x,p=nothing,d0 = 1) where {T,D1,D2}
                setrefval1(d,$df1(frefval1(d),frefval2(d)))
                setrefval2(d,$df2(frefval1(d),frefval2(d)))
                $fname(d.d1,y,x,p,d0*refval1(d))
                $fname(d.d2,y,x,p,d0*refval2(d))
                return
            end
        end
    end
end


for fname in [:non_caching_eval, :non_caching_eval_threaded]
    @eval begin
        @inline $fname(e::Variable{T},x,p=nothing) where {T <: AbstractFloat}  = @inbounds getindex(x,index(e))
        @inline $fname(e::Parameter{T},x,p=nothing) where {T <: AbstractFloat} = @inbounds getindex(p,index(e))
        @inline $fname(laghess::LagrangianHessian,z,x,p,l,s) = $fname(laghess.f,z,x,p,l,s)
        @inline $fname(h::H,z,x,p,l,s) where H <: Hessian = $fname(h,z,x,p,s)
        @inline $fname(lagentry::LagrangianEntry{T,E},z,x,p,l) where {T,E} = $fname(lagentry.e,z,x,p,l[index(lagentry)])



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
    end

    for (f,df,ddf) in f_nargs_1
        @eval begin
            @inline $fname(e::Expression1{T,typeof($f),E},x,p=nothing) where {T,E} = $f($fname(e.e1,x,p))
            @inline function $fname(d::Gradient1{T,typeof($f),D1},y,x,p=nothing,d0=1.)  where {T,D1}
                $fname(d.d1,y,x,p,d0*$df(frefval1(d)))
                return 
            end
            @inline function $fname(h::Hessian11{T,typeof($f),H1,H11},z,x,p=nothing,h0=1) where {T,H1,H11}
                $fname(h.h1,z,x,p,h0*refval(h))
                $fname(h.h11,z,x,p,h0*$ddf(frefval(h)))
                return
            end

        end
    end


    for (f,df1,df2,ddf11,ddf12,ddf22) in f_nargs_2
        @eval begin
            @inline $fname(e::Expression2{T,typeof($f),F1,F2},x,p=nothing) where {T,F1,F2} = $f($fname(e.e1,x,p),$fname(e.e2,x,p))
            @inline $fname(e::Expression2{T,typeof($f),F1,F2},x,p=nothing) where {T,F1<:Real,F2} = $f(e.e1,$fname(e.e2,x,p))
            @inline $fname(e::Expression2{T,typeof($f),F1,F2},x,p=nothing) where {T,F1,F2<:Real} = $f($fname(e.e1,x,p),e.e2)
            
            @inline function $fname(d::Gradient2F1{T,typeof($f),D1,R},y,x,p=nothing,d0=1.) where {T,D1,R}
                $fname(d.d1,y,x,p,d0*$df2(d.a,frefval1(d)))
                return 
            end
            @inline function $fname(d::Gradient2F2{T,typeof($f),D1,R},y,x,p=nothing,d0=1.) where {T,D1,R}
                $fname(d.d1,y,x,p,d0*$df1(frefval1(d),d.a))
                return 
            end
            @inline function $fname(d::Gradient2{T,typeof($f),D1,D2},y,x,p=nothing,d0 = 1) where {T,D1,D2}
                $fname(d.d1,y,x,p,d0*$df1(frefval1(d),frefval2(d)))
                $fname(d.d2,y,x,p,d0*$df2(frefval1(d),frefval2(d)))
                return
            end
            @inline function $fname(h::Hessian11F1{T,typeof($f),H1,H11,R},z,x,p=nothing,h0=1) where {T,H1,H11,R}
                $fname(h.h1,z,x,p,h0*refval(h))
                $fname(h.h11,z,x,p,h0*$ddf22(h.a,frefval(h)))
                return
            end
            @inline function $fname(h::Hessian11F2{T,typeof($f),H1,H11,R},z,x,p=nothing,h0=1) where {T,H1,H11,R}
                $fname(h.h1,z,x,p,h0*refval(h))
                $fname(h.h11,z,x,p,h0*$ddf11(frefval(h),h.a))
                return
            end
            @inline function $fname(h::Hessian22{T,typeof($f),H1,H2,H11,H12,H22},z,x,p=nothing,h0=1) where {T,H1,H2,H11,H12,H22}
                ddf12 = $ddf12(frefval1(h),frefval2(h))
                $fname(h.h1,z,x,p,h0*refval1(h))
                $fname(h.h2,z,x,p,h0*refval2(h))
                $fname(h.h11,z,x,p,h0*$ddf11(frefval1(h),frefval2(h)))
                $fname(h.h12,z,x,p,h0*ddf12)
                $fname(h.h21,z,x,p,h0*ddf12)
                $fname(h.h22,z,x,p,h0*$ddf22(frefval1(h),frefval2(h)))
                return
            end

        end
    end

end


for fname in [:default_eval, :default_eval_threaded, :non_caching_eval, :non_caching_eval_threaded]
    @eval begin 
        @inline $fname(e::Constant{T},x,p=nothing) where {T <: AbstractFloat}  = refval(e)
        @inline $fname(::GradientNull{T},z,x,p=nothing,d0=1) where T = nothing
        @inline $fname(d::Gradient0{T},y,x,p=nothing,d0=1) where T = (@inbounds y[d.offset] += d0; nothing)
        @inline $fname(d::Gradient0{T},(j,y)::Tuple{Int,M},x,p=nothing,d0=1) where {T, M <: AbstractMatrix} = (@inbounds y[j,d.offset] += d0; nothing)
        @inline $fname(d::Gradient0{T},(j,y)::Tuple{Int,M},x,p=nothing,d0=1) where {T, M <: AbstractVector} = (@inbounds y[d.offset] += d0; nothing)
        @inline $fname(e::FieldNull{T},y,x,p=nothing) where T = nothing
        @inline $fname(e::JacobianEntry{T,E},y,x,p=nothing) where {T,E} = $fname(e.e,(index(e),y),x,p)
        @inline $fname(e::IndexedExpression{T,E},y,x,p=nothing) where {T,E} = (y[e.index] = $fname(e.e,x,p))
    end
end

for (fname, fname_threaded) in [(:default_eval, :default_eval_threaded),
                                (:non_caching_eval, :non_caching_eval_threaded)]
    @eval begin
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

        @inline function $fname_threaded(d::GradientSum{T,D,I},y,x,p=nothing,d0=1) where {T,D,I}
            $fname_threaded(inner(d),y,x,p,d0)
            Threads.@threads for i in eachindex(d.ds)
                $fname_threaded(d.ds[i],y,x,p,d0)
            end
        end
        @inline function $fname_threaded(d::GradientSum{T,D,Nothing},y,x,p=nothing,d0=1)  where {T,D}
            Threads.@threads for i in eachindex(d.ds)
                $fname_threaded(d.ds[i],y,x,p,d0)
            end
        end
        @inline function $fname_threaded(f::Field1{T,E,I},y,x,p=nothing) where {T,E,I}
            $fname_threaded(inner(f),y,x,p)
            Threads.@threads for i in eachindex(f.es)
                $fname_threaded(f.es[i],y,x,p)
            end
        end
        @inline function $fname_threaded(f::Field1{T,E,Nothing},y,x,p=nothing) where {T,E}
            Threads.@threads for i in eachindex(f.es)
                $fname_threaded(f.es[i],y,x,p)
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
@inline function non_caching_eval_threaded(f::Field1{T,E,I},z,x,p,l,s) where {T,E,I}
    non_caching_eval_threaded(inner(f),z,x,p,l,s)
    Threads.@threads for i in eachindex(f.es)
        non_caching_eval_threaded(f.es[i],z,x,p,l)
    end
end
@inline function non_caching_eval_threaded(h::HessianSum{T,I,H},z,x,p=nothing,h0=1) where {T,I,H} 
    non_caching_eval_threaded(inner(h),z,x,p,h0)
    Threads.@threads for i in eachindex(h.hs)
        @inbounds non_caching_eval_threaded(h.hs[i],z,x,p,h0)
    end
end
@inline function non_caching_eval_threaded(h::HessianSum{T,Nothing,H},z,x,p=nothing,h0=1) where {T,H}
    Threads.@threads for i in eachindex(h.hs)
        @inbounds non_caching_eval_threaded(h.hs[i],z,x,p,h0)
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
@inline function default_eval_threaded(e::ExpressionSum{T,E,Nothing},x,p=nothing) where {T,E}
    setrefval(e,.0)
    Threads.@threads for i in eachindex(e.es)
        @inbounds addrefval(e,default_eval_threaded(e.es[i],x,p))
    end
    return refval(e)
end
@inline function non_caching_eval_threaded(e::ExpressionSum{T,E,Nothing},x,p=nothing) where {T,E}
    res = non_caching_eval_threaded(e.es[1],x,p)
    Threads.@threads for i in 2:length(e.es)
        @inbounds res = add_sum(res,non_caching_eval_threaded(e.es[i],x,p))
    end
    return res
end
@inline function default_eval_threaded(e::ExpressionSum{T,E,I},x,p=nothing) where {T,E,I}
    default_eval_threaded(inner(e),x,p)
    Threads.@threads for i in eachindex(e.es)
        @inbounds addrefval(e,default_eval_threaded(e.es[i],x,p))
    end
    return refval(e)
end
@inline function non_caching_eval_threaded(e::ExpressionSum{T,E,I},x,p=nothing) where {T,E,I}
    res = non_caching_eval_threaded(e.inner,x,p)
    Threads.@threads for i in eachindex(e.es)
        @inbounds res = add_sum(res,non_caching_eval_threaded(e.es[i],x,p))
    end
    return res
end

