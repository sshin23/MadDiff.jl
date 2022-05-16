# constructors definitions
# export non_caching_grad

for (f,df,ddf) in f_nargs_1
    @eval begin
        $f(e::E) where {E <: Expression} = Expression1($f,e)
        @inline (e::Expression1{T,typeof($f),E})(x,p=nothing) where {T,E} = setrefval(e,$f(e.e1(x,p)))
        
        @inline function (d::Gradient1{typeof($f),F1})(y,x,p=nothing,d0=1.) where {F1}
            setrefval(d,$df(frefval1(d)))
            d.d1(y,x,p,d0*refval(d))
            return 
        end
        
        @inline function (h::Hessian11{typeof($f),H1,H11})(z,x,p=nothing,h0=1) where {H1,H11}
            h.h1(z,x,p,h0*refval(h))
            h.h11(z,x,p,h0*$ddf(frefval(h)))
            return
        end
        
        @inline non_caching_eval(e::Expression1{T,typeof($f),E},x,p=nothing) where {T,E} = $f(non_caching_eval(e.e1,x,p))
        @inline function non_caching_eval(d::Gradient1{typeof($f),D1},y,x,p=nothing,d0=1.)  where D1
            non_caching_eval(d.d1,y,x,p,d0*$df(frefval1(d)))
            return 
        end
        # @inline non_caching_grad(e::Expression1{typeof($f),F},y,x,p=nothing,d0=1.) where F =
        #     non_caching_grad(e.e1,y,x,p,d0*$df(e.e1(x)))
    end
end


for (f,df1,df2,ddf11,ddf12,ddf22) in f_nargs_2
    @eval begin
        $f(e1::E1,e2::E2) where {E1 <: Expression,E2 <: Expression} =
            Expression2($f,e1,e2)
        $f(e1::E1,e2::E2) where {E1 <: Expression,E2 <: Real} =
            Expression2($f,e1,e2)
        $f(e1::E1,e2::E2) where {E1 <: Real,E2 <: Expression} =
            Expression2($f,e1,e2)

        @inline (e::Expression2{T,typeof($f),F1,F2})(x,p=nothing) where {T,F1,F2} = setrefval(e,$f(e.e1(x,p),e.e2(x,p)))
        @inline (e::Expression2{T,typeof($f),F1,F2})(x,p=nothing) where {T,F1<:Real,F2} = setrefval(e,$f(e.e1,e.e2(x,p)))
        @inline (e::Expression2{T,typeof($f),F1,F2})(x,p=nothing) where {T,F1,F2<:Real} = setrefval(e,$f(e.e1(x,p),e.e2))

        @inline function (d::Gradient2F1{typeof($f),F1,R})(y,x,p=nothing,d0=1.) where {F1,R}
            setrefval(d,$df2(d.a,frefval1(d)))
            d.d1(y,x,p,d0*refval(d))
            return 
        end
        @inline function (d::Gradient2F2{typeof($f),F1,R})(y,x,p=nothing,d0=1.) where {F1,R}
            setrefval(d,$df1(frefval1(d),d.a))
            d.d1(y,x,p,d0*refval(d))
            return 
        end
        @inline function (d::Gradient2{typeof($f),D1,D2})(y,x,p=nothing,d0 = 1) where {D1,D2}
            setrefval1(d,$df1(frefval1(d),frefval2(d)))
            setrefval2(d,$df2(frefval1(d),frefval2(d)))
            d.d1(y,x,p,d0*refval1(d))
            d.d2(y,x,p,d0*refval2(d))
            return
        end
        @inline function (h::Hessian11F1{typeof($f),H1,H11,R})(z,x,p=nothing,h0=1) where {H1,H11,R}
            h.h1(z,x,p,h0*refval(h))
            h.h11(z,x,p,h0*$ddf22(h.a,frefval(h)))
            return
        end
        @inline function (h::Hessian11F2{typeof($f),H1,H11,R})(z,x,p=nothing,h0=1) where {H1,H11,R}
            h.h1(z,x,p,h0*refval(h))
            h.h11(z,x,p,h0*$ddf11(frefval(h),h.a))
            return
        end
        @inline function (h::Hessian22{typeof($f),H1,H2,H11,H12,H22})(z,x,p=nothing,h0=1) where {H1,H2,H11,H12,H22}
            ddf12 = $ddf12(frefval1(h),frefval2(h))
            h.h1(z,x,p,h0*refval1(h))
            h.h2(z,x,p,h0*refval2(h))
            h.h11(z,x,p,h0*$ddf11(frefval1(h),frefval2(h)))
            h.h12(z,x,p,h0*ddf12)
            h.h21(z,x,p,h0*ddf12)
            h.h22(z,x,p,h0*$ddf22(frefval1(h),frefval2(h)))
            return
        end

        @inline non_caching_eval(e::Expression2{T,typeof($f),F1,F2},x,p=nothing) where {T,F1,F2} = $f(non_caching_eval(e.e1,x,p),non_caching_eval(e.e2,x,p))
        @inline non_caching_eval(e::Expression2{T,typeof($f),F1,F2},x,p=nothing) where {T,F1<:Real,F2} = $f(e.e1,non_caching_eval(e.e2,x,p))
        @inline non_caching_eval(e::Expression2{T,typeof($f),F1,F2},x,p=nothing) where {T,F1,F2<:Real} = $f(non_caching_eval(e.e1,x,p),e.e2)
        
        @inline function non_caching_eval(d::Gradient2F1{typeof($f),D1,R},y,x,p=nothing,d0=1.) where {D1,R}
            non_caching_eval(d.d1,y,x,p,d0*$df2(d.a,frefval1(d)))
            return 
        end
        @inline function non_caching_eval(d::Gradient2F2{typeof($f),D1,R},y,x,p=nothing,d0=1.) where {D1,R}
            non_caching_eval(d.d1,y,x,p,d0*$df1(frefval1(d),d.a))
            return 
        end
        @inline function non_caching_eval(d::Gradient2{typeof($f),D1,D2},y,x,p=nothing,d0 = 1) where {D1,D2}
            non_caching_eval(d.d1,y,x,p,d0*$df1(frefval1(d),frefval2(d)))
            non_caching_eval(d.d2,y,x,p,d0*$df2(frefval1(d),frefval2(d)))
            return
        end

        # @inline function non_caching_grad(e::Expression2{typeof($f),F1,F2},y,x,p=nothing,d0=1.) where {F1,F2}
        #     f1 = e.e1(x,p)
        #     f2 = e.e2(x,p)
        #     non_caching_grad(e.e1,y,x,p,d0*$df1(f1,f2))
        #     non_caching_grad(e.e2,y,x,p,d0*$df2(f1,f2))
        #     return
        # end
        # @inline non_caching_grad(e::Expression2{typeof($f),F1,F2},y,x,p=nothing,d0=1.) where {F1<:Real,F2} = non_caching_grad(e.e2,y,x,p,d0*$df2(e.e1,e.e2(x)))
        # @inline non_caching_grad(e::Expression2{typeof($f),F1,F2},y,x,p=nothing,d0=1.) where {F1,F2<:Real} = non_caching_grad(e.e1,y,x,p,d0*$df1(e.e1(x),e.e2))
    end
end

add_sum(e1::E,e2) where {T <: AbstractFloat, E <: Expression{T}} = add_sum(ExpressionSum([e1]),e2)
add_sum(e1::ExpressionSum{T,E,I},e2) where {T,E,I} = _add_sum(e1,e2) ? e1 : ExpressionSum(e1,[e2])
function _add_sum(e1::ExpressionSum{T,E,I},e2) where {T,E,I}
    if e2 isa eltype(e1.es)
        push!(e1.es,e2)
        return true
    else
        return _add_sum(inner(e1),e2)
    end
end
function _add_sum(e1::ExpressionSum{T,E,Nothing},e2) where {T,E}
    if e2 isa eltype(e1.es)
        push!(e1.es,e2)
        return true
    else
        return false
    end
end
