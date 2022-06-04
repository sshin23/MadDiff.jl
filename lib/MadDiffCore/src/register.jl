const _UNIVARIATE_FUNCTIONS = Function[]
const _BIVARIATE_FUNCTIONS = Function[]

macro register_univariate(f,df,ddf)
    return esc(
        quote
            $f(e::E) where {E <: MadDiffCore.Expression} = MadDiffCore.Expression1($f,e)
            
            @inline MadDiffCore.default_eval(e::MadDiffCore.Expression1{T,RT,typeof($f),E},x,p=nothing) where {T,RT,E} =
                MadDiffCore.setrefval(e,$f(MadDiffCore.default_eval(e.e1,x,p)))
            @inline function MadDiffCore.default_eval(d::MadDiffCore.Gradient1{T,RT,typeof($f),F1},y,x,p=nothing,d0=1.) where {T,RT,F1}
                MadDiffCore.setrefval(d,$df(MadDiffCore.frefval1(d)))
                MadDiffCore.default_eval(d.d1,y,x,p,d0*MadDiffCore.refval(d))
                return 
            end
            @inline MadDiffCore.non_caching_eval(e::MadDiffCore.Expression1{T,RT,typeof($f),E},x,p=nothing) where {T,RT,E} = $f(non_caching_eval(e.e1,x,p))
            @inline function MadDiffCore.non_caching_eval(d::MadDiffCore.Gradient1{T,RT,typeof($f),D1},y,x,p=nothing,d0=1.)  where {T,RT,D1}
                MadDiffCore.non_caching_eval(d.d1,y,x,p,d0*$df(MadDiffCore.frefval1(d)))
                return 
            end
            @inline function MadDiffCore.non_caching_eval(h::MadDiffCore.Hessian11{T,RT,typeof($f),H1,H11},z,x,p=nothing,h0=1) where {T,RT,H1,H11}
                MadDiffCore.non_caching_eval(h.h1,z,x,p,h0*MadDiffCore.refval(h))
                MadDiffCore.non_caching_eval(h.h11,z,x,p,h0*$ddf(MadDiffCore.frefval(h)))
                return
            end

            if (@__MODULE__() != Main) && (@isdefined _UNIVARIATE_FUNCTIONS)
                push!(_UNIVARIATE_FUNCTIONS,$f)
            end
            
            nothing
        end
    )
end

macro register_bivariate(f,df1,df2,ddf11,ddf12,ddf22)
    return esc(
        quote
            $f(e1::E1,e2::E2) where {T, RT, E1 <: MadDiffCore.Expression{T,RT},E2 <: MadDiffCore.Expression{T,RT}} =
                MadDiffCore.Expression2($f,e1,e2)
            $f(e1::E1,e2::E2) where {T, RT, E1 <: MadDiffCore.Expression{T,RT},E2 <: Real} =
                MadDiffCore.Expression2($f,e1,e2)
            $f(e1::E1,e2::E2) where {T, RT, E1 <: Real,E2 <: MadDiffCore.Expression{T,RT}} =
                MadDiffCore.Expression2($f,e1,e2)

            @inline MadDiffCore.default_eval(e::MadDiffCore.Expression2{T,RT,typeof($f),F1,F2},x,p=nothing) where {T,RT,F1,F2} =
                MadDiffCore.setrefval(e,$f(MadDiffCore.default_eval(e.e1,x,p),MadDiffCore.default_eval(e.e2,x,p)))
            @inline MadDiffCore.default_eval(e::MadDiffCore.Expression2{T,RT,typeof($f),F1,F2},x,p=nothing) where {T,RT,F1<:Real,F2} =
                MadDiffCore.setrefval(e,$f(e.e1,MadDiffCore.default_eval(e.e2,x,p)))
            @inline MadDiffCore.default_eval(e::MadDiffCore.Expression2{T,RT,typeof($f),F1,F2},x,p=nothing) where {T,RT,F1,F2<:Real} =
                MadDiffCore.setrefval(e,$f(MadDiffCore.default_eval(e.e1,x,p),e.e2))
            @inline function MadDiffCore.default_eval(d::MadDiffCore.Gradient2F1{T,RT,typeof($f),F1,R},y,x,p=nothing,d0=1.) where {T,RT,F1,R}
                MadDiffCore.setrefval(d,$df2(d.a,MadDiffCore.frefval1(d)))
                MadDiffCore.default_eval(d.d1,y,x,p,d0*MadDiffCore.refval(d))
                return 
            end
            @inline function MadDiffCore.default_eval(d::MadDiffCore.Gradient2F2{T,RT,typeof($f),F1,R},y,x,p=nothing,d0=1.) where {T,RT,F1,R}
                MadDiffCore.setrefval(d,$df1(MadDiffCore.frefval1(d),d.a))
                MadDiffCore.default_eval(d.d1,y,x,p,d0*MadDiffCore.refval(d))
                return 
            end
            @inline function MadDiffCore.default_eval(d::MadDiffCore.Gradient2{T,RT,typeof($f),D1,D2},y,x,p=nothing,d0 = 1) where {T,RT,D1,D2}
                MadDiffCore.setrefval1(d,$df1(MadDiffCore.frefval1(d),MadDiffCore.frefval2(d)))
                MadDiffCore.setrefval2(d,$df2(MadDiffCore.frefval1(d),MadDiffCore.frefval2(d)))
                MadDiffCore.default_eval(d.d1,y,x,p,d0*MadDiffCore.refval1(d))
                MadDiffCore.default_eval(d.d2,y,x,p,d0*MadDiffCore.refval2(d))
                return
            end
            
            @inline MadDiffCore.non_caching_eval(e::MadDiffCore.Expression2{T,RT,typeof($f),F1,F2},x,p=nothing) where {T,RT,F1,F2} = $f(MadDiffCore.non_caching_eval(e.e1,x,p),MadDiffCore.non_caching_eval(e.e2,x,p))
            @inline MadDiffCore.non_caching_eval(e::MadDiffCore.Expression2{T,RT,typeof($f),F1,F2},x,p=nothing) where {T,RT,F1<:Real,F2} = $f(e.e1,MadDiffCore.non_caching_eval(e.e2,x,p))
            @inline MadDiffCore.non_caching_eval(e::MadDiffCore.Expression2{T,RT,typeof($f),F1,F2},x,p=nothing) where {T,RT,F1,F2<:Real} = $f(MadDiffCore.non_caching_eval(e.e1,x,p),e.e2)
            
            @inline function MadDiffCore.non_caching_eval(d::MadDiffCore.Gradient2F1{T,RT,typeof($f),D1,R},y,x,p=nothing,d0=1.) where {T,RT,D1,R}
                MadDiffCore.non_caching_eval(d.d1,y,x,p,d0*$df2(d.a,MadDiffCore.frefval1(d)))
                return 
            end
            @inline function MadDiffCore.non_caching_eval(d::MadDiffCore.Gradient2F2{T,RT,typeof($f),D1,R},y,x,p=nothing,d0=1.) where {T,RT,D1,R}
                MadDiffCore.non_caching_eval(d.d1,y,x,p,d0*$df1(MadDiffCore.frefval1(d),d.a))
                return 
            end
            @inline function MadDiffCore.non_caching_eval(d::MadDiffCore.Gradient2{T,RT,typeof($f),D1,D2},y,x,p=nothing,d0 = 1) where {T,RT,D1,D2}
                MadDiffCore.non_caching_eval(d.d1,y,x,p,d0*$df1(MadDiffCore.frefval1(d),MadDiffCore.frefval2(d)))
                MadDiffCore.non_caching_eval(d.d2,y,x,p,d0*$df2(MadDiffCore.frefval1(d),MadDiffCore.frefval2(d)))
                return
            end
            @inline function MadDiffCore.non_caching_eval(h::MadDiffCore.Hessian11F1{T,RT,typeof($f),H1,H11,R},z,x,p=nothing,h0=1) where {T,RT,H1,H11,R}
                MadDiffCore.non_caching_eval(h.h1,z,x,p,h0*MadDiffCore.refval(h))
                MadDiffCore.non_caching_eval(h.h11,z,x,p,h0*$ddf22(h.a,MadDiffCore.frefval(h)))
                return
            end
            @inline function MadDiffCore.non_caching_eval(h::MadDiffCore.Hessian11F2{T,RT,typeof($f),H1,H11,R},z,x,p=nothing,h0=1) where {T,RT,H1,H11,R}
                MadDiffCore.non_caching_eval(h.h1,z,x,p,h0*MadDiffCore.refval(h))
                MadDiffCore.non_caching_eval(h.h11,z,x,p,h0*$ddf11(MadDiffCore.frefval(h),h.a))
                return
            end
            @inline function MadDiffCore.non_caching_eval(h::MadDiffCore.Hessian22{T,RT,typeof($f),H1,H2,H11,H12,H22},z,x,p=nothing,h0=1) where {T,RT,H1,H2,H11,H12,H22}
                ddf12 = $ddf12(MadDiffCore.frefval1(h),MadDiffCore.frefval2(h))
                MadDiffCore.non_caching_eval(h.h1,z,x,p,h0*MadDiffCore.refval1(h))
                MadDiffCore.non_caching_eval(h.h2,z,x,p,h0*MadDiffCore.refval2(h))
                MadDiffCore.non_caching_eval(h.h11,z,x,p,h0*$ddf11(MadDiffCore.frefval1(h),MadDiffCore.frefval2(h)))
                MadDiffCore.non_caching_eval(h.h12,z,x,p,h0*ddf12)
                MadDiffCore.non_caching_eval(h.h21,z,x,p,h0*ddf12)
                MadDiffCore.non_caching_eval(h.h22,z,x,p,h0*$ddf22(MadDiffCore.frefval1(h),MadDiffCore.frefval2(h)))
                return
            end
            
            if (@__MODULE__() != Main) && (@isdefined _BIVARIATE_FUNCTIONS)
                push!(_BIVARIATE_FUNCTIONS,$f)
            end
            
            nothing
        end
    )
end

