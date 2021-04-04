struct LagrangianEntry{E} <: Entry
    index::Int
    e::E
end
struct LagrangianHessian{F}
    f::F
end
@inline (laghess::LagrangianHessian)(z,x,p,l,s) = _eval_laghess(laghess.f,z,x,p,l,s)
@inline function _eval_laghess(f::Field1{E,I},z,x,p,l,s) where {E,I}
    _eval_laghess(inner(f),z,x,p,l,s)
    @simd for i in eachindex(f.es)
        f.es[i](z,x,p,l)
    end
end
@inline _eval_laghess(h::H,z,x,p,l,s) where H <: Hessian = h(z,x,p,s)
@inline (lagentry::LagrangianEntry{E})(z,x,p,l) where E = lagentry.e(z,x,p,l[index(lagentry)])

LagrangianHessian(obj,grad,con,jac,indexer=nothing) = LagrangianHessian(_LagrangianHessian(obj,grad,con,jac,indexer))

function _LagrangianHessian(obj,grad,con::Field1{E1,I1},jac::Field1{E2,I2},indexer = nothing) where {E1,E2,I1,I2}
    Field1(_LagrangianHessian(obj,grad,inner(con),inner(jac),indexer),
           [LagrangianEntry(index(con.es[i]),Hessian(con.es[i].e,jac.es[i].e,indexer)) for i in eachindex(con.es)])
end

function _LagrangianHessian(obj,grad,con::Field1{E1,Nothing},jac::Field1{E2,Nothing},indexer = nothing) where {E1,E2}
    Field1(Hessian(obj,grad,indexer),[LagrangianEntry(index(con.es[i]),Hessian(con.es[i].e,jac.es[i].e,indexer)) for i in eachindex(con.es)])
end

function _LagrangianHessian(obj,grad,con::FieldNull,jac::FieldNull,indexer = nothing) where {E1,E2}
    Hessian(obj,grad,indexer)
end
