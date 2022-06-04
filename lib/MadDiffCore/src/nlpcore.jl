# Lagrangian Hessian
struct LagrangianEntry{T,RT,E <: Hessian{T,RT}} <: Entry{T,RT}
    index::Int
    e::E
end
Field1(inner,::Vector{LagrangianEntry{HessianNull{T}}}) where T = inner

struct LagrangianHessian{T,RT,F <: AbstractExpression{T,RT}}
    f::F
end

LagrangianHessian(obj,grad,con,jac,indexer=nothing) = LagrangianHessian(_LagrangianHessian(obj,grad,con,jac,indexer))

function _LagrangianHessian(obj,grad,con::Field1{T,RT,E1,I1},jac::Field1{T,RT,E2,I2},indexer = nothing) where {T,RT,E1,E2,I1,I2}
    Field1(_LagrangianHessian(obj,grad,inner(con),inner(jac),indexer),
           [LagrangianEntry(index(con.es[i]),Hessian(con.es[i].e,jac.es[i].e,indexer)) for i in eachindex(con.es)])
end

function _LagrangianHessian(obj,grad,con::Field1{T,RT,E1,Nothing},jac::Field1{T,RT,E2,Nothing},indexer = nothing) where {T,RT,E1,E2}
    Field1(Hessian(obj,grad,indexer),[LagrangianEntry(index(con.es[i]),Hessian(con.es[i].e,jac.es[i].e,indexer)) for i in eachindex(con.es)])
end

function _LagrangianHessian(obj,grad,con::FieldNull{T,RT},jac::FieldNull{T,RT},indexer = nothing) where {T,RT,E1,E2}
    Hessian(obj,grad,indexer)
end

# NLPCore
abstract type AbstractNLPCore{T <: AbstractFloat, RT <: Ref{T}} end

struct NLPCore{T <: AbstractFloat, RT <: Ref{T}} <: AbstractNLPCore{T,RT}
    obj::Expression{T,RT}
    con::Field{T,RT}
    grad::Gradient{T,RT}
    jac::Field{T,RT}
    hess::LagrangianHessian{T,RT}
end

struct SparseNLPCore{T <: AbstractFloat, RT <: Ref{T}} <: AbstractNLPCore{T,RT}
    obj::Expression{T,RT}
    con::Field{T,RT}
    grad::Gradient{T,RT}
    jac::Field{T,RT}
    hess::LagrangianHessian{T,RT}
    jac_sparsity::Vector{Tuple{Int,Int}}
    hess_sparsity::Vector{Tuple{Int,Int}}
end

NLPCore(obj::Expression,con::Sink) = NLPCore(obj,con.inner)
function NLPCore(obj::Expression,con::Field)
    grad = Gradient(obj)
    jac = Jacobian(con)
    hess = LagrangianHessian(obj,grad,con,jac)
    return NLPCore(obj,con,grad,jac,hess)
end

SparseNLPCore(obj::Real,con::Sink{F}) where {T,RT,F <: Field{T}} = SparseNLPCore(Constant{T}(obj),con)
SparseNLPCore(obj::Expression,con::Sink{F})  where {T,RT,F <: Field{T}} = SparseNLPCore(obj,con.inner)
function SparseNLPCore(obj::Expression,con::Field)
    grad = Gradient(obj)
    jac,jac_sparsity = SparseJacobian(con)
    hess,hess_sparsity = SparseLagrangianHessian(obj,grad,con,jac)
    return SparseNLPCore(obj,con,grad,jac,hess,jac_sparsity,hess_sparsity)
end

@inline function obj(nlpcore::AbstractNLPCore,x,p)
    non_caching_eval(nlpcore.obj,x,p)
end
@inline function cons!(nlpcore::AbstractNLPCore,x,y,p)
    non_caching_eval(nlpcore.con,y,x,p)
end
@inline function grad!(nlpcore::AbstractNLPCore,x,y,p)
    y .= 0
    default_eval(nlpcore.obj,x,p)
    non_caching_eval(nlpcore.grad,y,x,p)
end
@inline function jac_coord!(nlpcore::AbstractNLPCore,x,J,p)
    J .= 0
    default_eval(nlpcore.con,DUMMY,x,p)
    non_caching_eval(nlpcore.jac,J,x,p)
end
@inline function hess_coord!(nlpcore::AbstractNLPCore,x,lag,z,p; obj_weight = 1.0)
    z .= 0
    default_eval(nlpcore.obj,x,p)
    default_eval(nlpcore.con,DUMMY,x,p)
    default_eval(nlpcore.grad,DUMMY,x,p)
    default_eval(nlpcore.jac,DUMMY,x,p)
    non_caching_eval(nlpcore.hess,z,x,p,lag, obj_weight)
end
