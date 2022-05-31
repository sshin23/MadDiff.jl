# Jacobian
struct JacobianEntry{T,E <: Gradient{T}} <: Entry{T}
    index::Int
    e::E
end
Jacobian(f::Field1{T,G,I},indexer = nothing) where {T,G,I} = Field1(Jacobian(inner(f),indexer),[JacobianEntry(index(ie),Gradient(ie.e,(index(ie),indexer))) for ie in f.es])
Jacobian(f::Field1{T,G,Nothing},indexer = nothing) where {T,G} = Field1(nothing,[JacobianEntry(index(ie),Gradient(ie.e,(index(ie),indexer))) for ie in f.es])
Jacobian(::MadDiffCore.FieldNull{T}) where T = FieldNull{T}()

# Lagrangian Hessian
struct LagrangianEntry{T,E <: Hessian{T}} <: Entry{T}
    index::Int
    e::E
end
Field1(inner,::Vector{LagrangianEntry{HessianNull{T}}}) where T = inner

struct LagrangianHessian{T,F <: AbstractExpression{T}}
    f::F
end

LagrangianHessian(obj,grad,con,jac,indexer=nothing) = LagrangianHessian(_LagrangianHessian(obj,grad,con,jac,indexer))

function _LagrangianHessian(obj,grad,con::Field1{T,E1,I1},jac::Field1{T,E2,I2},indexer = nothing) where {T,E1,E2,I1,I2}
    Field1(_LagrangianHessian(obj,grad,inner(con),inner(jac),indexer),
           [LagrangianEntry(index(con.es[i]),Hessian(con.es[i].e,jac.es[i].e,indexer)) for i in eachindex(con.es)])
end

function _LagrangianHessian(obj,grad,con::Field1{T,E1,Nothing},jac::Field1{T,E2,Nothing},indexer = nothing) where {T,E1,E2}
    Field1(Hessian(obj,grad,indexer),[LagrangianEntry(index(con.es[i]),Hessian(con.es[i].e,jac.es[i].e,indexer)) for i in eachindex(con.es)])
end

function _LagrangianHessian(obj,grad,con::FieldNull{T},jac::FieldNull{T},indexer = nothing) where {T,E1,E2}
    Hessian(obj,grad,indexer)
end

# NLPCore
abstract type AbstractNLPCore{T <: AbstractFloat} end

struct DenseNLPCore{T <: AbstractFloat} <: AbstractNLPCore{T}
    obj::Expression{T}
    con::Field{T}
    grad::Gradient{T}
    jac::Field{T}
    hess::LagrangianHessian{T}
end

struct SparseNLPCore{T <: AbstractFloat} <: AbstractNLPCore{T}
    obj::Expression{T}
    con::Field{T}
    grad::Gradient{T}
    jac::Field{T}
    hess::LagrangianHessian{T}
    jac_sparsity::Vector{Tuple{Int,Int}}
    hess_sparsity::Vector{Tuple{Int,Int}}
end

DenseNLPCore(obj::Expression,con::Sink) = DenseNLPCore(obj,con.inner)
function DenseNLPCore(obj::Expression,con::Field)
    grad = Gradient(obj)
    jac = Jacobian(con)
    hess = LagrangianHessian(obj,grad,con,jac)
    return DenseNLPCore(obj,con,grad,jac,hess)
end

SparseNLPCore(obj::Real,con::Sink{F}) where {T,F <: Field{T}} = SparseNLPCore(Constant{T}(obj),con)
SparseNLPCore(obj::Expression,con::Sink{F})  where {T,F <: Field{T}} = SparseNLPCore(obj,con.inner)
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
