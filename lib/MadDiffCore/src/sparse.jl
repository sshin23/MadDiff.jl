function SparseGradient(e::E) where E
    indexer = Dict{Int,Int}()
    d = Gradient(e,indexer)
    sparsity = Vector{Int}(undef,length(indexer))
    set_sparsity!(sparsity,indexer)
    return (d,sparsity)
end
SparseJacobian(s::Sink{Field}) = SparseJacobian(inner(s))
function SparseJacobian(f::Field1{G,I}) where {G,I}
    indexer = Dict{Tuple{Int,Int},Int}()
    jac = Jacobian(f,indexer)
    sparsity = Vector{Tuple{Int,Int}}(undef,length(indexer))
    set_sparsity!(sparsity,indexer)
    return (jac,sparsity)
end
SparseJacobian(f::FieldNull) = (f,Tuple{Int,Int}[])

function SparseHessian(e::E,d::G) where {E,G}
    indexer = Dict{Tuple{Int,Int},Int}()
    hess = Hessian(e,d,indexer)
    sparsity = Vector{Tuple{Int,Int}}(undef,length(indexer))
    set_sparsity!(sparsity,indexer)
    return (hess,sparsity)
end
function SparseLagrangianHessian(obj,grad,con,jac)
    indexer = Dict{Tuple{Int,Int},Int}()
    laghess = LagrangianHessian(obj,grad,con,jac,indexer)
    sparsity = Vector{Tuple{Int,Int}}(undef,length(indexer))
    set_sparsity!(sparsity,indexer)
    return (laghess,sparsity)
end

function set_sparsity!(sparsity,indexer)
    for (p,k) in indexer
        sparsity[k] = p
    end
end
function set_indexer!(indexer,index)
    if haskey(indexer,index)
        return indexer[index]
    else
        return indexer[index] = length(indexer)+1        
    end
end
function set_indexer!(indexer,index1,index2)
    if haskey(indexer,(index1,index2))
        return indexer[index1,index2]
    else
        return indexer[index1,index2] = length(indexer)+1        
    end
end
