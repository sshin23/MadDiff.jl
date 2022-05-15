const MOI = MathOptInterface
const X = MadDiffCore.Variable()
const P = MadDiffCore.Parameter()

import Base: sum, getindex
getindex(f::MadDiffCore.Sink,i) = f.inner[i]
function getindex(f::MadDiffCore.Field,i)
    for ie in f.es
        MadDiffCore.index(ie) == i && return ie.e
    end
    return f.inner[i]
end
sum(f::MadDiffCore.Sink) = sum(f.inner)
sum(f::MadDiffCore.Field) = sum(ie.e for ie in f.es) + (f.inner == nothing ? 0. : sum(f.inner))

is_another_child(ex,i,j) = j <= length(ex.nodes) && ex.nodes[j].parent == i
Expression(::Nothing; subex = nothing) = MadDiffCore.Constant(0.)
Expression(ex::MOI.Nonlinear.Expression; subex = nothing) =
    Expression(ex::MOI.Nonlinear.Expression, 1, -1; subex=subex)[1]
function Expression(ex::MOI.Nonlinear.Expression, i::Int, p::Int; subex = nothing)
    j = i + 1
    node = ex.nodes[i]
    typ = node.type
    
    if node.parent != p
        error("parent node does not match")
    end
    
    if typ == MOI.Nonlinear.NODE_CALL_MULTIVARIATE
        ex1, j = Expression(ex, j, i; subex=subex)
        ex2, j = Expression(ex, j, i; subex=subex)
        
        if is_another_child(ex,i,j)
            ex3, j = Expression(ex, j, i; subex=subex)
            if is_another_child(ex,i,j)
                ex4, j = Expression(ex, j, i; subex=subex)
                if is_another_child(ex,i,j)
                    exs = Any[ex1,ex2,ex3,ex4]
                    while is_another_child(ex,i,j)
                        ex1, j = Expression(ex, j, i; subex=subex)
                        push!(exs,ex1)
                    end

                    # if length(exs)<=6
                    #     return get_multivariate_fun(node.index)(exs...), j
                    # else
                    return node.index == 1 ? (sum(exs), j) :
                        node.index == 3 ? (prod(exs), j) :
                        (get_multivariate_fun(node.index)(exs...), j)
                    # end
                else
                    return get_multivariate_fun(node.index)(ex1,ex2,ex3,ex4), j
                end
            else
                return get_multivariate_fun(node.index)(ex1,ex2,ex3), j
            end
        else
            return get_multivariate_fun(node.index)(ex1,ex2), j
        end
            
    elseif typ == MOI.Nonlinear.NODE_CALL_UNIVARIATE
        ex1, j = Expression(ex,j,i; subex=subex)
        return get_univariate_fun(node.index)(ex1), j
        
    elseif typ == MOI.Nonlinear.NODE_MOI_VARIABLE
        return X[node.index], j
        
    elseif typ == MOI.Nonlinear.NODE_PARAMETER
        return P[node.index], j
        
    elseif typ == MOI.Nonlinear.NODE_VALUE
        return ex.values[node.index], j
        
    elseif typ == MOI.Nonlinear.NODE_SUBEXPRESSION
        return subex[node.index], j
        
    elseif typ == MOI.Nonlinear.NODE_COMPARISON
        ex1, j = Expression(ex, j, i; subex=subex)
        ex2, j = Expression(ex, j, i; subex=subex)
        return get_comparison(node.index, ex1, ex2), j
        
    elseif typ == MOI.Nonlinear.NODE_LOGIC
        ex1, j = Expression(ex, j, i; subex=subex)
        ex2, j = Expression(ex, j, i; subex=subex)
        if node.index == 1
            return MadDiffCore.and(ex1,ex2), j
        elseif node.index == 2
            return MadDiffCore.or(ex1,ex2), j
        else
            error("Unknown logic index")
        end
    else
        # NODE_VARIABLE is only used in SpasreReverseDiff
        error("Node type not supported")
    end

end

function get_comparison(index,ex1,ex2)

    if index == 1
        return ex1 <= ex2
    elseif index == 2
        return ex1 == ex2
    elseif index == 3
        return ex1 >= ex2
    elseif index == 4
        return ex1 < ex2
    elseif index == 5
        return ex1 > ex2
    else
        error("Unknown comparison index")
    end
end


function NLPCore(nlp_data::MathOptInterface.Nonlinear.Model)
    
    subex = MadDiffCore.Field()
    con = MadDiffCore.Field()
    @simd for i=1:length(nlp_data.expressions)
        @inbounds subex[i] = Expression(nlp_data.expressions[i])
    end

    for (key,val) in nlp_data.constraints
        @inbounds con[key.value] = Expression(val.expression; subex = subex)
    end
    
    obj = Expression(nlp_data.objective; subex = subex)

    return NLPCore(obj,con)
end

MOI.NLPBoundsPair(constraints) = [MOI.NLPBoundsPair(val.set) for (key,val) in constraints]
MOI.NLPBoundsPair(set::MOI.LessThan) = MOI.NLPBoundsPair(-Inf,set.upper)
MOI.NLPBoundsPair(set::MOI.GreaterThan) = MOI.NLPBoundsPair(set.lower, Inf)
MOI.NLPBoundsPair(set::MOI.EqualTo) = MOI.NLPBoundsPair(set.value,set.value)
MOI.NLPBoundsPair(set::MOI.Interval) = MOI.NLPBoundsPair(set.lower,set.upper)

struct MadDiffAD <: MOI.Nonlinear.AbstractAutomaticDifferentiation end

mutable struct MadDiffEvaluator <: MOI.AbstractNLPEvaluator
    backend
    bounds::Vector{MOI.NLPBoundsPair}
    parameters::Vector{Float64}
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64
end

function MOI.Nonlinear.Evaluator(
    model::MathOptInterface.Nonlinear.Model,
    ::MadDiffAD,
    ::Vector{MathOptInterface.VariableIndex}
    )
    backend = NLPCore(model)
    bounds  = MOI.NLPBoundsPair(model.constraints)
    MadDiffEvaluator(backend,bounds,model.parameters,0.,0.,0.,0.,0.)
end

function MOI.eval_objective(evaluator::MadDiffEvaluator, x)
    start = time()
    obj = MadDiffModels.obj(evaluator.backend, x, evaluator.parameters)
    evaluator.eval_objective_timer += time() - start
    return obj
end
function MOI.eval_objective_gradient(evaluator::MadDiffEvaluator, g, x)
    start = time()
    MadDiffModels.grad!(evaluator.backend, x, g, evaluator.parameters)
    evaluator.eval_objective_gradient_timer += time() - start
    return
end
function MOI.eval_constraint(evaluator::MadDiffEvaluator, g, x)
    start = time()
    MadDiffModels.cons!(evaluator.backend, x, g, evaluator.parameters)
    evaluator.eval_constraint_timer += time() - start
    return
end
function MOI.eval_constraint_jacobian(evaluator::MadDiffEvaluator, J, x)
    start = time()
    MadDiffModels.jac_coord!(evaluator.backend, x, J, evaluator.parameters)
    evaluator.eval_constraint_jacobian_timer += time() - start
    return
end
function MOI.eval_hessian_lagrangian(evaluator::MadDiffEvaluator, H, x, σ, μ)
    start = time()
    MadDiffModels.hess_coord!(evaluator.backend, x, μ, H, evaluator.parameters; obj_weight = σ)
    evaluator.eval_hessian_lagrangian_timer += time() - start
    return
end


function MOI.jacobian_structure(evaluator::MadDiffEvaluator)
    return evaluator.backend.jac_sparsity
end
function MOI.hessian_lagrangian_structure(evaluator::MadDiffEvaluator)
    return evaluator.backend.hess_sparsity
end

function MOI.NLPBlockData(evaluator::MadDiffEvaluator)
    return MOI.NLPBlockData(
        evaluator.bounds,
        evaluator,
        !(evaluator.backend.obj isa Constant && evaluator.backend.obj.ref.x == 0.), 
    )
end


MOI.features_available(::MadDiffEvaluator) = [:Grad,:Hess,:Jac]
MOI.initialize(evaluator::MadDiffEvaluator,::Vector{Symbol}) = nothing

const unidict = [eval(f) for f in MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS]
const multidict = [eval(f) for f in MOI.Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS]

get_univariate_fun(i) = unidict[i]
get_multivariate_fun(i) = multidict[i]
