"""
Expression(ex::MOI.Nonlinear.Expression; subex = nothing)

Create a `MadDiff.Expression` from `MOI.Expression`.

"""
MadDiffCore.Expression(ex::MOI.Nonlinear.Expression; subex = nothing) =
    _convert_expression(ex::MOI.Nonlinear.Expression, 1, -1; subex=subex)[1]
MadDiffCore.Expression(::Nothing; subex = nothing) = MadDiffCore.ExpressionNull()


function _convert_expression(ex::MOI.Nonlinear.Expression, i::Int, p::Int; subex = nothing)
    j = i + 1
    node = ex.nodes[i]
    typ = node.type
    
    if node.parent != p
        error("parent node does not match")
    end
    
    if typ == MOI.Nonlinear.NODE_CALL_MULTIVARIATE
        ex1, j = _convert_expression(ex, j, i; subex=subex)
        ex2, j = _convert_expression(ex, j, i; subex=subex)
        
        if _is_another_child(ex,i,j)
            ex3, j = _convert_expression(ex, j, i; subex=subex)
            if _is_another_child(ex,i,j)
                ex4, j = _convert_expression(ex, j, i; subex=subex)
                if _is_another_child(ex,i,j)
                    exs = Any[ex1,ex2,ex3,ex4]
                    while _is_another_child(ex,i,j)
                        ex1, j = _convert_expression(ex, j, i; subex=subex)
                        push!(exs,ex1)
                    end
                    return node.index == 1 ? (sum(e for e in exs), j) :
                        node.index == 3 ? (prod(e for e in exs), j) :
                        (_get_multivariate_fun(node.index)(exs...), j)
                else
                    return _get_multivariate_fun(node.index)(ex1,ex2,ex3,ex4), j
                end
            else
                return _get_multivariate_fun(node.index)(ex1,ex2,ex3), j
            end
        else
            return _get_multivariate_fun(node.index)(ex1,ex2), j
        end
            
    elseif typ == MOI.Nonlinear.NODE_CALL_UNIVARIATE
        ex1, j = _convert_expression(ex,j,i; subex=subex)
        return _get_univariate_fun(node.index)(ex1), j
        
    elseif typ == MOI.Nonlinear.NODE_MOI_VARIABLE
        return X[node.index], j
        
    elseif typ == MOI.Nonlinear.NODE_PARAMETER
        return P[node.index], j
        
    elseif typ == MOI.Nonlinear.NODE_VALUE
        v = ex.values[node.index]
        return isinteger(v) ? Int(v) : v, j
        
    elseif typ == MOI.Nonlinear.NODE_SUBEXPRESSION
        return subex[node.index], j
        
    elseif typ == MOI.Nonlinear.NODE_COMPARISON
        ex1, j = _convert_expression(ex, j, i; subex=subex)
        ex2, j = _convert_expression(ex, j, i; subex=subex)
        return _get_comparison(node.index, ex1, ex2), j
        
    elseif typ == MOI.Nonlinear.NODE_LOGIC
        ex1, j = _convert_expression(ex, j, i; subex=subex)
        ex2, j = _convert_expression(ex, j, i; subex=subex)
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

_is_another_child(ex,i,j) = (j <= length(ex.nodes) && ex.nodes[j].parent == i)

function _get_comparison(index,ex1,ex2)

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

"""
    MadDiffCore.SparseNLPCore(nlp_data::MOI.Nonlinear.Model)
Create `MadDiffCore.SparseNLPCore` from `MOI.Nonlinear.Model`.
"""
function MadDiffCore.SparseNLPCore(nlp_data::MOI.Nonlinear.Model)
    
    subex = MadDiffCore.Field()
    con = MadDiffCore.Field()
    @simd for i=1:length(nlp_data.expressions)
        @inbounds subex[i] = MadDiffCore.Expression(nlp_data.expressions[i])
    end

    for (key,val) in nlp_data.constraints
        @inbounds con[key.value] = MadDiffCore.Expression(val.expression; subex = subex)
    end
    
    obj = MadDiffCore.Expression(nlp_data.objective; subex = subex)

    return MadDiffCore.SparseNLPCore(obj,con)
end

MOI.NLPBoundsPair(constraints) = [MOI.NLPBoundsPair(val.set) for (key,val) in constraints]
MOI.NLPBoundsPair(set::MOI.LessThan) = MOI.NLPBoundsPair(-Inf,set.upper)
MOI.NLPBoundsPair(set::MOI.GreaterThan) = MOI.NLPBoundsPair(set.lower, Inf)
MOI.NLPBoundsPair(set::MOI.EqualTo) = MOI.NLPBoundsPair(set.value,set.value)
MOI.NLPBoundsPair(set::MOI.Interval) = MOI.NLPBoundsPair(set.lower,set.upper)

"""
    MadDiffAD() <: MOI.Nonlinear.AbstractAutomaticDifferentiation

A differentiation backend for MathOptInterface based on MadDiff
"""
Base.@kwdef struct MadDiffAD <: MOI.Nonlinear.AbstractAutomaticDifferentiation
    threaded::Bool = false
end

"""
    MadDiffEvaluator <: MOI.AbstractNLPEvaluator
A type for callbacks for `MathOptInterface`'s nonlinear model.
"""
mutable struct MadDiffEvaluator <: MOI.AbstractNLPEvaluator
    threaded::Bool
    backend::MadDiffCore.AbstractNLPCore
    bounds::Vector{MOI.NLPBoundsPair}
    parameters::Vector{Float64}
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64
end

"""
    MOI.Nonlinear.Evaluator(model::MOI.Nonlinear.Model, ::MadDiffAD, ::Vector{MOI.VariableIndex})
Create a `MOI.Nonlinear.Evaluator` from `MOI.Nonlinear.Model` using `MadDiff`'s AD capability.
"""
function MOI.Nonlinear.Evaluator(model::MOI.Nonlinear.Model, ad::MadDiffAD, ::Vector{MOI.VariableIndex})
    backend = MadDiffCore.SparseNLPCore(model)
    bounds  = MOI.NLPBoundsPair(model.constraints)
    MadDiffEvaluator(
        ad.threaded,backend::MadDiffCore.AbstractNLPCore,bounds,model.parameters,
        0.,0.,0.,0.,0.)
end


"""
    MOI.eval_objective(evaluator::MadDiffEvaluator, x)
Return the objective value of `evaluator` at `x`.
"""
function MOI.eval_objective(evaluator::MadDiffEvaluator, x)
    start = time()
    obj = MadDiffCore.obj(evaluator.backend, x, evaluator.parameters)
    evaluator.eval_objective_timer += time() - start
    return obj
end
"""
    MOI.eval_objective_gradient(evaluator::MadDiffEvaluator, g, x)
Evaluate the constraints of `evaluator` at `x` and store the result in the vector `g`. 
"""
function MOI.eval_objective_gradient(evaluator::MadDiffEvaluator, g, x)
    start = time()
    MadDiffCore.grad!(evaluator.backend, x, g, evaluator.parameters)
    evaluator.eval_objective_gradient_timer += time() - start
    return
end
"""
    MOI.eval_constraint(evaluator::MadDiffEvaluator, g, x)
Evaluate the gradient of `evaluator` at `x` and store the result in the vector `g`.
"""
function MOI.eval_constraint(evaluator::MadDiffEvaluator, g, x)
    start = time()
    MadDiffCore.cons!(evaluator.backend, x, g, evaluator.parameters)
    evaluator.eval_constraint_timer += time() - start
    return
end
"""
    MOI.eval_constraint_jacobian(evaluator::MadDiffEvaluator, J, x)
Evaluate the constraints Jacobian of `evaluator` at `x` and store the result in the vector `J` in sparse coordinate format.
"""
function MOI.eval_constraint_jacobian(evaluator::MadDiffEvaluator, J, x)
    start = time()
    MadDiffCore.jac_coord!(evaluator.backend, x, J, evaluator.parameters)
    evaluator.eval_constraint_jacobian_timer += time() - start
    return
end
"""
    MOI.eval_hessian_lagrangian(evaluator::MadDiffEvaluator, H, x, σ, μ)
Evaluate the Lagrangian Hessian of `evaluator` at primal `x`, dual `μ`, and objective weight `σ` and store the result in the vector `H` in sparse coordinate format.
"""
function MOI.eval_hessian_lagrangian(evaluator::MadDiffEvaluator, H, x, σ, μ)
    start = time()
    MadDiffCore.hess_coord!(evaluator.backend, x, μ, H, evaluator.parameters;
                            obj_weight = σ)
    evaluator.eval_hessian_lagrangian_timer += time() - start
    return
end

"""
    MOI.jacobian_structure(evaluator::MadDiffEvaluator)
Return the structure of the constraints Jacobian in Vector{Tuple{Int,Int}} format.
"""
function MOI.jacobian_structure(evaluator::MadDiffEvaluator)
    return evaluator.backend.jac_sparsity
end

"""
    MOI.hessian_lagrangian_structure(evaluator::MadDiffEvaluator)
Return the structure of the Lagrangian Hessian in Vector{Tuple{Int,Int}} format.
"""
function MOI.hessian_lagrangian_structure(evaluator::MadDiffEvaluator)
    return evaluator.backend.hess_sparsity
end

"""
    MOI.NLPBlockData(evaluator::MadDiffEvaluator)
Create MOI.NLPBlockData from MadDiffEvaluator
"""
function MOI.NLPBlockData(evaluator::MadDiffEvaluator)
    return MOI.NLPBlockData(
        evaluator.bounds,
        evaluator,
        !(evaluator.backend.obj isa MadDiffCore.ExpressionNull), 
    )
end

MOI.features_available(::MadDiffEvaluator) = [:Grad,:Hess,:Jac]
MOI.initialize(evaluator::MadDiffEvaluator,::Vector{Symbol}) = nothing

_get_univariate_fun(i) = _UNIDICT[i]
_get_multivariate_fun(i) = _MULTIDICT[i]
