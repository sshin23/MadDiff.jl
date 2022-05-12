const MOI = MathOptInterface
const X = MadDiffCore.Variable()
const P = MadDiffCore.Parameter()

Expression(ex::MOI.Nonlinear.Expression) = Expression(ex::MOI.Nonlinear.Expression, 1, -1)[1]
function Expression(ex::MOI.Nonlinear.Expression, i::Int, p::Int)
    j = i + 1
    node = ex.nodes[i]
    typ = node.type
    
    if node.parent != p
        error("parent node does not match")
    end
    
    if typ == MOI.Nonlinear.NODE_CALL_MULTIVARIATE
        if node.index == 1 || node.index == 3 # needs some performance check
            exs = []
            while j <= length(ex.nodes) && ex.nodes[j].parent == i
                ex1, j = Expression(ex, j, i)
                push!(exs,ex1)
            end
            if length(exs)<=4
                return node.index == 1 ? (+(exs...), j) : (*(exs...), j)
            else
                return node.index == 1 ? (sum(exs), j) : (prod(exs), j)
            end
        # if node.index == 1 || node.index == 2 # needs some performance check
        #     exs = []
        #     while j <= length(ex.nodes) && ex.nodes[j].parent == i
        #         ex1, j = Expression(ex, j, i)
        #         push!(exs,ex1)
        #     end
        #     return get_multivariate_fun(node.index)(exs...), j
        else
            ex1, j = Expression(ex, j, i)
            ex2, j = Expression(ex, j, i)
            return get_multivariate_fun(node.index)(ex1,ex2), j
        end
    elseif typ == MOI.Nonlinear.NODE_CALL_UNIVARIATE
        ex1, j = Expression(ex,j,i)
        return get_univariate_fun(node.index)(ex1), j
        
    elseif typ == MOI.Nonlinear.NODE_MOI_VARIABLE
        return X[node.index], j
        
    elseif typ == MOI.Nonlinear.NODE_PARAMETER
        return P[node.index], j
        
    elseif typ == MOI.Nonlinear.NODE_VALUE
        return ex.values[node.index], j
        
    else
        # TODO
        # NODE_LOGIC
        # NODE_COMPARISON
        # NODE_VARIABLE
        # NODE_SUBEXPRESSION
        error("node type not supported")
    end

end

function MadDiffModel(nlp_data::MathOptInterface.Nonlinear.Model)
    m = MadDiffModel()

    for i=1:10
        variable(m)
    end
    
    if nlp_data.objective != nothing
        objective(m,Expression(nlp_data.objective))
    end

    for (key,val) in nlp_data.constraints
        if val.set isa MOI.LessThan
            constraint(m,Expression(val.expression) <= val.set.upper)
        elseif val.set isa MOI.GreaterThan
            constraint(m,Expression(val.expression) >= val.set.lower)
        elseif val.set isa MOI.EqualTo
            constraint(m,Expression(val.expression) == val.set.value)
        else
            error("constraint type not supported")
        end
    end

    return m
end

struct MadDiffAD <: MOI.Nonlinear.AbstractAutomaticDifferentiation end

mutable struct MadDiffEvaluator <: MOI.AbstractNLPEvaluator
    backend
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
    
    MadDiffEvaluator(MadDiffModel(model),0.,0.,0.,0.,0.)
end

function MOI.eval_objective(evaluator::MadDiffEvaluator, x)
    start = time()
    obj = MadDiffModels.obj(evaluator.backend, x)
    evaluator.eval_objective_timer += time() - start
    return obj
end
function MOI.eval_objective_gradient(evaluator::MadDiffEvaluator, g, x)
    start = time()
    MadDiffModels.grad!(evaluator.backend, x, g)
    evaluator.eval_objective_gradient_timer += time() - start
    return
end
function MOI.eval_constraint(evaluator::MadDiffEvaluator, g, x)
    start = time()
    MadDiffModels.cons!(evaluator.backend, x, g)
    evaluator.eval_constraint_timer += time() - start
    return
end
function MOI.eval_constraint_jacobian(evaluator::MadDiffEvaluator, J, x)
    start = time()
    MadDiffModels.jac_coord!(evaluator.backend, x, J)
    evaluator.eval_constraint_jacobian_timer += time() - start
    return
end
function MOI.eval_hessian_lagrangian(evaluator::MadDiffEvaluator, H, x, σ, μ)
    start = time()
    MadDiffModels.hess_coord!(evaluator.backend, x, μ, H; obj_weight = σ)
    evaluator.eval_hessian_lagrangian_timer += time() - start
    return
end


function MOI.jacobian_structure(evaluator::MadDiffEvaluator)
    return evaluator.backend.evaluator.jac_sparsity
end
function MOI.hessian_lagrangian_structure(evaluator::MadDiffEvaluator)
    return evaluator.backend.evaluator.hess_sparsity
end

function MOI.NLPBlockData(evaluator::MadDiffEvaluator)
    return MOI.NLPBlockData(
        [MathOptInterface.NLPBoundsPair(l,u) for (l,u) in zip(evaluator.backend.gl,evaluator.backend.gu)],
        evaluator,
        !(evaluator.backend.obj isa Constant && evaluator.backend.obj.ref.x == 0.), 
    )
end


MOI.features_available(::MadDiffEvaluator) = [:Grad,:Hess,:Jac]
MOI.initialize(evaluator::MadDiffEvaluator,::Vector{Symbol}) = MadDiffModels.instantiate!(evaluator.backend)

const unidict = [eval(f) for f in MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS]
const multidict = [eval(f) for f in MOI.Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS]

get_univariate_fun(i) = unidict[i]
get_multivariate_fun(i) = multidict[i]
