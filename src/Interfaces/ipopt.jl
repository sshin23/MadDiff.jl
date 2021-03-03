module IpoptOptimizer

import ..SimpleNLModels: get_nlp_functions, Model, Expression
import Ipopt: createProblem, solveProblem, addOption

function createProblem(obj::Expression,cons::Vector{Expression},;
                       n = get_num_variables(obj,cons...),m = length(cons),
                       x = zeros(n), xl = -ones(n) * Inf, xu = ones(n) * Inf, gl = zeros(m), gu = zeros(m),
                       opt...)
    
    _obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess = get_nlp_functions(obj,cons)

    prob =  createProblem(
        n,xl,xu,m,gl,gu,nnz_jac,nnz_hess,
        _obj, (x,g)->_con!(g,x), (x,f)->_grad!(f,x),
        (x, mode, I,J, jac)-> mode == :Structure ? _jac_sparsity!(I,J) : _jac!(jac,x),
        (x, mode, I,J, sig, lag, hess)-> mode == :Structure ? _hess_sparsity!(I,J) : _hess!(hess,x,lag,sig)
    )

    prob.x .= x

    for (name,option) in opt
        addOption(prob,string(name),option)
    end

    return prob
end

function createProblem(m::Model) 
    prob = createProblem(m.obj,m.cons;n=m.n,m=m.m,x=m.x,xl=m.xl,xu=m.xu,gl=m.gl,gu=m.gu,m.opt...)

    m.x = prob.x
    m.l = prob.mult_g
    m.zl= prob.mult_x_L
    m.zu= prob.mult_x_U

    return prob
end

end
