module MadNLPOptimizer

import ..SimpleNLModels: get_nlp_functions, Model, Expression
import MadNLP: NonlinearProgram, Solver, INITIAL, optimize!

function create_problem(
    obj::Expression,cons::Vector{Expression};
    n = get_num_variables(obj,cons...),q = 0, m = length(cons),
    x = zeros(n), p = Float64[], g = zeros(m), l = zeros(m), zl = -ones(n), zu = ones(n), xl = -ones(n) * Inf, xu = ones(n) * Inf, gl = zeros(m), gu = zeros(m),opt...)
    
    _obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess = get_nlp_functions(obj,cons,p)

    return NonlinearProgram(
        n,m,nnz_hess,nnz_jac,
        0.,x,g,l,zl,zu,xl,xu,gl,gu,
        _obj,_grad!,_con!,_jac!,_hess!,
        _hess_sparsity!,_jac_sparsity!,
        INITIAL,Dict())
end

create_problem(m::Model;opt...) = create_problem(
    m.obj,m.cons;n=m.n,q=m.q,m=m.m,
    x=m.x,p=m.p,g=zeros(m.m),l=m.l,
    zl=m.zl,zu=m.zu,xl=m.xl,xu=m.xu,
    gl=m.gl,gu=m.gu)

function solve_problem(prob;opt...)
    prob.ext[:solver] = Solver(prob;opt...)
    optimize!(prob.ext[:solver])
end

end
