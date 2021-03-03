module MadNLPOptimizer

import ..SimpleNLModels: get_nlp_functions, Model, Expression
import MadNLP: NonlinearProgram, Solver, INITIAL, optimize!

function createProblem(
    obj::Expression,cons::Vector{Expression};
    n = get_num_variables(obj,cons...),m = length(cons),
    x = zeros(n), g = zeros(m), l = zeros(m), zl = -ones(n), zu = ones(n), xl = -ones(n) * Inf, xu = ones(n) * Inf, gl = zeros(m), gu = zeros(m))
    
    _obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess = get_nlp_functions(obj,cons)

    return NonlinearProgram(
        n,m,nnz_hess,nnz_jac,
        0.,x,g,l,zl,zu,xl,xu,gl,gu,
        _obj,_grad!,_con!,_jac!,_hess!,
        _hess_sparsity!,_jac_sparsity!,
        INITIAL,Dict())
end

createProblem(m::Model) = createProblem(
    m.obj,m.cons;n=m.n,m=m.m,
    x=m.x,g=zeros(m.m),l=m.l,
    zl=m.zl,zu=m.zu,xl=m.xl,xu=m.xu,
    gl=m.gl,gu=m.gu)

function solveProblem(prob)
    prob.ext[:solver] = Solver(prob)
    optimize!(prob.ext[:solver])
end

end
