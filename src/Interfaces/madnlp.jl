import ..MadNLP 

function instantiate!(m::Model,::Type{MadNLP.Optimizer};opt...) 

    _obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess = get_nlp_functions(m.objs,m.cons,m.p)
    m.prob = MadNLP.NonlinearProgram(
        m.n,m.m,nnz_hess,nnz_jac,
        0.,m.x,m.g,m.l,m.zl,m.zu,m.xl,m.xu,m.gl,m.gu,
        _obj,_grad!,_con!,_jac!,_hess!,
        _hess_sparsity!,_jac_sparsity!,
        MadNLP.INITIAL,Dict())
    
    m.prob.ext[:solver] = MadNLP.Solver(m.prob;opt...)
    
    return 
end


optimize!(prob,::Type{MadNLP.Optimizer};opt...) = MadNLP.optimize!(prob.ext[:solver])
