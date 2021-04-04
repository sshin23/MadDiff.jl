import ..MadNLP 

function MadNLP.Optimizer(m::Model) 
    _obj,_grad!,_con!,_jac!,jac_sparsity,_hess!,hess_sparsity = nlp_evaluator(m.obj,m.con)
    _jac_fill_sparsity! = @inline (I,J)->fill_sparsity!(I,J,jac_sparsity)
    _hess_fill_sparsity! = @inline (I,J)->fill_sparsity!(I,J,hess_sparsity)
    p = m.p

    prob = MadNLP.NonlinearProgram(
        m.n,m.m,length(hess_sparsity),length(jac_sparsity),
        0.,m.x,m.g,m.l,m.zl,m.zu,m.xl,m.xu,m.gl,m.gu,
        x->_obj(x,p),(g,x)->_grad!(g,x,p),(c,x)->_con!(c,x,p),(jac,x)->_jac!(jac,x,p),
        (hess,x,lag,sig)->_hess!(hess,x,p,lag,sig),
        _hess_fill_sparsity!,_jac_fill_sparsity!,
        MadNLP.INITIAL,Dict())
    
    prob.ext[:solver] = MadNLP.Solver(prob;m.opt...)
    
    return prob
end


optimize!(prob::MadNLP.NonlinearProgram) = MadNLP.optimize!(prob.ext[:solver])
