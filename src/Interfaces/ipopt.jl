import ..Ipopt

function instantiate!(m::Model,::Type{Ipopt.Optimizer};opt...)
    _obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess = get_nlp_functions(m.objs,m.cons,m.p)

    m.prob =  Ipopt.createProblem(
        m.n,m.xl,m.xu,m.m,m.gl,m.gu,nnz_jac,nnz_hess,
        _obj, (x,g)->_con!(g,x), (x,f)->_grad!(f,x),
        (x, mode, I,J, jac)-> mode == :Structure ? _jac_sparsity!(I,J) : _jac!(jac,x),
        (x, mode, I,J, sig, lag, hess)-> mode == :Structure ? _hess_sparsity!(I,J) : _hess!(hess,x,lag,sig)
    )

    m.prob.x .= m.x

    for (name,option) in opt
        Ipopt.addOption(m.prob,string(name),option)
    end

    m.x = m.prob.x
    m.l = m.prob.mult_g
    m.zl= m.prob.mult_x_L
    m.zu= m.prob.mult_x_U

    return
end

optimize!(prob,::Type{Ipopt.Optimizer};opt...) = Ipopt.solveProblem(prob)

