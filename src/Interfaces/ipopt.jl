import ..Ipopt

function create_problem(m::Model,::Type{Ipopt.Optimizer};opt...)
    _obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess = get_nlp_functions(m.obj,m.cons,m.p)

    prob =  Ipopt.createProblem(
        m.n,m.xl,m.xu,m.m,m.gl,m.gu,nnz_jac,nnz_hess,
        _obj, (x,g)->_con!(g,x), (x,f)->_grad!(f,x),
        (x, mode, I,J, jac)-> mode == :Structure ? _jac_sparsity!(I,J) : _jac!(jac,x),
        (x, mode, I,J, sig, lag, hess)-> mode == :Structure ? _hess_sparsity!(I,J) : _hess!(hess,x,lag,sig)
    )

    prob.x .= m.x

    for (name,option) in opt
        Ipopt.addOption(prob,string(name),option)
    end

    m.x = prob.x
    m.l = prob.mult_g
    m.zl= prob.mult_x_L
    m.zu= prob.mult_x_U

    return prob
end

solve_problem(m,::Type{Ipopt.Optimizer};opt...) = Ipopt.solveProblem(m)

