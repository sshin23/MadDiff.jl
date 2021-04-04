import ..Ipopt

DEFAULT_OPTIMIZER = Ipopt.Optimizer

function Ipopt.Optimizer(m::Model)
    _obj,_grad!,_con!,_jac!,jac_sparsity,_hess!,hess_sparsity = nlp_evaluator(m.obj,m.con)
    _jac_fill_sparsity! = @inline (I,J)->fill_sparsity!(I,J,jac_sparsity)
    _hess_fill_sparsity! = @inline (I,J)->fill_sparsity!(I,J,hess_sparsity)
    p = m.p
    
    prob =  Ipopt.createProblem(
        m.n,m.xl,m.xu,m.m,m.gl,m.gu,length(jac_sparsity),length(hess_sparsity),
        x->_obj(x,p), (x,g)->_con!(g,x,p), (x,f)->_grad!(f,x,p),
        (x, mode, I,J, jac)-> mode == :Structure ? _jac_fill_sparsity!(I,J) : _jac!(jac,x,p),
        (x, mode, I,J, sig, lag, hess)-> mode == :Structure ? _hess_fill_sparsity!(I,J) : _hess!(hess,x,p,lag,sig)
    )

    prob.x .= m.x

    for (name,option) in m.opt
        Ipopt.addOption(prob,string(name),option)
    end

    m.x = prob.x
    m.l = prob.mult_g
    m.zl= prob.mult_x_L
    m.zu= prob.mult_x_U

    return prob
end

optimize!(prob::Ipopt.IpoptProblem) = Ipopt.solveProblem(prob)
