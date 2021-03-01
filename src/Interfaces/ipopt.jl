function IpoptProblem(
    objective::Function,constraints::Function;
    obj=objective(Variable()), cons=constraints(Variable()),
    n = get_num_variables(obj,cons...),m = length(cons),
    x = zeros(n), xl = -ones(n) * Inf, xu = ones(n) * Inf, gl = zeros(m), gu = zeros(m),
)
    
    _obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess = get_nlp_functions(obj,cons,n,m)

    prob =  createProblem(
        n,xl,xu,m,gl,gu,nnz_jac,nnz_hess,
        _obj, (x,g)->_con!(g,x), (x,f)->_grad!(f,x),
        (x, mode, I,J, jac)-> mode == :Structure ? _jac_sparsity!(I,J) : _jac!(jac,x),
        (x, mode, I,J, sig, lag, hess)-> mode == :Structure ? _hess_sparsity!(I,J) : _hess!(hess,x,lag,sig)
    )

    prob.x .= x

    return prob
end

function IpoptProblem(m::Model)
   
    _obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess = get_nlp_functions(m.obj,m.cons,m.n,m.m)

    prob =  createProblem(
        m.n,m.xl,m.xu,m.m,m.gl,m.gu,nnz_jac,nnz_hess,
        _obj, (x,g)->_con!(g,x), (x,f)->_grad!(f,x),
        (x, mode, I,J, jac)-> mode == :Structure ? _jac_sparsity!(I,J) : _jac!(jac,x),
        (x, mode, I,J, sig, lag, hess)-> mode == :Structure ? _hess_sparsity!(I,J) : _hess!(hess,x,lag,sig)
    )
    
    prob.x .= m.x

    return prob
end
