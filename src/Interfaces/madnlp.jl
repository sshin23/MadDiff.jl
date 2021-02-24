function NonlinearProgram(
    objective::Function,constraints::Function;
    obj=objective(Expression()), cons=constraints(Expression()),
    n = get_num_variables(obj,cons...),m = length(cons),
    x = zeros(n), g = zeros(m), l = zeros(m), zl = zeros(n), zu = zeros(n),
    xl = -ones(n) * Inf, xu = ones(n) * Inf, gl = zeros(m), gu = zeros(m),
    obj_val = 0., status = INITIAL, ext = Dict{Symbol,Any}()
)

    _obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess = get_nlp_functions(obj,cons,n,m)
    
    return NonlinearProgram(
        n,m,nnz_hess,nnz_jac,
        obj_val,x,g,l,zl,zu,xl,xu,gl,gu,
        _obj,_grad!,_con!,_jac!,_hess!,
        _hess_sparsity!,_jac_sparsity!,
        status,ext
    )
end
