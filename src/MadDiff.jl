module MadDiff

for (mod,meths) in [
    (
        :MadDiffCore,
        [
            :Source, :Sink, :Constant, :Expression, :Variable, :Parameter, :Field, :Gradient, :Jacobian, :Hessian, :SparseHessian, :SparseJacobian, :SparseGradient, :function_evaluator, :gradient_evaluator, :sparse_gradient_evaluator, :hessian_evaluator, :sparse_hessian_evaluator, :field_evaluator, :jacobian_evaluator, :sparse_jacobian_evaluator, :obj, :cons!, :grad!, :jac_coord!, :hess_coord!, :non_caching_eval, :default_eval, :ifelse
        ]
    ),
    (
        :MadDiffModels,
        [
            :MadDiffModel, :variable, :constraint, :objective, :parameter, :instantiate!, :value, :setvalue, :set_lower_bound, :set_upper_bound, :lower_bound, :upper_bound, :set_lower_bound, :set_upper_bound, :lower_bound, :upper_bound, :dual
        ]
    ),
    (
        :MadDiffMOI,
        [
            :MadDiffAD
        ]
    )
    ]
    
    @eval begin
        import $mod
        export $mod
    end

    for meth in meths
        @eval begin
            import $mod: $meth
            export $meth
        end
    end
end

end # module
