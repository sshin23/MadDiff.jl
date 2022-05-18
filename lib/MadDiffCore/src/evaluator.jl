function_evaluator(f::E) where E <: Expression = @inline (x,p=nothing) -> non_caching_eval(f,x,p)
field_evaluator(f::F) where F <: Field = @inline (y,x,p=nothing) -> non_caching_eval(f,y,x,p)

function gradient_evaluator(f::Expression{T}) where T <: AbstractFloat
    d = Gradient(f)
    @inline function (y,x,p=nothing)
        y.=0
        f(x,p)
        non_caching_eval(d,y,x,p)
    end
end
function sparse_gradient_evaluator(f::Expression{T}) where T <: AbstractFloat
    d,sparsity = SparseGradient(f)
    (@inline function (y,x,p=nothing)
        y.=0
        f(x,p)
        non_caching_eval(d,y,x,p)
    end), sparsity
end
function hessian_evaluator(f::Expression{T}) where T <: AbstractFloat
    d = Gradient(f)
    h = Hessian(f,d)

    @inline function (z,x,p=nothing)
        z.=0
        f(x,p)
        d(DUMMY,x,p)
        h(z,x,p)
    end
end
function sparse_hessian_evaluator(f::Expression{T}) where T <: AbstractFloat
    d = Gradient(f)
    h,sparsity = SparseHessian(f,d)

    (@inline function (z,x,p=nothing)
        z.=0
        f(x,p)
        d(DUMMY,x,p)
        h(z,x,p)
    end), sparsity
end
function sparse_jacobian_evaluator(f::Field)
    jac,sparsity = SparseJacobian(f)
    
    (@inline function (y,x,p=nothing)
        y.=0
        f(DUMMY,x,p)
        non_caching_eval(jac,y,x,p)
    end), sparsity
end
function jacobian_evaluator(f::Field)
    jac = Jacobian(f)
    
    @inline function (y,x,p=nothing)
        y.=0
        f(DUMMY,x,p)
        non_caching_eval(jac,y,x,p)
    end
end

field_evaluator(f::Sink) = field_evaluator(inner(f))
jacobian_evaluator(f::Sink) = jacobian_evaluator(inner(f))
sparse_jacobian_evaluator(f::Sink) = sparse_jacobian_evaluator(inner(f))

