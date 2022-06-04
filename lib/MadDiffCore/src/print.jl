const ComplexExpression = Union{Expression2{T,typeof(+),E1,E2} where {T,E1,E2},Expression2{T,typeof(-),E1,E2} where {T,E1,E2}}

function stringsum(e::ExpressionSum{T,E,I}) where {T,E,I}
    res = string(e.inner)
    @simd for i in eachindex(e.es)
        res *= ("+" * string(e.es[i]))
    end
    res
end

function stringsum(e::ExpressionSum{T,E,Nothing}) where {T,E}
    res = string(e.es[1])
    @simd for i in 2:length(e.es)
        res *= ("+" * string(e.es[i]))
    end
    res
end

const printlist = [
    (:(Source{Parameter{T}} where T),e->"p"),
    (:(Source{Variable{T}} where T),e->"x"),
    (:(AbstractVariable{T} where T),e-> "x[$(index(e))]"),
    (:(AbstractParameter{T} where T),e-> "p[$(index(e))]"),
    # (:(Constant{T} where T),e-> "$(e.x)"),
    (:(Expression2{T,typeof(+),E1,E2} where {T,E1,E2}),e-> "$(string(e.e1)) + $(string(e.e2))"),
    (:(Expression2{T,typeof(-),E1,E2} where {T,E1,E2}),e-> "$(string(e.e1)) - $(string(e.e2))"),
    (:(Expression2{T,typeof(*),E1,E2} where {T,E1,E2}),e-> "$(string(e.e1))*$(string(e.e2))"),
    (:(Expression2{T,typeof(/),E1,E2} where {T,E1,E2}),e-> "$(string(e.e1))/$(string(e.e2))"),
    (:(Expression2{T,typeof(^),E1,E2} where {T,E1,E2}),e-> "$(string(e.e1))^$(string(e.e2))"),
    (:(Expression2{T,typeof(==),E1,E2} where {T,E1,E2}),e-> "$(string(e.e1)) == $(string(e.e2))"),
    (:(Expression2{T,typeof(>=),E1,E2} where {T,E1,E2}),e-> "$(string(e.e1)) >= $(string(e.e2))"),
    (:(Expression2{T,typeof(<=),E1,E2} where {T,E1,E2}),e-> "$(string(e.e1)) <= $(string(e.e2))"),
    (:(Expression2{T,typeof(+),E1,E2} where {T,E1<:ComplexExpression,E2}),e-> "$(string(e.e1)) + $(string(e.e2))"),
    (:(Expression2{T,typeof(-),E1,E2} where {T,E1<:ComplexExpression,E2}),e-> "$(string(e.e1)) - $(string(e.e2))"),
    (:(Expression2{T,typeof(*),E1,E2} where {T,E1<:ComplexExpression,E2}),e-> "($(string(e.e1)))*$(string(e.e2))"),
    (:(Expression2{T,typeof(/),E1,E2} where {T,E1<:ComplexExpression,E2}),e-> "($(string(e.e1)))/$(string(e.e2))"),
    (:(Expression2{T,typeof(^),E1,E2} where {T,E1<:ComplexExpression,E2}),e-> "($(string(e.e1)))^$(string(e.e2))"),
    (:(Expression2{T,typeof(+),E1,E2} where {T,E1,E2<:ComplexExpression}),e-> "$(string(e.e1)) + ($(string(e.e2)))"),
    (:(Expression2{T,typeof(-),E1,E2} where {T,E1,E2<:ComplexExpression}),e-> "$(string(e.e1)) - ($(string(e.e2)))"),
    (:(Expression2{T,typeof(*),E1,E2} where {T,E1,E2<:ComplexExpression}),e-> "$(string(e.e1))*($(string(e.e2)))"),
    (:(Expression2{T,typeof(/),E1,E2} where {T,E1,E2<:ComplexExpression}),e-> "$(string(e.e1))/($(string(e.e2)))"),
    (:(Expression2{T,typeof(^),E1,E2} where {T,E1,E2<:ComplexExpression}),e-> "$(string(e.e1))^($(string(e.e2)))"),
    (:(Expression2{T,typeof(+),E1,E2} where {T,E1<:ComplexExpression,E2<:ComplexExpression}),e-> "$(string(e.e1)) + $(string(e.e2))"),
    (:(Expression2{T,typeof(-),E1,E2} where {T,E1<:ComplexExpression,E2<:ComplexExpression}),e-> "$(string(e.e1)) - ($(string(e.e2)))"),
    (:(Expression2{T,typeof(*),E1,E2} where {T,E1<:ComplexExpression,E2<:ComplexExpression}),e-> "($(string(e.e1)))*($(string(e.e2)))"),
    (:(Expression2{T,typeof(/),E1,E2} where {T,E1<:ComplexExpression,E2<:ComplexExpression}),e-> "($(string(e.e1)))/($(string(e.e2)))"),
    (:(Expression2{T,typeof(^),E1,E2} where {T,E1<:ComplexExpression,E2<:ComplexExpression}),e-> "($(string(e.e1)))^($(string(e.e2)))"),
    (:(ExpressionSum{T,E,I} where {T,E,I}),stringsum),
]

show(io::IO,e::E) where {T, F, E <: Expression1{T,F}} = print(io,"$(string(F.instance))($(string(e.e1)))")
show(io::IO,e::E) where {T, F, E <: Expression2{T,F}} = print(io,"$(string(F.instance))($(string(e.e1)),$(string(e.e2)))")

for (T,f) in printlist
    @eval begin
        show(io::IO,e::$T) = print(io,$f(e))
    end
end

for fname in [:Evaluator, :GradientEvaluator, :SparseGradientEvaluator, :HessianEvaluator, :SparseHessianEvaluator]
    @eval begin
        function show(io::IO,e::T) where T <: $fname
            println(io,"$(string(nameof(T))):")
            show(io,e.e)
        end
    end
end
