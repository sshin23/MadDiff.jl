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
    (:(Variable{T} where T),e-> "x[$(index(e))]"),
    (:(Parameter{T} where T),e-> "p[$(index(e))]"),
    (:(Constant{T} where T),e-> "$(refval(e))"),
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


for (f,df,ddf) in f_nargs_1
    @eval begin
        push!(printlist,(:(Expression1{T,$(typeof($f)),E} where {T,E}),e->"$(string($f))($(string(e.e1)))"))
    end
end

for (f,df1,df2,ddf11,ddf12,ddf22) in f_nargs_2
    f in [:+,:-,:*,:^,:/,:(<=),:(==),:(>=)] && continue
    @eval begin
        push!(printlist,(:(Expression2{T,$(typeof($f)),E1,E2} where {T,E1,E2}),e->"$(string($f))($(string(e.e1)),$(string(e.e2)))"))
    end
end

for (T,f) in printlist
    @eval begin
        show(io::IO,e::$T) = print(io,$f(e))
    end
end


# (:(Model),e-> "NLP model with $(m.n) variables and $(m.m) constraints"
