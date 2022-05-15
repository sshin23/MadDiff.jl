const ComplexExpression = Union{Expression2{typeof(+),E1,E2} where {E1,E2},Expression2{typeof(-),E1,E2} where {E1,E2}}

function stringsum(e::ExpressionSum{E,I}) where {E,I}
    res = string(e.inner)
    @simd for i in eachindex(e.es)
        res *= ("+" * string(e.es[i]))
    end
    res
end

function stringsum(e::ExpressionSum{E,Nothing}) where E
    res = string(e.es[1])
    @simd for i in 2:length(e.es)
        res *= ("+" * string(e.es[i]))
    end
    res
end

const printlist = [
    (:(Source{Parameter}),e->"p"),
    (:(Source{Variable}),e->"x"),
    (:(Variable),e-> "x[$(index(e))]"),
    (:(Parameter),e-> "p[$(index(e))]"),
    (:(Constant),e-> "$(refval(e))"),
    (:(Expression2{typeof(+),E1,E2} where {E1,E2}),e-> "$(string(e.e1)) + $(string(e.e2))"),
    (:(Expression2{typeof(-),E1,E2} where {E1,E2}),e-> "$(string(e.e1)) - $(string(e.e2))"),
    (:(Expression2{typeof(*),E1,E2} where {E1,E2}),e-> "$(string(e.e1))*$(string(e.e2))"),
    (:(Expression2{typeof(/),E1,E2} where {E1,E2}),e-> "$(string(e.e1))/$(string(e.e2))"),
    (:(Expression2{typeof(^),E1,E2} where {E1,E2}),e-> "$(string(e.e1))^$(string(e.e2))"),
    (:(Expression2{typeof(==),E1,E2} where {E1,E2}),e-> "$(string(e.e1)) == $(string(e.e2))"),
    (:(Expression2{typeof(>=),E1,E2} where {E1,E2}),e-> "$(string(e.e1)) >= $(string(e.e2))"),
    (:(Expression2{typeof(<=),E1,E2} where {E1,E2}),e-> "$(string(e.e1)) <= $(string(e.e2))"),
    (:(Expression2{typeof(+),E1,E2} where {E1<:ComplexExpression,E2}),e-> "$(string(e.e1)) + $(string(e.e2))"),
    (:(Expression2{typeof(-),E1,E2} where {E1<:ComplexExpression,E2}),e-> "$(string(e.e1)) - $(string(e.e2))"),
    (:(Expression2{typeof(*),E1,E2} where {E1<:ComplexExpression,E2}),e-> "($(string(e.e1)))*$(string(e.e2))"),
    (:(Expression2{typeof(/),E1,E2} where {E1<:ComplexExpression,E2}),e-> "($(string(e.e1)))/$(string(e.e2))"),
    (:(Expression2{typeof(^),E1,E2} where {E1<:ComplexExpression,E2}),e-> "($(string(e.e1)))^$(string(e.e2))"),
    (:(Expression2{typeof(+),E1,E2} where {E1,E2<:ComplexExpression}),e-> "$(string(e.e1)) + ($(string(e.e2)))"),
    (:(Expression2{typeof(-),E1,E2} where {E1,E2<:ComplexExpression}),e-> "$(string(e.e1)) - ($(string(e.e2)))"),
    (:(Expression2{typeof(*),E1,E2} where {E1,E2<:ComplexExpression}),e-> "$(string(e.e1))*($(string(e.e2)))"),
    (:(Expression2{typeof(/),E1,E2} where {E1,E2<:ComplexExpression}),e-> "$(string(e.e1))/($(string(e.e2)))"),
    (:(Expression2{typeof(^),E1,E2} where {E1,E2<:ComplexExpression}),e-> "$(string(e.e1))^($(string(e.e2)))"),
    (:(Expression2{typeof(+),E1,E2} where {E1<:ComplexExpression,E2<:ComplexExpression}),e-> "$(string(e.e1)) + $(string(e.e2))"),
    (:(Expression2{typeof(-),E1,E2} where {E1<:ComplexExpression,E2<:ComplexExpression}),e-> "$(string(e.e1)) - ($(string(e.e2)))"),
    (:(Expression2{typeof(*),E1,E2} where {E1<:ComplexExpression,E2<:ComplexExpression}),e-> "($(string(e.e1)))*($(string(e.e2)))"),
    (:(Expression2{typeof(/),E1,E2} where {E1<:ComplexExpression,E2<:ComplexExpression}),e-> "($(string(e.e1)))/($(string(e.e2)))"),
    (:(Expression2{typeof(^),E1,E2} where {E1<:ComplexExpression,E2<:ComplexExpression}),e-> "($(string(e.e1)))^($(string(e.e2)))"),
    (:(ExpressionSum{E,I} where {E,I}),stringsum),
]


for (f,df,ddf) in f_nargs_1
    @eval begin
        push!(printlist,(:(Expression1{$(typeof($f)),E} where E),e->"$(string($f))($(string(e.e1)))"))
    end
end

for (f,df1,df2,ddf11,ddf12,ddf22) in f_nargs_2
    f in [:+,:-,:*,:^,:/,:(<=),:(==),:(>=)] && continue
    @eval begin
        push!(printlist,(:(Expression2{$(typeof($f)),E1,E2} where {E1,E2}),e->"$(string($f))($(string(e.e1)),$(string(e.e2)))"))
    end
end

for (T,f) in printlist
    @eval begin
        show(io::IO,e::$T) = print(io,$f(e))
    end
end


# (:(Model),e-> "NLP model with $(m.n) variables and $(m.m) constraints"
