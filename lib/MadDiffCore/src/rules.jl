for (f0,f,df,ddf) in _F_NARGS_1
    @eval $f0(e::E) where {E <: Expression} = Expression1($f,e)
end

for (f0,f,df,ddf) in _F_BASE
    @eval $f0(e::E) where {E <: Expression} = Expression1($f,e)
end

for (f0,f,df1,df2,ddf11,ddf12,ddf22) in _F_NARGS_2
    @eval begin
        $f0(e1::E1,e2::E2) where {E1 <: Expression,E2 <: Expression} =
            Expression2($f,e1,e2)
        $f0(e1::E1,e2::E2) where {E1 <: Expression,E2 <: Real} =
            Expression2($f,e1,e2)
        $f0(e1::E1,e2::E2) where {E1 <: Real,E2 <: Expression} =
            Expression2($f,e1,e2)
    end
end


add_sum(e1::E,e2) where {T <: AbstractFloat, E <: Expression{T}} = add_sum(ExpressionSum([e1]),e2)
add_sum(e1::ExpressionSum{T,E,I},e2) where {T,E,I} = _add_sum(e1,e2) ? e1 : ExpressionSum(e1,[e2])
add_sum(e1::ExpressionSum{T,E,I},e2::R) where {T,E,I,R <: Real} = add_sum(e1,Constant(e2)) # this is an ad-hoc solution. need to be revised later.
function _add_sum(e1::ExpressionSum{T,E,I},e2) where {T,E,I}
    if e2 isa eltype(e1.es)
        push!(e1.es,e2)
        return true
    else
        return _add_sum(inner(e1),e2)
    end
end
function _add_sum(e1::ExpressionSum{T,E,Nothing},e2) where {T,E}
    if e2 isa eltype(e1.es)
        push!(e1.es,e2)
        return true
    else
        return false
    end
end

ifelse(e0::E0,e1,e2) where {T,E0 <: Expression{T}} = ExpressionIfElse(e0,e1,e2)
