add_sum(e1::E,e2) where {T <: AbstractFloat, RT, E <: Expression{T,RT}} = add_sum(ExpressionSum([e1]),e2)
add_sum(e1::ExpressionSum{T,RT,E,I},e2) where {T,RT,E,I} = _add_sum(e1,e2) ? e1 : ExpressionSum(e1,[e2])
add_sum(e1::ExpressionSum{T,RT,E,I},e2::R) where {T,RT,E,I,R <: Real} = add_sum(e1,Constant{T}(e2)) # this is an ad-hoc solution. need to be revised later.
function _add_sum(e1::ExpressionSum{T,RT,E,I},e2) where {T,RT,E,I}
    if e2 isa eltype(e1.es)
        push!(e1.es,e2)
        return true
    else
        return _add_sum(inner(e1),e2)
    end
end
function _add_sum(e1::ExpressionSum{T,RT,E,Nothing},e2) where {T,RT,E}
    if e2 isa eltype(e1.es)
        push!(e1.es,e2)
        return true
    else
        return false
    end
end

ifelse(e0::Bool,e1,e2) = Base.ifelse(e0,e1,e2)
ifelse(e0::E0,e1,e2) where {T,RT,E0 <: Expression{T,RT}} = ExpressionIfElse(e0,e1,e2)
