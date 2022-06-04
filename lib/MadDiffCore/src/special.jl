acoth(x) = abs(x) < 1.0 ? NaN : Base.acoth(x) # NaNMath-like redefinition of coth

for f0 in [
    :sqrt, :log, :log2, :log1p, :log10, :sin, :cos, :tan, :asin, :acos, :atanh, :acoth
    ]
    @eval Base.$f0(e::E) where {E <: Expression} = Expression1($f0,e)
end

for (f0,id) in [
    (:+,0),
    (:-,0),
    (:/,1),
    (:^,1)
    ]
    @eval Expression2(::typeof($f0),e1::E1,e2::E2) where {T, RT, E1 <: Expression{T,RT}, E2 <: Real} = e2 == $id ? e1 : Expression2{T,RT,typeof($f0),typeof(e1),typeof(e2)}(e1,e2,RefValue{T}())
end

for (f0,id) in [
    (:+,0),
    ]
    @eval Expression2(::typeof($f0),e1::E1,e2::E2) where {T, RT, E1 <: Real, E2 <: Expression{T,RT}} = e1 == $id ? e2 : Expression2{T,RT,typeof($f0),typeof(e1),typeof(e2)}(e1,e2,RefValue{T}())
end

Expression2(::typeof(*),e1::E1,e2::E2) where {T, RT, E1 <: Expression{T}, E2 <: Real} = e2 == 0 ? 0 : e2 == 1 ? e2 : Expression2{T,RT,typeof(*),typeof(e1),typeof(e2)}(e1,e2,RefValue{T}())
Expression2(::typeof(*),e1::E1,e2::E2) where {T, RT, E1 <: Real, E2 <: Expression{T,RT}} = e1 == 0 ? 0 : e1 == 1 ? e1 : Expression2{T,RT,typeof(*),typeof(e1),typeof(e2)}(e1,e2,RefValue{T}())

