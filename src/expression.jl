abstract type Expression end

struct Source{T}
    str::String
end

struct Variable{T} <: Expression
    parent::T
    index::Int
end

struct Parameter{T} <: Expression
    parent::T
    index::Int
end

struct Constant{T,V<:Real} <: Expression
    parent::T
    val::V
end

struct Term{T,F<:Function,R<:AbstractDict{Int,Function}} <: Expression
    parent::T
    func::F
    deriv::R
end

Variable() = Source{Variable}(DEFAULT_VAR_STRING)
Variable(str::String) = Source{Variable}(str)
Variable(n::Int;parent=nothing) = Variable(parent,n)

Parameter() = Source{Parameter}(DEFAULT_PAR_STRING)
Parameter(str::String) = Source{Parameter}(str)
Parameter(n::Int;parent=nothing) = Parameter(parent,n)

index(e::Variable) = e.index
getindex(e::Source{Variable},n) = Variable(n;parent=e)
getindex(e::Source{Parameter},n) = Parameter(n;parent=e)

fxentry(n) = @inline (x,p=nothing)->x[n]
fpentry(n) = @inline (x,p=nothing)->p[n]

parent(e1,e2) = marry(parent(e1),parent(e2))
marry(p1::Source{Variable},p2::Source{Parameter}) = (p1,p2)
marry(p1::Source{Parameter},p2::Source{Variable}) = (p2,p1)
marry(p1::Tuple,p2::Tuple) = p1 == p2 ? p1 : error("Parents of the expressions are not compatible")
marry(p1::Tuple,p2) = p2 in p1 ? p1 : error("Parents of the expressions are not compatible")
marry(p1,p2::Tuple) = p1 in p2 ? p2 : error("Parents of the expressions are not compatible")
marry(p1,p2) = p1 == p2 ? p1 : error("Parents of the expressions are not compatible")

con_one(x,p=nothing) = 1.
con_zero(x,p=nothing) = 0.
TypeConOne=typeof(con_one)
TypeConZero=typeof(con_zero)
zero(e::E) where E <: Expression = 0.
one(e::E) where E <: Expression = 1.


parent(e::E) where E <: Expression = e.parent
func(e::Term) = e.func
func(e::Variable) = fxentry(e.index)
func(e::Parameter) = fpentry(e.index)
deriv(e::Term) = e.deriv
deriv(e::Variable) = ImmutableDict{Int,Function}(e.index,con_one)
deriv(e::Parameter) = ImmutableDict{Int,Function}()
func(a::T) where T <: Real = @inline (x,p=nothing)->a
deriv(e::T) where T <: Real = ImmutableDict{Int,Function}()

fsub(f::F) where F <: Function = @inline (x,p=nothing)->-f(x,p)
fadd(f1::F1,f2::F2) where {F1 <: Function, F2 <: Function} = @inline (x,p=nothing)->f1(x,p)+f2(x,p)
fmul(f1::F1,f2::F2) where {F1 <: Function, F2 <: Function} = @inline (x,p=nothing)->f1(x,p)*f2(x,p)
fdiv(f1::F1,f2::F2) where {F1 <: Function, F2 <: Function} = @inline (x,p=nothing)->f1(x,p)/f2(x,p)
fsub(f1::F1,f2::F2) where {F1 <: Function, F2 <: Function} = @inline (x,p=nothing)->f1(x,p)-f2(x,p)
fcom(f1::F1,f2::F2) where {F1 <: Function, F2 <: Function} = @inline (x,p=nothing)->f1(f2(x,p))
fcom(f1::F1,f2::F2,f3::F3) where {F1 <: Function,F2 <: Function,F3 <: Function} = @inline (x,p=nothing)->f1(f2(x,p),f3(x,p))

return_type(f,typ) = Base.return_types(f,(typ,))[1]
return_type(f,typ1,typ2) = Base.return_types(f,(typ1,typ2))[1]
fsub(f::return_type(fsub,Function)) = f.f
fsub(f::return_type(fsub,Function,Function)) = fsub(f.f2,f.f1)
fadd(f1::F,f2::TypeConZero) where F <: Function = f1
fadd(f1::TypeConZero,f2::F) where F <: Function = f2
fmul(f1::F,f2::TypeConOne) where F <: Function = f1
fmul(f1::TypeConOne,f2::F) where F <: Function = f2

fadd(f1::F,f2::T) where {T <: Real,F <: Function} = f2 == 0 ? f1 : @inline (x,p=nothing)->f1(x,p)+f2
fadd(f1::T,f2::F) where {T <: Real,F <: Function} = f1 == 0 ? f2 : @inline (x,p=nothing)->f1+f2(x,p)
fsub(f1::F,f2::T) where {T <: Real,F <: Function} = f2 == 0 ? f1 : @inline (x,p=nothing)->f1(x,p)-f2
fsub(f1::T,f2::F) where {T <: Real,F <: Function} = f1 == 0 ? fsub(f2) : @inline (x,p=nothing)->f1-f2(x,p)
fmul(f1::F,f2::T) where {T <: Real,F <: Function} = f2 == 1 ? f1 : f2 == 0 ? con_zero : @inline (x,p=nothing)->f1(x,p)*f2
fmul(f1::T,f2::F) where {T <: Real,F <: Function} = f1 == 1 ? f2 : f1 == 0 ? con_zero : @inline (x,p=nothing)->f1*f2(x,p)
fpow(f1::F,f2::T) where {T <: Real,F <: Function} = f2 == 1 ? f1 : f2 == 0 ? con_one  : @inline (x,p=nothing)->f1(x,p)^f2
fpow(f1::T,f2::F) where {T <: Real,F <: Function} = f1 == 1 ? f1 : f1 == 0 ? con_zero : @inline (x,p=nothing)->f1^f2(x,p)
fdiv(f1::F,f2::T) where {T <: Real,F <: Function} = f2 == 1 ? f1 : @inline (x,p=nothing)->f1(x,p)/f2
fdiv(f1::T,f2::F) where {T <: Real,F <: Function} = f1 == 0 ? con_zero : @inline (x,p=nothing)->f1/f2(x,p)
fcom(f1::F1,f2::F2,f3::T) where {T <: Real, F1 <: Function,F2 <: Function} = @inline (x,p=nothing)->f1(f2(x,p),f3)
fcom(f1::F1,f2::T,f3::F2) where {T <: Real, F1 <: Function,F2 <: Function} = @inline (x,p=nothing)->f1(f2,f3(x,p))
fsum(fs::Vector{F}) where F <: Function = @inline (x,p=nothing)->sum(f(x,p) for f in fs)

function f_add_sum(f1::return_type(fsum,Vector{Function}),f2::F) where F <: Function
    for f in f1.fs
        if f2 isa eltype(f.fs)
            push!(f.fs,f2)
            return f1
        end
    end
    push!(f1.fs,fsum([f2]))
    return f1
end
f_add_sum(f1::F,f2::F) where {F <: Function} = fsum(Function[fsum([f1,f2])])
f_add_sum(f1::F1,f2::F2) where {F1 <: Function,F2 <: Function} = fsum(Function[fsum([f1]),fsum([f2])])

deriv!(f,deriv) = (map!(f,values(deriv)); deriv)
deriv!(f,deriv::ImmutableDict{Int,Function}) = isempty(deriv) ? deriv : ImmutableDict{Int,Function}(deriv.key,f(deriv.value))
deriv!(f,deriv1,deriv2) = mergewith!(f,deriv1,deriv2)
deriv!(f,deriv1::ImmutableDict{Int,Function},deriv2) = mergewith(f,deriv1,deriv2)

+(e::E) where E <: Expression = e
-(e::E) where E <: Expression = Term(parent(e),fsub(func(e)),deriv!(fsub,deriv(e)))
+(e1::E1,e2::E2) where {E1<:Expression,E2<:Expression} = Term(parent(e1,e2),fadd(func(e1),func(e2)),deriv!(fadd,deriv(e1),deriv(e2)))
-(e1::E1,e2::E2) where {E1<:Expression,E2<:Expression} = Term(parent(e1,e2),fsub(func(e1),func(e2)),deriv!(fadd,deriv(e1),deriv!(fsub,deriv(e2))))
*(e1::E1,e2::E2) where {E1<:Expression,E2<:Expression} = Term(parent(e1,e2),fmul(func(e1),func(e2)),deriv!(fadd,deriv!(d->fmul(d,func(e2)),deriv(e1)),deriv!(d->fmul(d,func(e1)),deriv(e2))))
+(e::E,a::T) where {T <: Real, E <: Expression} = Term(parent(e),fadd(func(e),a),deriv(e))
+(a::T,e::E) where {T <: Real, E <: Expression} = Term(parent(e),fadd(a,func(e)),deriv(e))
-(e::E,a::T) where {T <: Real, E <: Expression} = Term(parent(e),fsub(func(e),a),deriv(e))
-(a::T,e::E) where {T <: Real, E <: Expression} = Term(parent(e),fsub(a,func(e)),deriv!(d->fsub(d),deriv(e)))
*(e::E,a::T) where {T <: Real, E <: Expression} = Term(parent(e),fmul(func(e),a),deriv!(d->fmul(d,a),deriv(e)))
*(a::T,e::E) where {T <: Real, E <: Expression} = Term(parent(e),fmul(a,func(e)),deriv!(d->fmul(a,d),deriv(e)))
^(e::E,a::T) where {T <: Real, E <: Expression} = Term(parent(e),fpow(func(e),a),deriv!(d->fmul(a,fmul(fpow(func(e),a-1),d)),deriv(e)))
^(a::T,e::E) where {T <: Real, E <: Expression} = Term(parent(e),fpow(a::T,func(e)),deriv!(d->fmul(fmul(fpow(a,func(e)),d),log(a)),deriv(e)))
/(e::E,a::T) where {T <: Real, E <: Expression} = Term(parent(e),fdiv(func(e),a),deriv!(d->fdiv(d,a),deriv(e)))
/(a::T,e::E) where {T <: Real, E <: Expression} = Term(parent(e),fdiv(a,func(e)),deriv!(d->fmul(-a,fdiv(d,fpow(func(e),2))),deriv(e)))

add_sum(e1::E1,e2::E2) where {E1<:Expression,E2<:Expression} = Term(parent(e1,e2),f_add_sum(func(e1),func(e2)),deriv!(fadd,deriv(e1),deriv(e2)))
f_caching(f1,f2) = d->((x,p=nothing)->d(x,p,f1(x,p),f2(x,p)))

for (M,f,nargs) in diffrules
    if nargs == 1
        df = :(@inline (x,p=nothing)->$(diffrule(M,f,:x)))
        @eval begin
            $f(e::E) where E <: Expression = Term(parent(e),fcom($f,func(e)),deriv!(d->fmul(fcom($df,func(e)),d),deriv(e)))
        end
    elseif nargs == 2
        dfsym,~ = diffrule(M,f,:x,:a)
        df = :(@inline (x,a,p=nothing)->$dfsym)
        @eval begin
            $f(e::E,a::T) where {T <: Real, E <: Expression} = Term(parent(e),fcom($f,func(e),a),deriv!(d->fmul(fcom($df,func(e),a),d),deriv(e)))
        end

        ~,dfsym = diffrule(M,f,:a,:x)
        df = :(@inline (a,x,p=nothing)->$dfsym)
        @eval begin
            $f(a::T,e::E) where {T <: Real, E <: Expression} = Term(parent(e),fcom($f,a,func(e)),deriv!(d->fmul(fcom($df,a,func(e)),d),deriv(e)))
        end
    end
end

for (M,f,nargs) in union(diffrules,[(:Base,:^,2),(:Base,:/,2)])
    if nargs == 2
        dfsym1,dfsym2 = diffrule(M,f,:x1,:x2)
        df1 = :(@inline (x1,x2,p=nothing)->$dfsym1)
        df2 = :(@inline (x1,x2,p=nothing)->$dfsym2)

        @eval begin
            $f(e1::E1,e2::E2) where {E1<:Expression,E2<:Expression} = Term(
                parent(e1,e2),fcom($f,func(e1),func(e2)),
                deriv!(
                    f_caching(func(e1),func(e2)),
                    deriv!(
                        (d1,d2)->((x,p,f1,f2)->d1(x,p,f1,f2)+d2(x,p,f1,f2)),
                        deriv!(d->((x,p,f1,f2) -> $df1(f1,f2)*d(x,p)),deriv(e1)),
                        deriv!(d->((x,p,f1,f2) -> $df2(f1,f2)*d(x,p)),deriv(e2))
                    )
                )
            )
        end
    end
end

