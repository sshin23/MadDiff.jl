abstract type Expression end

struct Source{T}
    str::String
end

struct Variable{T} <: Expression
    parent::T
    func::Function
    index::Int
end

struct Parameter{T} <: Expression
    parent::T
    func::Function
    index::Int
end

mutable struct Term{T} <: Expression
    parent::T
    func::Function
    deriv::Dict{Int,Function}
end

Variable() = Source{Variable}(DEFAULT_VAR_STRING)
Parameter() = Source{Parameter}(DEFAULT_PAR_STRING)
Variable(str::String) = Source{Variable}(str)
Parameter(str::String) = Source{Parameter}(str)
Variable(n::Int;parent=nothing) = Variable(parent,(x,p=nothing)->x[n],n)
Parameter(n::Int;parent=nothing) = Parameter(parent,(x,p=nothing)->p[n],n)
getindex(e::Source{Variable},n) = Variable(n;parent=e)
getindex(e::Source{Parameter},n) = Parameter(n;parent=e)

Term(;parent=nothing) = Term(parent,con_zero,Dict{Int,Function}())
parent(e1,e2) = marry(parent(e1),parent(e2))
marry(p1::Source{Variable},p2::Source{Parameter}) = (p1,p2)
marry(p1::Source{Parameter},p2::Source{Variable}) = (p2,p1)
marry(p1::Tuple,p2::Tuple) = p1 == p2 ? p1 : error("Parents of the expressions are not compatible")
marry(p1::Tuple,p2) = p2 in p1 ? p1 : error("Parents of the expressions are not compatible")
marry(p1,p2::Tuple) = p1 in p2 ? p2 : error("Parents of the expressions are not compatible")
marry(p1,p2) = p1 == p2 ? p1 : error("Parents of the expressions are not compatible")


parent(e::Expression) = e.parent
func(e::Expression) = e.func
deriv(e::Term) = e.deriv
deriv(e::Variable) = Dict{Int,Function}(e.index=>con_one)
deriv(e::Parameter) = Dict{Int,Function}()
deriv(e::Real) = Dict{Int,Function}()

con(a) = (x,p=nothing)->a
con_one(x,p=nothing) = 1.
con_zero(x,p=nothing) = 0.
    
fsub(f::Function) = (x,p=nothing)->-f(x,p)
fadd(f1::Function,f2::Function) = (x,p=nothing)->f1(x,p)+f2(x,p)
fmul(f1::Function,f2::Function) = (x,p=nothing)->f1(x,p)*f2(x,p)
fmul(f1::Function,f2::typeof(con_one)) = f1
fmul(f1::typeof(con_one),f2::Function) = f2
fsub(f1::Function,f2::Function) = (x,p=nothing)->f1(x,p)-f2(x,p)
fpow(f1::Function,f2::Function) = (x,p=nothing)->f1(x,p)^f2(x,p)
fcom(f1::Function,f2::Function) = (x,p=nothing)->f1(f2(x,p))
fcom(f1::Function,f2::Function,f3::Function) = (x,p=nothing)->f1(f2(x,p),f3(x,p))

for T in [typeof(con_one),typeof(con_zero),typeof(con(0)),typeof(con(0.))]
    @eval begin
        fcom(f1::$T,f2::Function) = f1
        fcom(f1::$T,f2::Function,f3::Function) = f1
    end
end

for T in Reals
    @eval begin 
        fadd(f1::Function,f2::$T) = f2 == 0 ? f1 : (x,p=nothing)->f1(x,p)+f2
        fadd(f1::$T,f2::Function) = f1 == 0 ? f2 : (x,p=nothing)->f1+f2(x,p)
        fsub(f1::Function,f2::$T) = f2 == 0 ? f1 : (x,p=nothing)->f1(x,p)-f2
        fsub(f1::$T,f2::Function) = f1 == 0 ? fsub(f2) : (x,p=nothing)->f1-f2(x,p)
        fmul(f1::Function,f2::$T) = f2 == 1 ? f1 : f2 == 0 ? con_zero : (x,p=nothing)->f1(x,p)*f2
        fmul(f1::$T,f2::Function) = f1 == 1 ? f2 : f1 == 0 ? con_zero : (x,p=nothing)->f1*f2(x,p)
        fpow(f1::Function,f2::$T) = f2 == 1 ? f1 : f2 == 0 ? con_one  : (x,p=nothing)->f1(x,p)^f2
        fpow(f1::$T,f2::Function) = f1 == 1 ? f1 : f1 == 0 ? con_zero : (x,p=nothing)->f1^f2(x,p)
        fcom(f1::Function,f2::Function,a::$T) = (x,p=nothing)->f1(f2(x,p),a)
        fcom(f1::Function,a::$T,f3::Function) = (x,p=nothing)->f1(a,f3(x,p))
    end
end
fsum(fs) = (x,p=nothing)->sum(f(x,p) for f in fs)

function f_add_sum(f1,f2)
    if hasfield(typeof(f1),:fs)
        i = findfirst(f->f2 isa eltype(f.fs),f1.fs)
        if i==nothing            
            push!(f1.fs,fsum([f2]))
        else
            push!(f1.fs[i].fs,f2)
        end
        return f1
    else
        if typeof(f1)==typeof(f2)
            return fsum([fsum([f1,f2])])
        else
            return fsum([
                fsum([f1]),
                fsum([f2])
            ])
        end
    end
end

function deriv!(deriv,f)
    for (i,d) in deriv
        deriv[i] = f(d)
    end
    return deriv
end
function derivadd!(deriv1,deriv2)
    for (i,d) in deriv2
        deriv1[i] = haskey(deriv1,i) ? fadd(deriv1[i],d) : d
    end
    return deriv1
end
function derivsub!(deriv1,deriv2)
    for (i,d) in deriv2
        deriv1[i] = haskey(deriv1,i) ? fsub(deriv1[i],d) : fsub(d)
    end
    return deriv1
end
function derivmul!(deriv,f)
    for (i,d) in deriv
        deriv[i] = fmul(d,f)
    end
    return deriv
end
derivmul!(deriv1,deriv2,f1,f2) = derivadd!(derivmul!(deriv1,f2),derivmul!(deriv2,f1))

+(e::Expression) = e
-(e::Expression) = Term(parent(e),fsub(func(e)),deriv!(deriv(e),fsub))
+(e1::Expression,e2::Expression) = Term(parent(e1,e2),fadd(func(e1),func(e2)),derivadd!(deriv(e1),deriv(e2)))
-(e1::Expression,e2::Expression) = Term(parent(e1,e2),fsub(func(e1),func(e2)),derivsub!(deriv(e1),deriv(e2)))
*(e1::Expression,e2::Expression) = Term(parent(e1,e2),fmul(func(e1),func(e2)),derivmul!(deriv(e1),deriv(e2),func(e1),func(e2)))

for T in Reals
    @eval begin
        +(e::Expression,a::$T) = Term(parent(e),fadd(func(e),a),deriv(e))
        +(a::$T,e::Expression) = Term(parent(e),fadd(a,func(e)),deriv(e))
        -(e::Expression,a::$T) = Term(parent(e),fsub(func(e),a),deriv(e))
        -(a::$T,e::Expression) = Term(parent(e),fsub(a,func(e)),deriv!(deriv(e),d->fsub(d)))
        *(e::Expression,a::$T) = Term(parent(e),fmul(func(e),a),derivmul!(deriv(e),a))
        *(a::$T,e::Expression) = Term(parent(e),fmul(a,func(e)),derivmul!(deriv(e),a))
        ^(e::Expression,a::$T) = Term(parent(e),fpow(func(e),a),deriv!(deriv(e),d->fmul(a,fmul(fpow(func(e),a-1),d))))
        ^(a::$T,e::Expression) = Term(parent(e),fpow(a,func(e)),deriv!(deriv(e),d->fmul(fmul(fpow(a,fsub(func(e),1)),d),log(a))))
    end
end
for (T1,T2) in [(Expression,Expression),[(Expression,T) for T in Reals]...,[(T,Expression) for T in Reals]...]
    @eval begin
        /(e1::$T1,e2::$T2) = e1*inv(e2)
    end
end

add_sum(e1::Expression,e2::Expression) = Term(parent(e1,e2),f_add_sum(func(e1),func(e2)),derivadd!(deriv(e1),deriv(e2)))


for (M,f,nargs) in diffrules
    if nargs == 1
        df = :((x,p=nothing)->$(diffrule(M,f,:x)))
        @eval begin
            $f(e::Expression)=Term(
                parent(e),
                fcom($f,func(e)),
                deriv!(deriv(e),d->fmul(fcom($df,func(e)),d))
            )
        end
    elseif nargs == 2
        dfsym,~ = diffrule(M,f,:x,:a)
        df = :((x,a,p=nothing)->$dfsym)
        for T in Reals
            @eval begin
                $f(e::Expression,a::$T) = Term(
                    parent(e),
                    fcom($f,func(e),a),
                    deriv!(deriv(e),d->fmul(fcom($df,func(e),a),d))
                )
            end
        end

        ~,dfsym = diffrule(M,f,:a,:x)
        df = :((a,x,p=nothing)->$dfsym)

        for T in Reals
            @eval begin
                $f(a::$T,e::Expression) = Term(
                    parent(e),
                    fcom($f,a,func(e)),
                    deriv!(deriv(e),d->fmul(fcom($df,a,func(e)),d))
                )
            end
        end

        dfsym1,dfsym2 = diffrule(M,f,:x1,:x2)
        df1 = :((x1,x2,p=nothing)->$dfsym1)
        df2 = :((x1,x2,p=nothing)->$dfsym2)

        @eval begin
            $f(e1::Expression,e2::Expression) = Term(
                parent(e1,e2),
                fcom($f,func(e1),func(e2)),
                derivadd!(
                    deriv!(deriv(e1),d->fmul(fcom($df1,func(e1),func(e2)),d)),
                    deriv!(deriv(e2),d->fmul(fcom($df2,func(e1),func(e2)),d))
                )
            )
        end
    end
end

