abstract type Expression end

struct Source
    str::String
end

struct Variable <: Expression
    func::Function
    n::Int
end

mutable struct Term <: Expression
    func::Function
    deriv::Dict{Int,Function}
end

Source() = Source("")
Variable(n;name=nothing) = Variable(name ==nothing ? x->x[n] : x->x isa Number ? x[n] : PrintVariable(name),n)
Term() = Term(con_zero,Dict{Int,Function}())
Term!(e::Variable,f,d) = Term(f,d)
Term!(e::Term,f,d) = begin
    e.func = f
    e.deriv = d
    return e
end

string(e::Variable) = string(func(e)(PrintSource()))
print(io::IO,e::Variable) = print(io, string(e))

string(e::Term) = raw(func(e)(PrintSource()))
print(io::IO,e::Term) = print(io, string(e))

string(e::Source) = e.str
print(io::IO,e::Source) = print(io, string(e))

show(io::IO,::MIME"text/plain",e::Expression) = print(io,e)
show(io::IO,::MIME"text/plain",e::Source) = print(io,e)

func(e::Expression) = e.func
deriv(e::Term) = e.deriv
deriv(e::Variable) = Dict{Int,Function}(e.n=>con_one)
deriv(e::Real) = Dict{Int,Function}()

con(a) = x->a
con_one(x) = 1.
con_zero(x) = 0.
    
getindex(e::Source,n) = Variable(n;name = e.str == "" ? nothing : e.str*"[$n]")

fsub(f::Function) = x->-f(x)
fadd(f1::Function,f2::Function) = x->f1(x)+f2(x)
fmul(f1::Function,f2::Function) = x->f1(x)*f2(x)
fmul(f1::Function,f2::typeof(con_one)) = f1
fmul(f1::typeof(con_one),f2::Function) = f2
fsub(f1::Function,f2::Function) = x->f1(x)-f2(x)
fpow(f1::Function,f2::Function) = x->f1(x)^f2(x)
fcom(f1::Function,f2::Function) = x->f1(f2(x))
fcom(f1::Function,f2::Function,f3::Function) = x->f1(f2(x),f3(x))

for T in [typeof(con_one),typeof(con_zero),typeof(con(0)),typeof(con(0.))]
    @eval begin
        fcom(f1::$T,f2::Function) = f1
        fcom(f1::$T,f2::Function,f3::Function) = f1
    end
end

for T in Reals
    @eval begin 
        fadd(f1::Function,f2::$T) = f2 == 0 ? f1 : x->f1(x)+f2
        fadd(f1::$T,f2::Function) = f1 == 0 ? f2 : x->f1+f2(x)
        fsub(f1::Function,f2::$T) = f2 == 0 ? f1 : x->f1(x)-f2
        fsub(f1::$T,f2::Function) = f1 == 0 ? -f2 : x->f1-f2(x)
        fmul(f1::Function,f2::$T) = f2 == 1 ? f1 : f2 == 0 ? con_zero : x->f1(x)*f2
        fmul(f1::$T,f2::Function) = f1 == 1 ? f2 : f1 == 0 ? con_zero : x->f1*f2(x)
        fpow(f1::Function,f2::$T) = f2 == 1 ? f1 : f2 == 0 ? con_one  : x->f1(x)^f2
        fpow(f1::$T,f2::Function) = f1 == 1 ? f1 : f1 == 0 ? con_zero : x->f1^f2(x)
        fcom(f1::Function,f2::Function,a::$T) = x->f1(f2(x),a)
        fcom(f1::Function,a::$T,f3::Function) = x->f1(a,f3(x))
    end
end
fsum(fs) = x->sum(f(x) for f in fs)

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
-(e::Expression) = Term!(e,fsub(func(e)),deriv!(deriv(e),fsub))
+(e1::Expression,e2::Expression) = Term!(e1,fadd(func(e1),func(e2)),derivadd!(deriv(e1),deriv(e2)))
-(e1::Expression,e2::Expression) = Term!(e1,fsub(func(e1),func(e2)),derivsub!(deriv(e1),deriv(e2)))
*(e1::Expression,e2::Expression) = Term!(e1,fmul(func(e1),func(e2)),derivmul!(deriv(e1),deriv(e2),func(e1),func(e2)))

for T in Reals
    @eval begin
        +(e::Expression,a::$T) = Term!(e,fadd(func(e),a),deriv(e))
        +(a::$T,e::Expression) = Term!(e,fadd(a,func(e)),deriv(e))
        -(e::Expression,a::$T) = Term!(e,fsub(func(e),a),deriv(e))
        -(a::$T,e::Expression) = Term!(e,fsub(a,func(e)),deriv!(e,d->fsub(d)))
        *(e::Expression,a::$T) = Term!(e,fmul(func(e),a),derivmul!(deriv(e),a))
        *(a::$T,e::Expression) = Term!(e,fmul(a,func(e)),derivmul!(deriv(e),a))
        ^(e::Expression,a::$T) = Term!(e,fpow(func(e),a),deriv!(deriv(e),d->fmul(a,fmul(fpow(func(e),a-1),d))))
        ^(a::$T,e::Expression) = Term!(e,fpow(a,func(e)),deriv!(deriv(e),d->fmul(fmul(a,fpow(func(e)-1),d),log(a))))
    end
end
for (T1,T2) in [(Expression,Expression),[(Expression,T) for T in Reals]...,[(T,Expression) for T in Reals]...]
    @eval begin
        -(e1::$T1,e2::$T2) = e1+(-e2)
        /(e1::$T1,e2::$T2) = e1*inv(e2)
    end
end

add_sum(e1::Expression,e2::Expression) = Term!(e1,f_add_sum(func(e1),func(e2)),derivadd!(deriv(e1),deriv(e2)))


for (M,f,nargs) in diffrules
    if nargs == 1
        df = :(x->$(diffrule(M,f,:x)))
        @eval begin
            $f(e::Expression)=Term!(
                e,
                fcom($f,func(e)),
                deriv!(deriv(e),d->fmul(fcom($df,func(e)),d))
            )
        end
    elseif nargs == 2
        dfsym,~ = diffrule(M,f,:x,:a)
        df = :((x,a)->$dfsym)
        for T in Reals
            @eval begin
                $f(e::Expression,a::$T) = Term!(
                    e,
                    fcom($f,func(e),a),
                    deriv!(deriv(e),d->fmul(fcom($df,func(e),a),d))
                )
            end
        end

        ~,dfsym = diffrule(M,f,:a,:x)
        df = :((a,x)->$dfsym)

        for T in Reals
            @eval begin
                $f(a::$T,e::Expression) = Term!(
                    e,
                    fcom($f,a,func(e)),
                    deriv!(deriv(e),d->fmul(fcom($df,a,func(e)),d))
                )
            end
        end

        dfsym1,dfsym2 = diffrule(M,f,:x1,:x2)
        df1 = :((x1,x2)->$dfsym1)
        df2 = :((x1,x2)->$dfsym2)

        @eval begin
            $f(e1::Expression,e2::Expression) = Term!(
                e1,
                fcom($f,func(e1),func(e2)),
                derivadd!(
                    deriv!(deriv(e1),d->fmul(fcom($df1,func(e1),func(e2)),d)),
                    deriv!(deriv(e2),d->fmul(fcom($df2,func(e1),func(e2)),d))
                )
            )
        end
    end
end

