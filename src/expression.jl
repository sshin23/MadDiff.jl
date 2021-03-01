const Reals = [Int,Float64]

abstract type Expression end

struct Variable <: Expression
    fun::Function
    der::Dict{Int,Function}
end

mutable struct Term <: Expression
    fun::Function
    der::Dict{Int,Function}
end

struct PrintSymbol
    str::String
end

Variable() = Variable(con_zero,Dict{Int,Function}())
Variable(n) = Variable(fentry(n),Dict{Int,Function}(n=>con_one))
Term() = Variable(con_zero,Dict{Int,Function}())
Term!(e::Variable,f,d) = Term(f,d)
Term!(e::Term,f,d) = begin
    e.fun = f
    e.der = d
    return e
end
PrintSymbol() = PrintSymbol("x")

fun(e::Expression) = e.fun
der(e::Term) = e.der
der(e::Variable) = copy(e.der)
der(e::Real) = Dict{Int,Function}()
string(p::PrintSymbol) = p.str

con(a) = x->a
con_one(x) = 1.
con_zero(x) = 0.
const TypeCons = [typeof(con_one),typeof(con_zero),typeof(con(0)),typeof(con(0.))]
    
getindex(e::Expression,n) = Variable(n)
getindex(p::PrintSymbol,n) = PrintSymbol(string(p)*"[$n]")

string(e::Expression) = string(fun(e)(PrintSymbol()))
print(io::IO,e::Expression) = print(io, string(e))
show(io::IO,::MIME"text/plain",e::Expression) = print(io,e)

fsub(f::Function) = x->-f(x)
fentry(n) = x->x[n]
fadd(f1::Function,f2::Function) = x->f1(x)+f2(x)
fmul(f1::Function,f2::Function) = x->f1(x)*f2(x)
fmul(f1::Function,f2::typeof(con_one)) = f1
fmul(f1::typeof(con_one),f2::Function) = f2
fsub(f1::Function,f2::Function) = x->f1(x)-f2(x)
fpow(f1::Function,f2::Function) = x->f1(x)^f2(x)
fcom(f1::Function,f2::Function) = x->f1(f2(x))
fcom(f1::Function,f2::Function,f3::Function) = x->f1(f2(x),f3(x))

for T in TypeCons
    @eval begin
        fcom(f1::$T,f2::Function) = f1
        fcom(f1::$T,f2::Function,f3::Function) = f1
    end
end

for T in Reals
    @eval begin 
        fadd(f1::Function,f2::$T) = f2 == 0 ? f1 : x->f1(x)+f2
        fadd(f1::$T,f2::Function) = f1 == 0 ? f2 : x->f1+f2(x)
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

function der!(der,f)
    for (i,d) in der
        der[i] = f(d)
    end
    return der
end
function deradd!(der1,der2)
    for (i,d) in der2
        der1[i] = haskey(der1,i) ? fadd(der1[i],d) : d
    end
    return der1
end
function dersub!(der1,der2)
    for (i,d) in der2
        der1[i] = haskey(der1,i) ? fsub(der1[i],d) : fsub(d)
    end
    return der1
end
function dermul!(der,f)
    for (i,d) in der
        der[i] = fmul(d,f)
    end
    return der
end
dermul!(der1,der2,f1,f2) = deradd!(dermul!(der1,f2),dermul!(der2,f1))

+(e::Expression) = e
-(e::Expression) = Term!(e,fsub(fun(e)),der!(der(e),fsub))
+(e1::Expression,e2::Expression) = Term!(e1,fadd(fun(e1),fun(e2)),deradd!(der(e1),der(e2)))
*(e1::Expression,e2::Expression) = Term!(e1,fmul(fun(e1),fun(e2)),dermul!(der(e1),der(e2),fun(e1),fun(e2)))

for T in Reals
    @eval begin
        +(e::Expression,a::$T) = Term!(e,fadd(fun(e),a),der(e))
        +(a::$T,e::Expression) = Term!(e,fadd(a,fun(e)),der(e))
        *(e::Expression,a::$T) = Term!(e,fmul(fun(e),a),dermul!(der(e),a))
        *(a::$T,e::Expression) = Term!(e,fmul(a,fun(e)),dermul!(der(e),a))
        ^(e::Expression,a::$T) = Term!(e,fpow(fun(e),a),der!(der(e),d->fmul(a,fmul(fpow(fun(e),a-1),d))))
        ^(a::$T,e::Expression) = Term!(e,fpow(a,fun(e)),der!(der(e),d->fmul(fmul(a,fpow(fun(e)-1),d),log(a))))
    end
end
for (T1,T2) in [(Expression,Expression),[(Expression,T) for T in Reals]...,[(T,Expression) for T in Reals]...]
    @eval begin
        -(e1::$T1,e2::$T2) = e1+(-e2)
        /(e1::$T1,e2::$T2) = e1*inv(e2)
    end
end

-(p::PrintSymbol) = PrintSymbol("(-"*string(p)*")")
for (T1,T2) in [(PrintSymbol,PrintSymbol),[(PrintSymbol,T) for T in Reals]...,[(T,PrintSymbol) for T in Reals]...]
    @eval begin
        +(p1::$T1,p2::$T2) = PrintSymbol(string(p1)*" + "*string(p2))
        -(p1::$T1,p2::$T2) = PrintSymbol(string(p1)*" - "*string(p2))
        *(p1::$T1,p2::$T2) = PrintSymbol("("*string(p1)*")*("*string(p2)*")")
        /(p1::$T1,p2::$T2) = PrintSymbol("("*string(p1)*")/("*string(p2)*")")
        ^(p1::$T1,p2::$T2) = PrintSymbol("("*string(p1)*")^("*string(p2)*")")
    end
end

add_sum(e1::Expression,e2::Expression) = Term!(e1,f_add_sum(fun(e1),fun(e2)),deradd!(der(e1),der(e2)))


for (M,f,nargs) in diffrules()
    if M in [ :Base ] && f in [:inv,:abs2,:exp,:exp2,:exp10,:log,:log2,:log10,
                               :sin,:cos,:tan,:csc,:sec,:cot,:asin,:acos,:atan,:acsc,:asec,:acot,
                               :sind,:cosd,:tand,:cscd,:secd,:cotd,:asind,:acosd,:atand,:acscd,:asecd,:acotd,
                               :sinh,:cosh,:tanh,:csch,:sech,:coth,:asinh,:acosh,:atanh,:acsch,:asech,:acoth]
        @eval import $M.$f
        if nargs == 1
            df = :(x->$(diffrule(M,f,:x)))
            @eval begin
                $f(e::Expression)=Term!(
                    e,
                    fcom($f,fun(e)),
                    der!(der(e),d->fmul(fcom($df,fun(e)),d))
                )
                $f(p::PrintSymbol) = PrintSymbol(string($f)*"( "*string(p)*" )")
            end
        elseif nargs == 2
            dfsym,~ = diffrule(M,f,:x,:a)
            df = :((x,a)->$dfsym)
            for T in Reals
                @eval begin
                    $f(e::Expression,a::$T) = Term!(
                        e,
                        fcom($f,fun(e),a),
                        der!(der(e),d->fmul(fcom($df,fun(e),a),d))
                    )
                    $f(p::PrintSymbol,a::$T) = PrintSymbol(string($f)*"( "*string(p)*" , "*string(a)*" )")
                end
            end

            ~,dfsym = diffrule(M,f,:a,:x)
            df = :((a,x)->$dfsym)

            for T in Reals
                @eval begin
                    $f(a::$T,e::Expression) = Term!(
                        e,
                        fcom($f,a,fun(e)),
                        der!(der(e),d->fmul(fcom($df,a,fun(e)),d))
                    )
                    $f(a::$T,p::PrintSymbol) = PrintSymbol(string($f)*"( "*string(a)*" , "*string(pa)*" )")
                end
            end

            dfsym1,dfsym2 = diffrule(M,f,:x1,:x2)
            df1 = :((x1,x2)->$dfsym1)
            df2 = :((x1,x2)->$dfsym2)

            @eval begin
                $f(e1::Expression,e2::Expression) = Term!(
                    e1,
                    fcom($f,fun(e1),fun(e2)),
                    deradd!(
                        der!(der(e1),d->fmul(fcom($df1,fun(e1),fun(e2)),d)),
                        der!(der(e2),d->fmul(fcom($df2,fun(e1),fun(e2)),d))
                    )
                )
                $f(p1::PrintSymbol,p2::PrintSymbol) = PrintSymbol(string($f)*"( "*string(p1)*" , "*string(p2)*" )")
            end
        end
    end
end

