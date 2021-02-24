struct Expression
    fun::Function
    der::Dict{Int,Function}
end

const Reals = [Int,Float64]
const derivatives = Dict(
    :exp => :exp,
    :sin => :cos,
    :cos => :dcos,
    :log => :inv,
    :inv => :dinv,
    :tan => :dtan
)

dcos(x) = -sin(x)
dinv(x) = -1/x^2
dtan(x) = (cos(x))^-2

const derivatives_2_1 = Dict(
    :^ => :d1pow,
)
d1pow(x,a) = a*x^(a-1)

const derivatives_2_2 = Dict(
    :^ => :d2pow
)
d2pow(a,x) = a^x*log(a)

fentry(n) = x->x[n]
fcon(a) = x->a
fminus(f) = x->-f(x)
fadd(f1::Function,f2::Function) = x->f1(x)+f2(x)
fmul(f1::Function,f2::Function) = x->f1(x)*f2(x)
fsub(f1::Function,f2::Function) = x->f1(x)-f2(x)
fdiv(f1::Function,f2::Function) = x->f1(x)/f2(x)
fcom(f1::Function,f2::Function) = x->f2(f1(x))
for T in Reals
    @eval begin 
        fadd(f1::Function,f2::$T) = x->f1(x)+f2
        fadd(f1::$T,f2::Function) = x->f1+f2(x)
        fmul(f1::Function,f2::$T) = x->f1(x)*f2
        fmul(f1::$T,f2::Function) = x->f1*f2(x)
    end
end
fcom_2_1(f1,f2,a) = x->f2(f1(x),a)
fcom_2_2(f1,f2,a) = x->f2(a,f1(x))
fsum(fs) = x->sum(f(x) for f in fs)


function f_add_sum(f1,f2)
    if hasfield(typeof(f1),:fs)
        i = findfirst(x->f2 isa eltype(x.fs),f1.fs)
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

const fone  = (x)-> 1.
const fzero = (x)-> 0.

Expression() = Expression(identity,Dict{Int,Function}())

fun(e::Expression) = e.fun
der(e::Expression) = e.der
der(e::Real) = Dict{Int,Function}()
nz(e::Expression) = keys(e.der)
nz(e1,e2) = union(nz(e1),nz(e2))
nz(es::Vector{Expression}) = union(nz.(es)...)
nz(es::Expression...) = union(nz.(es)...)

getindex(e::Expression,n) = Expression(fentry(n),Dict{Int,Function}(n=>fone))

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
function dermul!(der,f)
    for (i,d) in der
        der[i] = fmul(d,f)
    end
    return der
end
function dermul!(der1,der2,f1,f2)
    dermul!(der1,f2)
    dermul!(der2,f1)
    return deradd!(der1,der2)
end

+(e::Expression) = e
-(e::Expression) = Expression(fminus(fun(e)),der!(der(e),fminus))
+(e1::Expression,e2::Expression) = Expression(fadd(fun(e1),fun(e2)),deradd!(der(e1),der(e2)))
*(e1::Expression,e2::Expression) = Expression(fmul(fun(e1),fun(e2)),dermul!(der(e1),der(e2),fun(e1),fun(e2)))
for T in Reals
    @eval begin
        +(e1::Expression,e2::$T) = Expression(fadd(fun(e1),e2),der(e1))
        +(e1::$T,e2::Expression) = Expression(fadd(e1,fun(e2)),der(e2))
        *(e1::Expression,e2::$T) = Expression(fmul(fun(e1),e2),dermul!(der(e1),e2))
        *(e1::$T,e2::Expression) = Expression(fmul(e1,fun(e2)),dermul!(der(e2),e1))
    end
end
add_sum(e1::Expression,e2::Expression) = Expression(f_add_sum(fun(e1),fun(e2)),deradd!(der(e1),der(e2)))

for T in Reals
    for (T1,T2) in [(Expression,Expression),(Expression,T),(T,Expression)]
        @eval begin
            -(e1::$T1,e2::$T2) = +(e1,-e2)
            /(e1::$T1,e2::$T2) = *(e1,inv(e2))
        end
    end
end

for (f,df) in derivatives
    @eval begin
        $f(e::Expression)=Expression(
            fcom(fun(e),$f),
            der!(der(e),d->fmul(fcom(fun(e),$df),d))
        )
    end
end

for (f,df) in derivatives_2_1
    for T in Reals
        @eval begin
            $f(e::Expression,a::$T) = Expression(
                fcom_2_1(fun(e),$f,a),
                der!(der(e),d->fmul(fcom_2_1(fun(e),$df,a),d))
            )
        end
    end
end

for (f,df) in derivatives_2_2
    for T in Reals
        @eval begin
            $f(a::$T,e::Expression) = Expression(
                fcom_2_2(fun(e),$f,a),
                der!(der(e),d->fmul(fcom_2_2(fun(e),$df,a),d))
            )
        end
    end
end

struct PrintSymbol
    str::String
end

string(p::PrintSymbol) = p.str
PrintSymbol() = PrintSymbol("x")
getindex(p::PrintSymbol,n) = PrintSymbol(string(p)*"[$n]")

-(p::PrintSymbol) = PrintSymbol("(-"*string(p)*")")
for T in Reals
    for (T1,T2) in [(PrintSymbol,T),(T,PrintSymbol),(PrintSymbol,PrintSymbol)]
        @eval begin
            +(p1::$T1,p2::$T2) = PrintSymbol(string(p1)*" + "*string(p2))
            -(p1::$T1,p2::$T2) = PrintSymbol(string(p1)*" - "*string(p2))
            *(p1::$T1,p2::$T2) = PrintSymbol("("*string(p1)*")*("*string(p2)*")")
            /(p1::$T1,p2::$T2) = PrintSymbol("("*string(p1)*")/("*string(p2)*")")
            ^(p1::$T1,p2::$T2) = PrintSymbol("("*string(p1)*")^("*string(p2)*")")
        end
    end
end
for f in keys(derivatives)
    @eval begin
        $f(p::PrintSymbol) = PrintSymbol(string($f)*"( "*string(p)*" )")
    end
end

string(e::Expression) = string(fun(e)(PrintSymbol()))
print(io::IO,e::Expression) = print(io, string(e))
show(io::IO,::MIME"text/plain",e::Expression) = print(io,e)
