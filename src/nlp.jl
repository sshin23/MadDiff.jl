mutable struct Model
    vars::Vector{PrintVariable}
    pars::Vector{PrintVariable}
    
    cons::Vector{Expression}
    obj::Expression

    n::Int # num vars
    q::Int # num pars
    m::Int # num cons

    x::Vector{Float64}
    p::Vector{Float64}
    l::Vector{Float64}
    zl::Vector{Float64}
    zu::Vector{Float64}
    xl::Vector{Float64}
    xu::Vector{Float64}
    gl::Vector{Float64}
    gu::Vector{Float64}

    prob
    optimizer::Module
    ext::Dict{Symbol,Any}
    opt::Dict{Symbol,Any}
end

struct Constraint
    parent::Model
    index::Int
end

struct Equality e::Expression end
struct Inequality e::Expression end

==(e1::Expression,e2) = Equality(e1-e2)
==(e1,e2::Expression) = Equality(e1-e2)
>=(e1::Expression,e2) = Inequality(e1-e2)
>=(e1,e2::Expression) = Inequality(e1-e2)
<=(e1::Expression,e2) = Inequality(e2-e1)
<=(e1,e2::Expression) = Inequality(e2-e1)

constraint(eq::Equality) = constraint(eq.e)
constraint(eq::Inequality) = constraint(eq.e;ub=Inf)

parent(c::Constraint) = c.parent
index(c::Constraint) = c.index
func(c::Constraint) = func(parent(c).cons[index(c)])

string(c::Constraint) = begin
    gl = parent(c).gl[index(c)]
    gu = parent(c).gu[index(c)]
    str = string(parent(c).cons[index(c)])
    return gl==gu ? str * " == $gl" : (gl > -Inf ? "$gl <= " : "") * str * (gl < Inf ? " <= $gu" : "")
end
print(io::IO,e::Constraint) = print(io, string(e))
show(io::IO,::MIME"text/plain",e::Constraint) = print(io,e)


string(m::Model) = "NLP model with $(m.n) variables and $(m.m) constraints"
print(io::IO,e::Model) = print(io, string(e))
show(io::IO,::MIME"text/plain",e::Model) = print(io,e)


Model(optimizer::Module;opt...) = Model(
    PrintVariable[],PrintVariable[],Expression[],Term(),0,0,0,Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],nothing,optimizer,
    Dict{Symbol,Any}(),Dict{Symbol,Any}(name=>option for (name,option) in opt))
Model() = Model(IpoptOptimizer)

PrintSource(m::Model) = (m.vars,m.pars)

function variable(m::Model;lb=-Inf,ub=Inf,start=0.,name="$DEFAULT_VAR_STRING[$(m.n+1)]")
    m.n += 1
    push!(m.vars,PrintVariable("$name"))
    push!(m.x,start)
    push!(m.xl,lb)
    push!(m.xu,ub)
    push!(m.zl,1.)
    push!(m.zu,1.)
    Variable(m.n;parent=m)
end

function parameter(m::Model;value=0.,name="$DEFAULT_PAR_STRING[$(m.q+1)]")
    m.q += 1
    push!(m.pars,PrintVariable("$name"))
    push!(m.p,value)
    Parameter(m.q;parent=m)
end

function constraint(m::Model,e::Expression;lb=0.,ub=0.,name=nothing)
    m.m += 1
    push!(m.cons,e)
    push!(m.l,0.)
    push!(m.gl,lb)
    push!(m.gu,ub)
    Constraint(m,m.m)
end

function objective(m::Model,e::Expression)
    m.obj = e
    nothing
end

function instantiate!(m::Model)
    m.prob = m.optimizer.create_problem(m;m.opt...)
    m
end

function optimize!(m::Model)
    m.prob == nothing && instantiate!(m)
    m.optimizer.solve_problem(m.prob;m.opt...)
end

value(e::Term{Model}) = func(e)(parent(e).x,parent(e).p)
for (T,var) in [(Variable{Model},:x), (Parameter{Model},:p), (Constraint,:g)]
    @eval value(e::$T) = func(e)(parent(e).x,parent(e).p)
    @eval setvalue(e::$T,val) = parent(e).$var[e.n] = val
end

for (T,ub,lb) in [(Variable{Model},:xl,:xu), (Constraint,:gl,:gu)]
    @eval set_lower_bound(e::$T,val) = parent(e).$lb[e.n] = val
    @eval set_upper_bound(e::$T,val) = parent(e).$ub[e.n] = val
end

dual(c::Constraint) = parent(c).l[index(c)]
getobjectivevalue(m::Model) = value(m.obj)

