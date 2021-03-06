mutable struct Model
    vars::Vector{PrintVariable}
    pars::Vector{PrintVariable}
    
    cons::Vector{Expression}
    objs::Vector{Expression}

    n::Int # num vars
    q::Int # num pars
    m::Int # num cons

    x::Vector{Float64}
    p::Vector{Float64}
    g::Vector{Float64}
    l::Vector{Float64}
    zl::Vector{Float64}
    zu::Vector{Float64}
    xl::Vector{Float64}
    xu::Vector{Float64}
    gl::Vector{Float64}
    gu::Vector{Float64}

    prob
    optimizer
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

constraint(m,eq::Equality) = constraint(m,eq.e)
constraint(m,eq::Inequality) = constraint(m,eq.e;ub=Inf)

parent(c::Constraint) = c.parent
index(c::Constraint) = c.index
func(c::Constraint) = func(parent(c).cons[index(c)])

SimpleModel(optimizer=nothing;opt...) = Model(optimizer;opt...)
Model(optimizer=nothing;opt...) =Model(
    PrintVariable[],PrintVariable[],Expression[],Expression[],0,0,0,Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],nothing,optimizer,
    Dict{Symbol,Any}(),Dict{Symbol,Any}(name=>option for (name,option) in opt))

set_optimizer(m::Model,optimizer) = m.optimizer = optimizer

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

function parameter(m::Model,val=0.;name="$DEFAULT_PAR_STRING[$(m.q+1)]")
    m.q += 1
    push!(m.pars,PrintVariable("$name"))
    push!(m.p,val)
    Parameter(m.q;parent=m)
end

function constraint(m::Model,e::Expression;lb=0.,ub=0.,name=nothing)
    m.m += 1
    push!(m.cons,e)
    push!(m.l,0.)
    push!(m.g,0.)
    push!(m.gl,lb)
    push!(m.gu,ub)
    Constraint(m,m.m)
end

function objective(m::Model,e::Expression)
    push!(m.objs,e)
    nothing
end

instantiate!(m::Model) = instantiate!(m,m.optimizer;m.opt...)
function optimize!(m::Model)
    m.prob == nothing && instantiate!(m)
    optimize!(m.prob,m.optimizer;m.opt...)
end

value(e::Term{Model}) = func(e)(parent(e).x,parent(e).p)
for (T,var) in [(Variable{Model},:x), (Parameter{Model},:p), (Constraint,:g)]
    @eval value(e::$T) = func(e)(parent(e).x,parent(e).p)
    @eval setvalue(e::$T,val) = parent(e).$var[e.index] = val
end

for (T,ub,lb) in [(Variable{Model},:xl,:xu), (Constraint,:gl,:gu)]
    @eval set_lower_bound(e::$T,val) = parent(e).$lb[e.index] = val
    @eval set_upper_bound(e::$T,val) = parent(e).$ub[e.index] = val
    @eval lower_bound(e::$T) = parent(e).$lb[e.index]
    @eval upper_bound(e::$T) = parent(e).$ub[e.index]
end

dual(c::Constraint) = parent(c).l[index(c)]
objective_value(m::Model) = sum_init_0(value(obj) for obj in m.objs)

num_variables(m::Model) = length(m.vars)
num_constraints(m::Model) = length(m.cons)
