mutable struct Model
    vars::Vector{PrintVariable}
    cons::Vector{Expression}
    obj::Expression

    n::Int
    m::Int

    x::Vector{Float64}
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
    PrintVariable[],Expression[],Term(),0,0,Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],nothing,optimizer,
    Dict{Symbol,Any}(),Dict{Symbol,Any}(name=>option for (name,option) in opt))
Model() = Model(IpoptOptimizer)

PrintSource(m::Model) = m.vars

variable(m::Model;lb=-Inf,ub=Inf,start=0.,name="x[$(m.n+1)]") = begin
    m.n += 1
    push!(m.vars,PrintVariable("$name"))
    push!(m.x,start)
    push!(m.xl,lb)
    push!(m.xu,ub)
    push!(m.zl,1.)
    push!(m.zu,1.)
    Variable(m.n;parent=m)
end

constraint(m::Model,e::Expression;lb=0.,ub=0.,name=nothing) = begin
    m.m += 1
    push!(m.cons,e)
    push!(m.l,0.)
    push!(m.gl,lb)
    push!(m.gu,ub)
    Constraint(m,m.m)
end
objective(m::Model,e::Expression) = begin
    m.obj = e
    nothing
end

function instantiate!(m::Model)
    m.prob = m.optimizer.createProblem(m)
    m
end

function optimize!(m::Model)
    m.prob == nothing && instantiate!(m)
    m.optimizer.solveProblem(m.prob)
end

for T in [Variable{Model}, Term{Model}, Constraint]
    @eval value(e::$T) = func(e)(parent(e).x)
end
dual(c::Constraint) = parent(c).l[index(c)]
getobjectivevalue(m::Model) = value(m.obj)
