mutable struct Model
    cons::Vector{Expression}
    obj::Expression

    n::Int
    m::Int

    x::Vector{Float64}
    xl::Vector{Float64}
    xu::Vector{Float64}
    
    gl::Vector{Float64}
    gu::Vector{Float64}

    prob
    optimizer::Module
    ext::Dict{Symbol,Any}
    opt::Dict{Symbol,Any}
end

Model(optimizer::Module;opt...) = Model(
    Expression[],Term(),0,0,Float64[],Float64[],Float64[],Float64[],Float64[],nothing,optimizer,
    Dict{Symbol,Any}(),Dict{Symbol,Any}(name=>option for (name,option) in opt))
Model() = Model(IpoptOptimizer)

variable(m::Model;lb=-Inf,ub=Inf,start=0.,name=nothing) = begin
    m.n += 1
    push!(m.x,start)
    push!(m.xl,lb)
    push!(m.xu,ub)
    Variable(m.n;name=name)
end
constraint(m::Model,e::Expression;lb=0.,ub=0.,name=nothing) = begin
    m.m += 1
    push!(m.cons,e)
    push!(m.gl,lb)
    push!(m.gu,ub)
    Variable(m.m;name=name)
end
objective(m::Model,e::Expression) = begin
    m.obj = e
    nothing
end

function instantiate!(m::Model)
    m.prob = m.optimizer.createProblem(m)
end

function optimize!(m::Model)
    m.prob == nothing && instantiate!(m)
    m.optimizer.solveProblem(m.prob)
end
