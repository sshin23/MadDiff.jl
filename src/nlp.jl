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
end

Model() = Model(Expression[],Term(),0,0,Float64[],Float64[],Float64[],Float64[],Float64[],nothing)

variable(m::Model;lb=-Inf,ub=Inf,start=0.) = begin
    m.n += 1
    push!(m.x,start)
    push!(m.xl,lb)
    push!(m.xu,ub)
    Variable(m.n)
end
constraint(m::Model,e::Expression;lb=0.,ub=0.) = begin
    m.m += 1
    push!(m.cons,e)
    push!(m.gl,lb)
    push!(m.gu,ub)
    Variable(m.m)
end
objective(m::Model,e::Expression) = begin
    m.obj = e
    nothing
end

export variable,constraint,objective
