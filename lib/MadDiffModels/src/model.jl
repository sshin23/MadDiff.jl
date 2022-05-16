"""
    MadDiffModel{T <: Real}

A mathematical model of a nonlinaer program of the following form:

```
min  f(x)
s.t. xl <=   x  <= xu
     gl <= g(x) <= gu,
```
where `f:R^n -> R` and `g(·):R^n -> R^m`.
"""
Base.@kwdef mutable struct MadDiffModel{T <: Real} <: AbstractNLPModel{T,Vector{T}}
    con::Sink{Field} = Field()
    obj::Union{Nothing,Expression} = ExpressionNull()

    n::Int = 0 # num vars
    q::Int = 0 # num pars
    m::Int = 0 # num cons

    
    x::Vector{T} = Float64[]
    p::Vector{T} = Float64[]
    l::Vector{T} = Float64[]
    zl::Vector{T} = Float64[]
    zu::Vector{T} = Float64[]
    xl::Vector{T} = Float64[]
    xu::Vector{T} = Float64[]
    gl::Vector{T} = Float64[]
    gu::Vector{T} = Float64[]

    instantiated::Bool = false
    meta::Union{Nothing,NLPModelMeta{T, Vector{T}}} = nothing
    nlpcore::Union{Nothing,NLPCore} = nothing
    counters::Union{Nothing,NLPModelsCounters} = nothing
    
    ext::Dict{Symbol,Any} = Dict{Symbol,Any}()
    opt::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

"""
    Constraint
A constraint index of MadDiffModel.
"""
struct Constraint
    parent::MadDiffModel
    index::Int
end

"""
    ModelComponent
A model component (variable or parameter) of MadDiffModel.
"""
struct ModelComponent{C <: Union{Variable,Parameter}}
    parent::MadDiffModel
    inner::C
end

for (f,df,ddf) in f_nargs_1
    @eval $f(e1::T) where T <: ModelComponent = $f(inner(e1)) 
end
for (f,df1,df2,ddf1,ddf12,ddf22) in f_nargs_2
    @eval $f(e1::T1,e2::T2) where {T1 <: ModelComponent,T2 <: ModelComponent} = $f(inner(e1),e2)
    @eval $f(e1::T,e2) where T <: ModelComponent = $f(inner(e1),e2)
    @eval $f(e1,e2::T) where T <: ModelComponent = $f(e1,inner(e2))
end

add_sum(a::E,b::ModelComponent{C}) where {C,E<:ExpressionSum} = add_sum(a,b.inner)
add_sum(a::ModelComponent{C},b::E) where {C,E<:ExpressionSum} = add_sum(a.inner,b)
add_sum(a::ModelComponent{C},b::ModelComponent{C}) where C = add_sum(a.inner,b.inner)

parent(c) = c.parent
index(c::ModelComponent{C}) where C = index(inner(c))
convert(::Type{Expression},e::ModelComponent{C}) where C = inner(e)

getindex(m::MadDiffModel,idx::Symbol) = m.ext[idx]
setindex!(m::MadDiffModel,val,idx::Symbol) = setindex!(m.ext,val,idx)

"""
    MadDiffModel(; options...)

Creates a `MadDiffModel` with optional solver options.

**Example**
`MadDiffModel(linear_solver = "ma27")`
"""
MadDiffModel(; options...)=
    MadDiffModel{Float64}(opt = Dict{Symbol,Any}(name=>option for (name,option) in options))

function variable(m::MadDiffModel;lb=-Inf,ub=Inf,start=0.,name=nothing)
    m.instantiated = false
    m.n += 1
    push!(m.x,start)
    push!(m.xl,lb)
    push!(m.xu,ub)
    push!(m.zl,1.)
    push!(m.zu,1.)
    ModelComponent(m,Variable(m.n))
end

function parameter(m::MadDiffModel,val=0.;name=nothing)
    m.instantiated = false
    m.q += 1
    push!(m.p,val)
    ModelComponent(m,Parameter(m.q))
end

function constraint(m::MadDiffModel,e::E;lb=0.,ub=0.) where E <: Expression
    m.instantiated = false
    m.m += 1
    m.con[m.m] = e
    push!(m.l,0.)
    push!(m.gl,lb)
    push!(m.gu,ub)
    Constraint(m,m.m)
end

function objective(m::MadDiffModel,e::E) where E <: Expression
    m.instantiated = false
    m.obj = e
    nothing
end

constraint(m::MadDiffModel,eq::Expression2{typeof(==),E1,E2}) where {E1<:Expression, E2<:Expression} =
    constraint(m,eq.e1-eq.e2)
constraint(m::MadDiffModel,eq::Expression2{typeof(<=),E1,E2}) where {E1<:Expression, E2<:Expression} =
    constraint(m,eq.e1-eq.e2;ub=Inf)
constraint(m::MadDiffModel,eq::Expression2{typeof(>=),E1,E2}) where {E1<:Expression, E2<:Expression} =
    constraint(m,eq.e1-eq.e2;lb=-Inf)
constraint(m::MadDiffModel,eq::Expression2{typeof(==),E1,E2}) where {E1<:Expression, E2<:Real} =
    constraint(m,eq.e1; lb=eq.e2, ub=eq.e2)
constraint(m::MadDiffModel,eq::Expression2{typeof(<=),E1,E2}) where {E1<:Expression, E2<:Real} =
    constraint(m,eq.e1; lb=-Inf, ub=eq.e2)
constraint(m::MadDiffModel,eq::Expression2{typeof(>=),E1,E2}) where {E1<:Expression, E2<:Real} =
    constraint(m,eq.e1; lb=eq.e2, ub=Inf)
constraint(m::MadDiffModel,eq::Expression2{typeof(==),E1,E2}) where {E1<:Real, E2<:Expression} =
    constraint(m,eq.e2; lb=eq.e1, ub=eq.e1)
constraint(m::MadDiffModel,eq::Expression2{typeof(<=),E1,E2}) where {E1<:Real, E2<:Expression} =
    constraint(m,eq.e2; lb=eq.e1, ub=Inf)
constraint(m::MadDiffModel,eq::Expression2{typeof(>=),E1,E2}) where {E1<:Real, E2<:Expression} =
    constraint(m,eq.e2; lb=-Inf, ub=eq.e1)

value(e::ModelComponent{C}) where C = non_caching_eval(inner(e),parent(e).x,parent(e).p)
setvalue(e::ModelComponent{V},val) where V <: Variable = parent(e).x[index(e)] = val
setvalue(e::ModelComponent{P},val) where P <: Parameter= parent(e).p[index(e)] = val
set_lower_bound(e::ModelComponent{V},val) where V <: Variable = parent(e).xl[index(e)] = val
set_upper_bound(e::ModelComponent{V},val) where V <: Variable = parent(e).xu[index(e)] = val
lower_bound(e::ModelComponent{V}) where V <: Variable = parent(e).xl[index(e)]
upper_bound(e::ModelComponent{V}) where V <: Variable = parent(e).xu[index(e)]

value(e::Constraint) = parent(e).g[index(e)]
# value(e::Constraint) = non_caching_eval(parent(e).con[e.index],parent(e).x,parent(e).p)
# setvalue(e::Constraint,val) = parent(e).g[index(e)] = val
set_lower_bound(e::Constraint,val) = parent(e).gl[index(e)] = val
set_upper_bound(e::Constraint,val) = parent(e).gu[index(e)] = val
lower_bound(e::Constraint) = parent(e).gl[index(e)]
upper_bound(e::Constraint) = parent(e).gu[index(e)]


dual(c::Constraint) = parent(c).l[index(c)]
objective_value(m::MadDiffModel) = non_caching_eval(m.obj,m.x,m.p)

# num_variables(m::MadDiffModel) = m.n
# num_constraints(m::MadDiffModel) = m.m

function instantiate!(m::MadDiffModel)
    
    # grad = Gradient(m.obj)
    # jac,jac_sparsity = SparseJacobian(m.con)
    # hess,hess_sparsity = SparseLagrangianHessian(m.obj,grad,inner(m.con),jac)
    
    m.nlpcore = NLPCore(
        m.obj,m.con
    )
    m.meta = NLPModelMeta(
        m.n,
        x0 = m.x,
        lvar = m.xl,
        uvar = m.xu,
        ncon = m.m,
        y0 = m.l,
        lcon = m.gl,
        ucon = m.gu,
        nnzj = length(m.nlpcore.jac_sparsity),
        nnzh = length(m.nlpcore.hess_sparsity),
        minimize = true
    )
    m.counters = NLPModelsCounters()
    
    m.instantiated=true
    
    return m
end

@inline function NLPModels.jac_structure!(m::MadDiffModel,I::AbstractVector{T},J::AbstractVector{T}) where T
    fill_sparsity!(I,J,m.nlpcore.jac_sparsity)
    return 
end
@inline function NLPModels.hess_structure!(m::MadDiffModel,I::AbstractVector{T},J::AbstractVector{T}) where T
    fill_sparsity!(I,J,m.nlpcore.hess_sparsity)
    return 
end

@inline function NLPModels.obj(m::MadDiffModel,x::AbstractVector)
    increment!(m, :neval_obj)
    obj(m.nlpcore,x,m.p)
end
@inline function NLPModels.cons!(m::MadDiffModel,x::AbstractVector,y::AbstractVector)
    increment!(m, :neval_cons)
    cons!(m.nlpcore,x,y,m.p)
    return 
end
@inline function NLPModels.grad!(m::MadDiffModel,x::AbstractVector,y::AbstractVector)
    increment!(m, :neval_grad)
    grad!(m.nlpcore,x,y,m.p)
    return 
end
@inline function NLPModels.jac_coord!(m::MadDiffModel,x::AbstractVector,J::AbstractVector)
    increment!(m, :neval_jac)
    jac_coord!(m.nlpcore,x,J,m.p)
    return 
end
@inline function NLPModels.hess_coord!(m::MadDiffModel,x::AbstractVector,lag::AbstractVector,z::AbstractVector; obj_weight = 1.0)
    increment!(m, :neval_hess)
    hess_coord!(m.nlpcore,x,lag,z,m.p; obj_weight = obj_weight)
    return 
end

@inline function fill_sparsity!(I,J,tuples)
    @simd for l in eachindex(tuples)
        (i,j) = tuples[l]
        @inbounds I[l] = i
        @inbounds J[l] = j
    end
end

function show(io::IO, m::MadDiffModel)
    if m.instantiated 
        println(io, "A MadDiffModel (instantiated).")
        show(io, m.meta)
        show(io, m.counters)
    else
        println(io, "A MadDiffModel (not instantiated).")
    end
end

show(io::IO, m::ModelComponent{C}) where C = show(io,inner(m))

