struct Equality{E<:Expression} e::E end
struct Inequality{E<:Expression} e::E end
for (T1,T2) in [(Expression,Expression),(Expression,Real),(Real,Expression)]
    @eval begin
        ==(e1::T1,e2::T2) where {T1 <: $T1, T2 <: $T2} = Equality(e1-e2)
        >=(e1::T1,e2::T2) where {T1 <: $T1, T2 <: $T2} = Inequality(e1-e2)
        <=(e1::T1,e2::T2) where {T1 <: $T1, T2 <: $T2} = Inequality(e2-e1)
    end
end


struct Evaluator
    grad
    jac
    hess
    jac_sparsity::Vector{Tuple{Int,Int}}
    hess_sparsity::Vector{Tuple{Int,Int}}
end
mutable struct Model{T} <: AbstractNLPModel{T,Vector{T}}
    con::Sink{Field}
    obj::Union{Nothing,Expression}

    n::Int # num vars
    q::Int # num pars
    m::Int # num cons

    x::Vector{T}
    p::Vector{T}
    l::Vector{T}
    zl::Vector{T}
    zu::Vector{T}
    xl::Vector{T}
    xu::Vector{T}
    gl::Vector{T}
    gu::Vector{T}

    meta::Union{Nothing,NLPModelMeta{T, Vector{T}}}
    evaluator::Union{Nothing,Evaluator}
    counters::Union{Nothing,NLPModelsCounters}
    
    ext::Dict{Symbol,Any}
    opt::Dict{Symbol,Any}
end

struct Constraint
    parent::Model
    index::Int
end

struct ModelComponent{C <: Union{Variable,Parameter}}
    parent::Model
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
for o in [:(==),:(>=),:(<=)]
    @eval $o(e1::T1,e2::T2) where {T1 <: ModelComponent, T2 <: ModelComponent} = $o(inner(e1),inner(e2))
    @eval $o(e1::T1,e2) where {T1 <: ModelComponent} = $o(inner(e1),e2)
    @eval $o(e1,e2::T2) where {T2 <: ModelComponent} = $o(e1,inner(e2))
end
add_sum(a::E,b::ModelComponent{C}) where {C,E<:ExpressionSum} = add_sum(a,b.inner)
add_sum(a::ModelComponent{C},b::E) where {C,E<:ExpressionSum} = add_sum(a.inner,b)
add_sum(a::ModelComponent{C},b::ModelComponent{C}) where C = add_sum(a.inner,b.inner)

parent(c) = c.parent
index(c::ModelComponent{C}) where C = index(inner(c))
convert(::Type{Expression},e::ModelComponent{C}) where C = inner(e)

getindex(m::Model,idx::Symbol) = m.ext[idx]
setindex!(m::Model,val,idx::Symbol) = setindex!(m.ext,val,idx)

Model(;opt...) =Model(
    Field(),Constant(0.),0,0,0,
    Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],
    nothing,nothing,nothing,
    Dict{Symbol,Any}(),Dict{Symbol,Any}(name=>option for (name,option) in opt))

function variable(m::Model;lb=-Inf,ub=Inf,start=0.,name=nothing)
    m.n += 1
    push!(m.x,start)
    push!(m.xl,lb)
    push!(m.xu,ub)
    push!(m.zl,1.)
    push!(m.zu,1.)
    ModelComponent(m,Variable(m.n))
end

function parameter(m::Model,val=0.;name=nothing)
    m.q += 1
    push!(m.p,val)
    ModelComponent(m,Parameter(m.q))
end

function constraint(m::Model,e::E;lb=0.,ub=0.) where E <: Expression
    m.m += 1
    m.con[m.m] = e
    push!(m.l,0.)
    push!(m.gl,lb)
    push!(m.gu,ub)
    Constraint(m,m.m)
end

constraint(m,eq::Equality{E}) where {E<:Expression} = constraint(m,eq.e)
constraint(m,eq::Inequality{E}) where {E<:Expression} = constraint(m,eq.e;ub=Inf)

function objective(m::Model,e::E) where E <: Expression
    m.obj = e
    nothing
end

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
objective_value(m::Model) = non_caching_eval(m.obj,m.x,m.p)

num_variables(m::Model) = m.n
num_constraints(m::Model) = m.m

function instantiate!(m::Model)
    
    grad = Gradient(m.obj)
    jac,jac_sparsity = SparseJacobian(m.con)
    hess,hess_sparsity = SparseLagrangianHessian(m.obj,grad,inner(m.con),jac)
    
    m.evaluator = Evaluator(
        grad,jac,hess,jac_sparsity,hess_sparsity
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
        nnzj = length(jac_sparsity),
        nnzh = length(hess_sparsity),
        minimize = true
    )
    m.counters = NLPModelsCounters()
    
    return 
end
@inline function jac_structure!(m::Model,I::AbstractVector{T},J::AbstractVector{T}) where T
    fill_sparsity!(I,J,m.evaluator.jac_sparsity)
    return 
end
@inline function hess_structure!(m::Model,I::AbstractVector{T},J::AbstractVector{T}) where T
    fill_sparsity!(I,J,m.evaluator.hess_sparsity)
    return 
end
@inline function obj(m::Model,x::AbstractVector)
    increment!(m, :neval_obj)
    non_caching_eval(m.obj,x,m.p)
end
@inline function cons!(m::Model,x::AbstractVector,y::AbstractVector)
    increment!(m, :neval_cons)
    non_caching_eval(m.con,y,x,m.p)
    return 
end
@inline function grad!(m::Model,x::AbstractVector,y::AbstractVector)
    increment!(m, :neval_grad)
    y .= 0
    m.obj(x,m.p)
    non_caching_eval(m.evaluator.grad,y,x,m.p)
    return 
end
@inline function jac_coord!(m::Model,x::AbstractVector,J::AbstractVector)
    increment!(m, :neval_jac)
    J .= 0
    m.con(DUMMY,x,m.p)
    non_caching_eval(m.evaluator.jac,J,x,m.p)
    return 
end
@inline function hess_coord!(m::Model,x::AbstractVector,lag::AbstractVector,z::AbstractVector; obj_weight = 1.0)
    increment!(m, :neval_hess)
    z .= 0
    m.obj(x,m.p)
    m.con(DUMMY,x,m.p)
    m.evaluator.grad(DUMMY,x,m.p)
    m.evaluator.jac(DUMMY,x,m.p)
    m.evaluator.hess(z,x,m.p,lag, obj_weight) ###
    return 
end


# etc
getindex(::Dummy,key::Int) = 0.
getindex(::Tuple{Int,Dummy},key::Int) = 0.
setindex!(::Dummy,val,key) = nothing
setindex!(::Tuple{Int,Dummy},val,key) = nothing
# ndims(::Type{Dummy}) = 0
# ndims(::Type{T}) where T <: HessianSum = 0
# size(::Dummy) = 0
# copyto!(::Dummy,val) = nothing

@inline function fill_sparsity!(I,J,tuples)
    @simd for l in eachindex(tuples)
        (i,j) = tuples[l]
        @inbounds I[l] = i
        @inbounds J[l] = j
    end
end

