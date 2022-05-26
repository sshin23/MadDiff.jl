struct NullNLPModelMeta{T,V <: Vector{T}} <: NLPModels.AbstractNLPModelMeta{T, V} end
struct NullNLPCore{T} <: MadDiffCore.AbstractNLPCore{T} end 
getproperty(::NullNLPModelMeta,::Symbol) = error("The model is not instantiated yet.")

"""
    MadDiffModel{T <: Real}

A mathematical model of a nonlinaer program.
"""
mutable struct MadDiffModel{T <: AbstractFloat} <: NLPModels.AbstractNLPModel{T,Vector{T}}
    con::MadDiffCore.Sink
    obj::MadDiffCore.Expression

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

    threaded::Bool
    instantiated::Bool
    meta::NLPModels.AbstractNLPModelMeta{T, Vector{T}}
    nlpcore::MadDiffCore.AbstractNLPCore{T}
    counters::NLPModels.Counters
    
    ext::Dict{Symbol,Any}
end

"""
    MadDiffModel{T}()

Creates an empty `MadDiffModel{T}`.

**Example**
`m = MadDiffModel{Float32}()`
"""
MadDiffModel{T}(;threaded = false) where T = MadDiffModel{T}(
    MadDiffCore.Field{T}(), MadDiffCore.Constant(zero(T)),
    0, 0, 0,
    T[], T[], T[], T[], T[], T[], T[], T[], T[], 
    threaded, false, NullNLPModelMeta{T,Vector{T}}(), NullNLPCore{T}(), NLPModels.Counters(),
    Dict{Symbol,Any}()
)

"""
    MadDiffModel()

Creates an empty `MadDiffModel{Float64}`.

**Example**
`m = MadDiffModel(linear_solver = "ma27")`
"""
MadDiffModel(;threaded = false)= MadDiffModel{Float64}(; threaded = threaded)


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
A model component (eitehr a variable or a parameter) of MadDiffModel.
"""
struct ModelComponent{C <: Union{MadDiffCore.Variable,MadDiffCore.Parameter}}
    parent::MadDiffModel
    inner::C
end

for (f,df,ddf) in MadDiffCore.f_nargs_1
    @eval $f(e1::T) where T <: ModelComponent = $f(e1.inner) 
end
for (f,df1,df2,ddf1,ddf12,ddf22) in MadDiffCore.f_nargs_2
    @eval $f(e1::T1,e2::T2) where {T1 <: ModelComponent,T2 <: ModelComponent} = $f(e1.inner,e2.inner)
    @eval $f(e1::T,e2) where T <: ModelComponent = $f(e1.inner,e2)
    @eval $f(e1,e2::T) where T <: ModelComponent = $f(e1,e2.inner)
end
add_sum(a::E,b::ModelComponent{C}) where {C,E<:MadDiffCore.ExpressionSum} = add_sum(a,b.inner)
add_sum(a::ModelComponent{C},b::E) where {C,E<:MadDiffCore.ExpressionSum} = add_sum(a.inner,b)
add_sum(a::ModelComponent{C},b::ModelComponent{C}) where C = add_sum(a.inner,b.inner)
convert(::Type{MadDiffCore.Expression},e::ModelComponent{C}) where C = e

getindex(m::MadDiffModel,idx::Symbol) = m.ext[idx]
setindex!(m::MadDiffModel,val,idx::Symbol) = setindex!(m.ext,val,idx)

"""
    variable(m::MadDiffModel{T}; lb=-Inf, ub=Inf, start=0.)
Creates a variable for MadDiffModel.

**Example**
```
m = MadDiffModel()

x = variable(m; lb = -1, ub = 1, start = 0.5)
```
"""
function variable(m::MadDiffModel{T}; lb=-Inf, ub=Inf, start=0.) where T
    m.instantiated = false
    m.n += 1
    push!(m.x,start)
    push!(m.xl,lb)
    push!(m.xu,ub)
    push!(m.zl,1.)
    push!(m.zu,1.)
    ModelComponent(m,MadDiffCore.Variable{T}(m.n))
end

"""
    parameter(m::MadDiffModel{T}, val)

Creates a parameter for MadDiffModel with value `val`.
**Example**
```
m = MadDiffModel()

p = parameter(m, 0.5)
```
"""
function parameter(m::MadDiffModel{T}, val) where T
    m.instantiated = false
    m.q += 1
    push!(m.p,val)
    ModelComponent(m,MadDiffCore.Parameter{T}(m.q))
end

"""
    objective(m::MadDiffModel, e::MadDiffCore.Expression
Sets the objective function for MadDiffModel. Only minimization is supported.
**Example**
```
m = MadDiffModel()

x = [variable(m) for i=1:3]
objective(m, x[1]^2 + x[2] + sin(x[3]))
```
"""
function objective(m::MadDiffModel, e::MadDiffCore.Expression)
    m.instantiated = false
    m.obj = e
    return 
end

"""
    constraint(m::MadDiffModel, e::MadDiffCore.Expression; lb=0., ub=0.)
Adds a constraint to MadDiffModel.
**Example**
```
m = MadDiffModel()

x = [variable(m) for i=1:3]
constraint(m, x[1]^2 + 2*sin(x[2]) - exp(x[3]) >= 0)
constraint(m, x[1]^4+ x[2]^4 x[3]^4; lb = 0.1, ub = 1.)
"""
function constraint(m::MadDiffModel, e::MadDiffCore.Expression; lb=0., ub=0.)
    m.instantiated = false
    m.m += 1
    m.con[m.m] = e
    push!(m.l,0.)
    push!(m.gl,lb)
    push!(m.gu,ub)
    Constraint(m,m.m)
end
constraint(m::MadDiffModel,eq::MadDiffCore.Expression2{T,typeof(==),E1,E2}) where {T <: AbstractFloat, E1<:MadDiffCore.Expression, E2<:MadDiffCore.Expression} =
    constraint(m,eq.e1-eq.e2)
constraint(m::MadDiffModel,eq::MadDiffCore.Expression2{T,typeof(<=),E1,E2}) where {T <: AbstractFloat, E1<:MadDiffCore.Expression, E2<:MadDiffCore.Expression} =
    constraint(m,eq.e1-eq.e2;ub=Inf)
constraint(m::MadDiffModel,eq::MadDiffCore.Expression2{T,typeof(>=),E1,E2}) where {T <: AbstractFloat, E1<:MadDiffCore.Expression, E2<:MadDiffCore.Expression} =
    constraint(m,eq.e1-eq.e2;lb=-Inf)
constraint(m::MadDiffModel,eq::MadDiffCore.Expression2{T,typeof(==),E1,E2}) where {T <: AbstractFloat, E1<:MadDiffCore.Expression, E2<:Real} =
    constraint(m,eq.e1; lb=eq.e2, ub=eq.e2)
constraint(m::MadDiffModel,eq::MadDiffCore.Expression2{T,typeof(<=),E1,E2}) where {T <: AbstractFloat, E1<:MadDiffCore.Expression, E2<:Real} =
    constraint(m,eq.e1; lb=-Inf, ub=eq.e2)
constraint(m::MadDiffModel,eq::MadDiffCore.Expression2{T,typeof(>=),E1,E2}) where {T <: AbstractFloat, E1<:MadDiffCore.Expression, E2<:Real} =
    constraint(m,eq.e1; lb=eq.e2, ub=Inf)
constraint(m::MadDiffModel,eq::MadDiffCore.Expression2{T,typeof(==),E1,E2}) where {T <: AbstractFloat, E1<:Real, E2<:MadDiffCore.Expression} =
    constraint(m,eq.e2; lb=eq.e1, ub=eq.e1)
constraint(m::MadDiffModel,eq::MadDiffCore.Expression2{T,typeof(<=),E1,E2}) where {T <: AbstractFloat, E1<:Real, E2<:MadDiffCore.Expression} =
    constraint(m,eq.e2; lb=eq.e1, ub=Inf)
constraint(m::MadDiffModel,eq::MadDiffCore.Expression2{T,typeof(>=),E1,E2}) where {T <: AbstractFloat, E1<:Real, E2<:MadDiffCore.Expression} =
    constraint(m,eq.e2; lb=-Inf, ub=eq.e1)

"""
    instantiate!(m::MadDiffModel; sparse = true)
Instantiates the model `m`. The model *must be* instantiated before solving.
# Example
```julia
using MadDiff, NLPModelsIpopt

m = MadDiffModel() 

x = [variable(m) for i=1:3]
objective(m, x[1]^2 + x[2]^2 + sin(x[3]))
constraint(m, 3x[2]^2 <= 1.)

instantiate!(m)
ipopt(m)
```
"""
function instantiate!(m::MadDiffModel; sparse = true)
    if m.n == 0
        error("The model is empty")
    end

    if sparse 
        m.nlpcore = MadDiffCore.SparseNLPCore(
            m.obj,m.con
        )
        m.meta = NLPModels.NLPModelMeta(
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
    else
        m.nlpcore = MadDiffCore.DenseNLPCore(
            m.obj,m.con
        )
        m.meta = NLPModels.NLPModelMeta(
            m.n,
            x0 = m.x,
            lvar = m.xl,
            uvar = m.xu,
            ncon = m.m,
            y0 = m.l,
            lcon = m.gl,
            ucon = m.gu,
            minimize = true
        )
    end
    
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
    NLPModels.increment!(m, :neval_obj)
    MadDiffCore.obj(m.nlpcore,x,m.p; threaded = m.threaded)
end
@inline function NLPModels.cons!(m::MadDiffModel,x::AbstractVector,y::AbstractVector)
    NLPModels.increment!(m, :neval_cons)
    MadDiffCore.cons!(m.nlpcore,x,y,m.p; threaded = m.threaded)
    return 
end
@inline function NLPModels.grad!(m::MadDiffModel,x::AbstractVector,y::AbstractVector)
    NLPModels.increment!(m, :neval_grad)
    MadDiffCore.grad!(m.nlpcore,x,y,m.p; threaded = m.threaded)
    return 
end
@inline function NLPModels.jac_coord!(m::MadDiffModel,x::AbstractVector,J::AbstractVector)
    NLPModels.increment!(m, :neval_jac)
    MadDiffCore.jac_coord!(m.nlpcore,x,J,m.p; threaded = m.threaded)
    return 
end
@inline function NLPModels.hess_coord!(m::MadDiffModel,x::AbstractVector,lag::AbstractVector,z::AbstractVector; obj_weight = 1.0)
    NLPModels.increment!(m, :neval_hess)
    MadDiffCore.hess_coord!(m.nlpcore,x,lag,z,m.p; obj_weight = obj_weight, threaded = m.threaded)
    return 
end

@inline function fill_sparsity!(I,J,tuples)
    @simd for l in eachindex(tuples)
        (i,j) = tuples[l]
        @inbounds I[l] = i
        @inbounds J[l] = j
    end
end

function show(io::IO, m::MadDiffModel{T}) where T
    if m.instantiated 
        println(io, "MadDiffModel{$T} (instantiated).")
        show(io, m.meta)
        show(io, m.counters)
    else
        println(io, "MadDiffModel{$T} (not instantiated).")
    end
end
