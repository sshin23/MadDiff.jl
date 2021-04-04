struct Dummy end

const DUMMY = Dummy()

mutable struct Model{Optimizer}
    con::Sink{Field}
    obj::Union{Nothing,Expression}

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
    ext::Dict{Symbol,Any}
    opt::Dict{Symbol,Any}
end

struct ModelComponent{C <: Union{Variable,Parameter}}
    parent::Model
    inner::C
end

struct Constraint
    parent::Model
    index::Int
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

parent(c) = c.parent
index(c::ModelComponent{C}) where C = index(inner(c))
convert(::Type{Expression},e::ModelComponent{C}) where C = inner(e)

getindex(m::Model,idx::Symbol) = m.ext[idx]
setindex!(m::Model,val,idx::Symbol) = setindex!(m.ext,val,idx)

default_optimizer() = @isdefined(DEFAULT_OPTIMIZER) ? DEFAULT_OPTIMIZER : error("DEFAULT_OPTIMIZER is not defined. To use Ipopt as a default optimizer, do: using Ipopt")
Model(Optimizer=default_optimizer();opt...) =Model{Optimizer}(
    Field(),Constant(0.),0,0,0,Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],nothing,Dict{Symbol,Any}(),Dict{Symbol,Any}(name=>option for (name,option) in opt))

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
    push!(m.g,0.)
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

function objective!(m::Model,e::E) where E <: Expression
    m.obj = add_sum(m.obj,e)
    nothing
end

instantiate!(m::Model{Optimizer}) where Optimizer = (m.prob = Optimizer(m))
function optimize!(m::Model)
    m.prob == nothing && instantiate!(m)
    optimize!(m.prob)
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
setvalue(e::Constraint,val) = parent(e).g[index(e)] = val
set_lower_bound(e::Constraint,val) = parent(e).gl[index(e)] = val
set_upper_bound(e::Constraint,val) = parent(e).gu[index(e)] = val
lower_bound(e::Constraint) = parent(e).gl[index(e)]
upper_bound(e::Constraint) = parent(e).gu[index(e)]


dual(c::Constraint) = parent(c).l[index(c)]
objective_value(m::Model) = non_caching_eval(m.obj,m.x,m.p)

num_variables(m::Model) = m.n
num_constraints(m::Model) = m.m


function nlp_evaluator(obj,con)
    
    grad = Gradient(obj)
    jac,jac_sparsity = SparseJacobian(con)
    hess,hess_sparsity = SparseLagrangianHessian(obj,grad,inner(con),jac)
    
    _obj = @inline (x,p)->non_caching_eval(obj,x,p)
    _con = @inline (y,x,p)->non_caching_eval(con,y,x,p)
    _grad! = @inline function (y,x,p)
        y .= 0
        obj(x,p)
        non_caching_eval(grad,y,x,p)
    end
    _jac! = @inline function (J,x,p)
        J .= 0
        con(DUMMY,x,p)
        non_caching_eval(jac,J,x,p)
    end
    _hess! = @inline function (z,x,p,lag,sig)
        z .= 0
        obj(x,p)
        con(DUMMY,x,p)
        grad(DUMMY,x,p)
        jac(DUMMY,x,p)
        hess(z,x,p,lag,sig)
    end
        
    nnz_jac = length(jac_sparsity)
    nnz_hess = length(hess_sparsity)
    
    return _obj,_grad!,_con,_jac!,jac_sparsity,_hess!,hess_sparsity
end

# etc
getindex(::Dummy,key::Int) = 0.
getindex(::Tuple{Int,Dummy},key::Int) = 0.
setindex!(::Dummy,val,key) = nothing
setindex!(::Tuple{Int,Dummy},val,key) = nothing
ndims(::Type{Dummy}) = 0
ndims(::Type{T}) where T <: HessianSum = 0
size(::Dummy) = 0
copyto!(::Dummy,val) = nothing

@inline function fill_sparsity!(I,J,tuples)
    @simd for l in eachindex(tuples)
        (i,j) = tuples[l]
        @inbounds I[l] = i
        @inbounds J[l] = j
    end
end

mutable struct KKTErrorEvaluator{G,C,J}
    f::Vector{Float64}
    c::Vector{Float64}
    
    jac::SparseMatrixCSC{Float64,Int}
    jac_view::SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}
    
    grad!::G
    con!::C
    jac!::J
end

function KKTErrorEvaluator(m::Model)
    
    obj = m.obj
    con = m.con
    p = m.p
    
    grad = Gradient(obj)
    jacobian,jac_sparsity = SparseJacobian(con)
    
    _obj = @inline x->non_caching_eval(obj,x,p)
    _con! = @inline (y,x)->non_caching_eval(con,y,x,p)
    _grad! = @inline function (y,x)
        y .= 0
        obj(x,p)
        non_caching_eval(grad,y,x,p)
    end
    _jac! = @inline function (J,x)
        J .= 0
        con(DUMMY,x,p)
        non_caching_eval(jacobian,J,x,p)
    end
    _jac_fill_sparsity! = @inline (I,J)->fill_sparsity!(I,J,jac_sparsity)
    
    nnz_jac = length(jac_sparsity)

    I = Vector{Int}(undef,nnz_jac)
    J = Vector{Int}(undef,nnz_jac)
    _jac_fill_sparsity!(I,J)

    f = similar(m.x)
    c = similar(m.l)
    
    tmp = sparse(J,I,collect(1:nnz_jac),num_variables(m),num_constraints(m))
    jac = SparseMatrixCSC{Float64,Int64}(tmp.m,tmp.n,tmp.colptr,tmp.rowval,Vector{Float64}(undef,nnz_jac))

    idx = Vector{Int}(undef,length(tmp.nzval))
    for i=eachindex(tmp.nzval)
        idx[tmp.nzval[i]] = i
    end
    jac_view = view(jac.nzval,idx)
    
    return KKTErrorEvaluator(f,c,jac,jac_view,_grad!,_con!,_jac!)
end

function (evaluator::KKTErrorEvaluator{G,C,J})(x,l,gl) where {G,C,J}
    
    evaluator.grad!(evaluator.f,x)
    evaluator.con!(evaluator.c,x)
    evaluator.jac!(evaluator.jac_view,x)

    mul!(evaluator.f,evaluator.jac,l,1.,1.)
    evaluator.c .-= gl

    return max(norm(evaluator.f,Inf),norm(evaluator.c,Inf))
end
