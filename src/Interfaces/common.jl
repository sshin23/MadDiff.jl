
const VAR = Variable()
const PAR = Parameter()
struct HessianSource
    index::Int
end

getindex(e::HessianSource,n) = n <= e.index ? Variable(n;parent=VAR) : Term(VAR,fxentry(n),ImmutableDict{Int,Function}())
sum_init_0(itr) = isempty(itr) ? 0 : sum(itr)

@inline function add!(obj,funcs,x,p)
    @simd for func in funcs
        obj[] += func(x,p)
    end
end

@inline function fill!(f,tuples,x,p)
    @simd for l in eachindex(tuples)
        i,~,d = tuples[l]
        @inbounds f[i] += d(x,p)
    end
end

@inline function fill_jacobian!(J,tuples,x,p)
    @simd for l in eachindex(tuples)
        ~,k,d = tuples[l]
        @inbounds J[k] = d(x,p)
    end
end

@inline function fill_hessian!(H,tuples,x,p,sig,lag)
    @simd for l in eachindex(tuples)
        ~,k,d = tuples[l]
        @inbounds H[k] = d(x,p,sig,lag)
    end
end

@inline function sparsity!(I,J,tuples)
    @simd for l in eachindex(tuples)
        (i,j),k,~ = tuples[l]
        @inbounds I[k] = i
        @inbounds J[k] = j
    end
end

function get_evaluators_obj(objs,p)
    funcss = get_obj_funcss(objs)
    obj = [.0]
    @inline function (x)
        obj[] = 0.
        @simd for l in eachindex(funcss)
            @inbounds add!(obj,funcss[l],x,p)
        end
        return obj[]
    end
end

function get_evaluators_grad!(objs,p)
    tupless = get_grad_tupless(objs)
    @inline function (f,x)
        f.=0
        @simd for l in eachindex(tupless)
            @inbounds fill!(f,tupless[l],x,p)
        end
    end
end

function get_evaluators_con!(cons,p)
    tupless = get_con_tupless(cons)
    @inline function (c,x)
        c.= 0
        @simd for l in eachindex(tupless)
            @inbounds fill!(c,tupless[l],x,p)
        end
    end
end

function get_evaluators_jac!(cons,p)
    tupless = get_jac_tupless(cons)
    
    @inline function jac!(J,x)
        @simd for l in eachindex(tupless)
            @inbounds fill_jacobian!(J,tupless[l],x,p)
        end
    end
    @inline function jac_sparsity!(I,J)
        @simd for l in eachindex(tupless)
            @inbounds sparsity!(I,J,tupless[l])
        end
    end
    nnz_jac = sum_init_0(length(tupless[l]) for l in eachindex(tupless))
    
    return jac!,jac_sparsity!,nnz_jac
end

function get_evaluators_hess!(grad!,jac!,p)
    tupless = get_hess_tupless(get_hess_dict(get_obj_2nds(grad!.tupless),get_con_2nds(jac!.tupless)))

    @inline function hess!(H,x,lag,sig)
        @simd for l in eachindex(tupless)
            @inbounds fill_hessian!(H,tupless[l],x,p,sig,lag)
        end
    end

    @inline function hess_sparsity!(I,J)
        @simd for l in eachindex(tupless)
            @inbounds sparsity!(I,J,tupless[l])
        end
    end

    nnz_hess = sum_init_0(length(tupless[l]) for l in eachindex(tupless)) 

    return hess!,hess_sparsity!,nnz_hess
end

get_obj_2nds(tupless) = [[(i,deriv(d(HessianSource(i),Parameter()))) for (i,j,d) in tuples] for tuples in tupless]
get_con_2nds(tupless) = [[(k,i,deriv(d(HessianSource(i),Parameter()))) for ((k,i),j,d) in tuples] for tuples in tupless]

push_dict!(d,dict) = haskey(dict,typeof(d)) ? push!(dict[typeof(d)],d) : (dict[typeof(d)] = [d])
push_dict!(d,dict,i) = haskey(dict,typeof(d)) ? push!(dict[typeof(d)],(i,nothing,d)) : (dict[typeof(d)] = [(i,nothing,d)])
push_dict!(d,dict,i,j,k) = haskey(dict,typeof(d)) ? push!(dict[typeof(d)],((i,j),k[]+=1,d)) : (dict[typeof(d)] = [((i,j),k[]+=1,d)])

function get_obj_funcss(objs)
    dict = Dict{DataType,Vector{F} where F <: Function}()
    @simd for i in eachindex(objs)
        push_dict!(func(objs[i]),dict)
    end
    return collect(values(dict))
end
function get_grad_tupless(objs)
    dict = Dict{DataType,Vector{Tuple{Int,Nothing,F}} where F <: Function}()
    @simd for i in eachindex(objs)
        get_grad_tupless(objs[i],dict,i)
    end
    return collect(values(dict))
end
function get_grad_tupless(obj,dict,i)
    for (j,d) in deriv(obj)
        push_dict!(d,dict,j)
    end
end
function get_con_tupless(cons)
    dict = Dict{DataType,Vector{Tuple{Int,Nothing,F}} where F <: Function}()
    @simd for i in eachindex(cons)
        push_dict!(func(cons[i]),dict,i)
    end
    return collect(values(dict))
end
function get_jac_tupless(cons)
    dict = Dict{DataType,Vector{Tuple{Tuple{Int,Int},Int,F}} where F <: Function}()
    k = [0]
    @simd for i in eachindex(cons)
        get_jac_tupless(cons[i],dict,i,k)
    end
    return collect(values(dict))
end
function get_jac_tupless(con,dict,i,k)
    for (j,d) in deriv(con)
        push_dict!(d,dict,i,j,k)
    end
end
function get_hess_tupless(hess)
    dict = Dict{DataType,Vector{Tuple{Tuple{Int,Int},Int,F}} where F <: Function}()
    k = [0]
    for ((i,j),d) in hess
        push_dict!(d,dict,i,j,k)
    end
    return collect(values(dict))
end


function get_hess_dict(obj_2nds,con_2nds)
    hess_dict = Dict{Tuple{Int,Int},Function}()

    for obj_2nd in obj_2nds
        for (i,dobj) in obj_2nd
            hess_dict_obj!(hess_dict,dobj,i)
        end
    end

    for con_2nd in con_2nds
        for (k,i,dcon) in con_2nd
            hess_dict_con!(hess_dict,dcon,i,k)
        end
    end
    
    return hess_dict
end
function hess_dict_obj!(hess_dict,dobj,i)
    for (j,d) in dobj
        if haskey(hess_dict,(i,j))
            hess_dict[i,j] = f_hess_obj_entry(d,hess_dict[i,j])
        else
            hess_dict[i,j] = f_hess_obj_entry(d)
        end
    end
end
function hess_dict_con!(hess_dict,dcon,i,k)
    for (j,d) in dcon
        if haskey(hess_dict,(i,j))
            hess_dict[i,j] = f_hess_con_entry(d,hess_dict[i,j],k)
        else
            hess_dict[i,j] = f_hess_con_entry(d,k)
        end
    end
end
f_hess_obj_entry(d) = @inline (x,p,sig,lag)->d(x,p)*sig
f_hess_obj_entry(d,dold) = @inline (x,p,sig,lag)->dold(x,p,sig,lag)+d(x,p)*sig
f_hess_con_entry(d,k) = @inline (x,p,sig,lag)->d(x,p)*lag[k]
f_hess_con_entry(d,dold,k) = @inline (x,p,sig,lag)->dold(x,p,sig,lag)+d(x,p)*lag[k]

function get_nlp_functions(objs,cons,p)
    _obj = get_evaluators_obj(objs,p)
    _grad! = get_evaluators_grad!(objs,p)
    _con! = get_evaluators_con!(cons,p)
    _jac!,_jac_sparsity!,nnz_jac = get_evaluators_jac!(cons,p)
    _hess!,_hess_sparsity!,nnz_hess = get_evaluators_hess!(_grad!,_jac!,p)
    
    return _obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess
end


