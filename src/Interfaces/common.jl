sum_init_0(itr) = isempty(itr) ? 0 : sum(itr)

const VAR = Variable()
const PAR = Parameter()
struct HessianSource
    index::Int
end
getindex(e::HessianSource,n) = n <= e.index ? Variable(n;parent=VAR) : Term(VAR,fxentry(n),ImmutableDict{Int,Function}())

@inline function add!(obj,funcs,x,p)
    @simd for func in funcs
        obj[] += func(x,p)
    end
end

@inline function fill!(f,pairs,x,p)
    @simd for l in eachindex(pairs)
        i,~,d = pairs[l]
        @inbounds f[i] += d(x,p)
    end
end

@inline function fill_jacobian!(J,pairs,x,p)
    @simd for l in eachindex(pairs)
        ~,k,d = pairs[l]
        @inbounds J[k] = d(x,p)
    end
end

@inline function fill_hessian!(H,pairs,x,p,sig,lag)
    @simd for l in eachindex(pairs)
        ~,k,d = pairs[l]
        @inbounds H[k] = d(x,p,sig,lag)
    end
end

@inline function sparsity!(I,J,pairs)
    @simd for l in eachindex(pairs)
        (i,j),k,~ = pairs[l]
        @inbounds I[k] = i
        @inbounds J[k] = j
    end
end

function eval_obj(objs,p)
    funcss = get_obj_funcss(objs)
    @inline function (x)
        obj = [.0]
        @simd for l in eachindex(funcss)
            @inbounds add!(obj,funcss[l],x,p)
        end
        return obj[]
    end
end

function eval_grad!(objs,p)
    pairss = get_pairss(get_grad_dict(objs))
    @inline function (f,x)
        f.=0
        @simd for l in eachindex(pairss)
            @inbounds fill!(f,pairss[l],x,p)
        end
    end
end

function eval_con!(cons,p)
    pairss = get_con_pairss(cons)
    @inline function (c,x)
        c.= 0
        @simd for l in eachindex(pairss)
            @inbounds fill!(c,pairss[l],x,p)
        end
    end
end

function eval_jac!(cons,p)
    pairss = get_pairss(get_jac_dict(cons))
    
    @inline function jac!(J,x)
        @simd for l in eachindex(pairss)
            @inbounds fill_jacobian!(J,pairss[l],x,p)
        end
    end
    @inline function jac_sparsity!(I,J)
        offset = 0
        @simd for l in eachindex(pairss)
            @inbounds sparsity!(I,J,pairss[l])
        end
    end
    nnz_jac = sum_init_0(length(pairss[l]) for l in eachindex(pairss))
    
    (jac!,jac_sparsity!,nnz_jac)
end

function eval_hess!(objs,cons,p)
    pairss = get_pairss(get_hess_dict(get_2nds(objs),get_2nds(cons)))

    @inline function hess!(H,x,lag,sig)
        @simd for l in eachindex(pairss)
            @inbounds fill_hessian!(H,pairss[l],x,p,sig,lag)
        end
    end

    @inline function hess_sparsity!(I,J)
        @simd for l in eachindex(pairss)
            @inbounds sparsity!(I,J,pairss[l])
        end
    end

    nnz_hess = sum_init_0(length(pairss[l]) for l in eachindex(pairss)) 

    (hess!,hess_sparsity!,nnz_hess)
end

get_pairss(dict) = (k=0; [[(i,k+=1,d) for (i,d) in dict if typeof(d)==type] for type in union(typeof(d) for (i,d) in dict)])
get_con_pairss(cons) = [[(i,nothing,func(cons[i])) for i=1:length(cons) if typeof(func(cons[i]))==type] for type in union(typeof(func(con)) for con in cons)]
get_obj_funcss(objs) = [[func(objs[i]) for i in eachindex(objs) if typeof(func(objs[i]))==type] for type in union(typeof(func(obj)) for obj in objs)]
get_2nds(exprs) = [[(i,deriv(d(HessianSource(i),PAR))) for (i,d) in deriv(expr)] for expr in exprs]
get_grad_dict(objs) = [(i,d) for obj in objs for (i,d) in deriv(obj)]
get_jac_dict(cons)= [((i,j),d) for i=1:length(cons) for (j,d) in deriv(cons[i])]

function get_hess_dict(obj_2nds,con_2nds)
    hess_dict = Dict{Tuple{Int,Int},Function}()

    for obj_2nd in obj_2nds
        for (i,dobj) in obj_2nd
            hess_dict_obj!(hess_dict,dobj,i)
        end
    end

    for k=1:length(con_2nds)
        for (i,dcon) in con_2nds[k]
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
    _obj = eval_obj(objs,p)
    _grad! = eval_grad!(objs,p)
    _con! = eval_con!(cons,p)
    _jac!,_jac_sparsity!,nnz_jac = eval_jac!(cons,p)
    _hess!,_hess_sparsity!,nnz_hess = eval_hess!(objs,cons,p)
    
    (_obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess)
end


