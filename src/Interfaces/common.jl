sum_init_0(itr) = isempty(itr) ? 0 : sum(itr)
    

@inline function fill!(f,pairs,x,p)
    @simd for l in eachindex(pairs)
        (i,k,d) = pairs[l]
        @inbounds f[i] += d(x,p)
    end
end

@inline function fill_jacobian!(J,pairs,x,p)
    @simd for l in eachindex(pairs)
        (i,k,d) = pairs[l]
        @inbounds J[k] = d(x,p)
    end
end

@inline function fill_hessian!(H,pairs,x,p,sig,lag)
    for ((i,j),k,d) in pairs
        @inbounds H[k] = d(x,p,sig,lag)
    end
end


get_pairss(dict) =  begin
    k = 0
    [
        [(i,k+=1,d) for (i,d) in dict if typeof(d)==type]
        for type in union(typeof(d) for (i,d) in dict)
    ]
end
get_con_pairss(cons) = [
    [(i,nothing,func(cons[i])) for i=1:length(cons) if typeof(func(cons[i]))==type]
    for type in union(typeof(func(con)) for con in cons)
]

function eval_grad!(obj,p)
    pairss = get_pairss(deriv(obj))
    @inline function (f,x)
        f.=0
        for pairs in pairss
            fill!(f,pairs,x,p)
        end
    end
end

function eval_con!(cons,p)
    pairss = get_con_pairss(cons)
    @inline function (c,x)
        c.= 0
        for pairs in pairss
            fill!(c,pairs,x,p)
        end
    end
end

get_jac_dict(cons)= [((i,j),d) for i=1:length(cons) for (j,d) in deriv(cons[i])]

function eval_jac!(cons,p)
    jac_dict = get_jac_dict(cons)
    pairss = get_pairss(jac_dict)
    
    @inline function jac!(J,x)
        for pairs in pairss
            fill_jacobian!(J,pairs,x,p)
        end
    end
    @inline function jac_sparsity!(I,J)
        offset = 0
        for pairs in pairss
            sparsity!(I,J,pairs)
        end
    end
    nnz_jac = sum_init_0(length(pairs) for pairs in pairss)
    
    (jac!,jac_sparsity!,nnz_jac)
end

function get_hess_dict(obj_2nd,con_2nds)
    hess_dict = Dict{Tuple{Int,Int},Function}()
    
    for (i,dobj) in obj_2nd
        for (j,d) in dobj
            # if haskey(hess_dict,(i,j))
            #     dold = hess_dict[i,j]
            #     hess_dict[i,j] = (x,p,sig,lag)->dold(x,p,sig,lag)+d(x,p)*sig
            # else
            i>=j && (hess_dict[i,j] = (x,p,sig,lag)->d(x,p)*sig)
            # end
        end
    end

    for k=1:length(con_2nds)
        for (i,dcon) in con_2nds[k]
            for (j,d) in dcon
                if haskey(hess_dict,(i,j))
                    dold = hess_dict[i,j]
                    hess_dict[i,j] = (x,p,sig,lag)->dold(x,p,sig,lag)+d(x,p)*lag[k]
                else
                    i>=j && (hess_dict[i,j] = (x,p,sig,lag)->d(x,p)*lag[k])
                end
            end
        end
    end
    return hess_dict
end

get_obj_2nd(obj) = [(i,deriv(d(Variable(),Parameter()))) for (i,d) in deriv(obj)]
get_con_2nds(cons) = [[(i,deriv(d(Variable(),Parameter()))) for (i,d) in deriv(con)] for con in cons]

# TODO: need to think about multi-threading in model creation

# function get_obj_2nd(obj)
#     ids = [(i,d) for (i,d) in deriv(obj)]
#     result = Vector{Tuple{Int,Dict{Int,Function}}}(undef,length(ids))
#     Threads.@threads for j=1:length(ids)
#         (i,d) = ids[j]
#         result[j] = (i, deriv(d(Variable(),Parameter())))
#     end
#     return result
# end
# function get_con_2nds(cons)
#     result = Vector{Vector{Tuple{Int,Dict{Int,Function}}}}(undef,length(cons))
#     Threads.@threads for j=1:length(cons)
#         [(i,deriv(d(Variable(),Parameter()))) for (i,d) in deriv(cons[j])]
#     end
#     return result
# end

function eval_hess!(obj,cons,p)
    obj_2nd = get_obj_2nd(obj)
    con_2nds = get_con_2nds(cons)
    
    hess_dict = get_hess_dict(obj_2nd,con_2nds)
    pairss = get_pairss(hess_dict)
    
    @inline function hess!(H,x,lag,sig)
        for pairs in pairss
            fill_hessian!(H,pairs,x,p,sig,lag)
        end
    end

    @inline function hess_sparsity!(I,J)
        for pairs in pairss
            sparsity!(I,J,pairs)
        end
    end

    nnz_hess = sum_init_0(length(pairs) for pairs in pairss) 

    (hess!,hess_sparsity!,nnz_hess)
end

function sparsity!(I,J,pairs)
    for ((i,j),k,d) in pairs
        I[k] = i
        J[k] = j
    end
end


function get_nlp_functions(objs,cons,p)
    obj = sum_init_0(obj for obj in objs)
    
    _obj = x->func(obj)(x,p)
    _grad! = eval_grad!(obj,p)
    _con! = eval_con!(cons,p)
    _jac!,_jac_sparsity!,nnz_jac = eval_jac!(cons,p)
    _hess!,_hess_sparsity!,nnz_hess = eval_hess!(obj,cons,p)
    
    (_obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess)
end
