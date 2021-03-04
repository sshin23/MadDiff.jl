sum_init_0(itr) = isempty(itr) ? 0 : sum(itr)
    

function fill!(f,pairs,x,p)
    for (i,d) in pairs
        f[i] += d(x,p)
    end
end

function fill!(J,pairs,x,p,offset)
    for (j,d) in pairs
        J[offset += 1] = d(x,p)
    end
    return offset
end

function fill_hessian!(H,pairs,x,p,sig,lag,offset)
    for ((i,j),d) in pairs
        H[offset += 1] = d(x,p,sig,lag)
    end
    return offset
end


get_pairss(dict) =  [
    [(i,d) for (i,d) in dict if typeof(d)==type]
    for type in union(typeof(d) for (i,d) in dict)
]
get_con_pairss(cons) = [
    [(i,func(cons[i])) for i=1:length(cons) if typeof(func(cons[i]))==type]
    for type in union(typeof(func(con)) for con in cons)
]

function eval_grad!(obj,p)
    pairss = get_pairss(deriv(obj))
    function (f,x)
        f.=0
        for pairs in pairss
            fill!(f,pairs,x,p)
        end
    end
end

function eval_con!(cons,p)
    pairss = get_con_pairss(cons)
    function (c,x)
        c.= 0
        for pairs in pairss
            fill!(c,pairs,x,p)
        end
    end
end

function get_jac_dict(cons)
    jac_dict = Dict{Tuple{Int,Int},Function}()
    for i=1:length(cons)
        for (j,d) in deriv(cons[i])
            # if haskey(jac_dict,(i,j))
            #     dold = jac_dict[i,j]
            #     jac_dict[i,j] = (x,p)->dold(x,p)+d(x,p)
            # else
            jac_dict[i,j] = (x,p)->d(x,p)
            # end
        end
    end
    return jac_dict
end

function eval_jac!(cons,p)
    jac_dict = get_jac_dict(cons)
    pairss = get_pairss(jac_dict)
    
    function jac!(J,x)
        offset = 0
        for pairs in pairss
            offset = fill!(J,pairs,x,p,offset)
        end
    end
    function jac_sparsity!(I,J)
        offset = 0
        for pairs in pairss
            offset = sparsity!(I,J,pairs,offset)
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
    
    function hess!(H,x,lag,sig)
        offset = 0
        for pairs in pairss
            offset = fill_hessian!(H,pairs,x,p,sig,lag,offset)
        end
    end

    function hess_sparsity!(I,J)
        offset = 0
        for pairs in pairss
            offset = sparsity!(I,J,pairs,offset)
        end
    end

    nnz_hess = sum_init_0(length(pairs) for pairs in pairss) 

    (hess!,hess_sparsity!,nnz_hess)
end

function sparsity!(I,J,pairs,offset)
    for ((i,j),d) in pairs
        offset += 1
        I[offset] = i
        J[offset] = j
    end
    return offset
end


function get_nlp_functions(obj,cons,p)

    _obj = x->func(obj)(x,p)
    _grad! = eval_grad!(obj,p)
    _con! = eval_con!(cons,p)
    _jac!,_jac_sparsity!,nnz_jac = eval_jac!(cons,p)
    _hess!,_hess_sparsity!,nnz_hess = eval_hess!(obj,cons,p)
    
    (_obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess)
end
